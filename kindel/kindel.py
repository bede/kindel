# Author: Bede Constantinides - b|at|bede|dot|im, @beconstant
# License: GPL V3

import io
import os
import sys
import argh
import tqdm
import simplesam
import subprocess
import scipy.stats

from collections import OrderedDict, defaultdict, namedtuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import numpy as np
import pandas as pd


def parse_records(ref_id, ref_len, records):
    '''
    Iterate over records, returning namedtuple of base frequencies, indels and soft clipping info
    weights: lists of dicts of base frequencies after CIGAR reconciliation
    insertions, deletions: list of dicts of base frequencies
    clip_starts, clip_ends: lists of clipping start and end frequencies in an LTR direction
    clip_start_weights, clip_end_weights: base frequencies from left- and right-clipped regions
    '''
    weights = [{'A':0, 'T':0, 'G':0, 'C':0, 'N':0} for p in range(ref_len)]
    clip_start_weights = [{'A':0, 'T':0, 'G':0, 'C':0, 'N':0} for p in range(ref_len)]
    clip_end_weights = [{'A':0, 'T':0, 'G':0, 'C':0, 'N':0} for p in range(ref_len)]
    clip_starts = [0] * (ref_len + 1)
    clip_ends = [0] * (ref_len + 1)
    insertions = [defaultdict(int) for p in range(ref_len + 1)]
    deletions = [0] * (ref_len + 1)
    for record in tqdm.tqdm(records, desc='loading sequences'): # Progress bar
        q_pos = 0 
        r_pos = record.pos-1  # Zero indexed genome coordinates
        for i, cigarette in enumerate(record.cigars):
            length, operation = cigarette
            if operation == 'M':
                for _ in range(length):
                    q_nt = record.seq[q_pos].upper()
                    weights[r_pos][q_nt] += 1
                    r_pos += 1
                    q_pos += 1
            elif operation == 'I':
                nts = record.seq[q_pos:q_pos+length].upper()
                insertions[r_pos][nts] += 1
                q_pos += length
            elif operation == 'D':
                deletions[r_pos] += 1
                r_pos += length
            elif operation == 'S':
                if i == 0:  # Count right-of-gap / l-clipped start positions (e.g. start of ref)
                    clip_ends[r_pos] += 1
                    for gap_i in range(length):
                        q_nt = record.seq[gap_i].upper()
                        rel_r_pos = r_pos - length + gap_i
                        if rel_r_pos >= 0:
                            clip_end_weights[rel_r_pos][q_nt] += 1
                    q_pos += length
                else:  # Count left-of-gap / r-clipped start position (e.g. end of ref)
                    clip_starts[r_pos] += 1 
                    for pos in range(length):
                        q_nt = record.seq[q_pos].upper()
                        if r_pos < ref_len:
                            clip_start_weights[r_pos][q_nt] += 1
                            r_pos += 1
                            q_pos += 1
    clip_start_cov = [sum({nt:w[nt] for nt in list('ACGT')}.values()) for w in clip_start_weights]
    clip_end_cov = [sum({nt:w[nt] for nt in list('ACGT')}.values()) for w in clip_end_weights]
    clip_cov = list(map(lambda x, y: x+y, clip_start_cov, clip_end_cov))
    alignment = namedtuple('alignment', ['ref_id', 'weights', 'insertions', 'deletions',
                           'clip_starts', 'clip_ends', 'clip_start_weights', 'clip_end_weights',
                           'clip_start_cov', 'clip_end_cov', 'clip_cov'])
    return alignment(ref_id, weights, insertions, deletions, clip_starts, clip_ends,
                     clip_start_weights, clip_end_weights, clip_start_cov, clip_end_cov, clip_cov)


def parse_bam(bam_path):
    '''
    Returns alignment information for each reference sequence as an OrderedDict
    '''
    alignments = OrderedDict()
    # print('\n>>>>>>>>>>')
    # print(bam_path, file=sys.stderr)
    with open(bam_path, 'r') as bam_fh:
        bam = simplesam.Reader(bam_fh)
        refs_lens = {n.replace('SN:', ''): int(l[0].replace('LN:', '')) for n, l in bam.header['@SQ'].items()}
        refs_records = {id: (r for r in bam if r.rname == id) for id in refs_lens} # {ref:(records)}
        for ref_id, records in refs_records.items():
            alignments[ref_id] = parse_records(ref_id, refs_lens[ref_id], records)
    return alignments


def clip_dominant_consensuses(weights, clip_start_cov, clip_end_cov, min_depth, terminal_mask=50):
    masked_positions = (list(range(0, terminal_mask)) +
                       list(range(len(ten.weights)-terminal_mask, len(ten.weights))))
    regions_r = [] # read right-clipped / left of low coverage region
    regions_l = []
    for pos, (w, s, s_cov) in enumerate(zip(weights, clip_start_cov)):
        if pos not in masked_positions and s/sum(w.values()) > 0.5:
            regions_r.append(pos)


def clip_start_consensuses(weights, clip_start_weights, clip_start_cov, clip_cov_decay_threshold=0.5):
    # Add min cov requirement??
    positions = list(range(len(weights)))
    masked_positions = positions[:50] + positions[-50:]
    regions = []
    region_active = False
    for pos, sc, w in zip(positions, clip_start_cov, weights):
        sum_weights = sum(w.values())
        if sum_weights and sc/sum_weights > 0.5 and pos not in masked_positions and not region_active:
            region = {'start': pos, 'clip_cns': ''}
            start_pos = pos
            region_active = True
            for pos_, (sc_, sw_, w_) in enumerate(zip(clip_start_cov[pos:],
                                                      clip_start_weights[pos:], weights[pos:])):
                if sc_ > sum(w_.values()) * clip_cov_decay_threshold:
                    region['clip_cns'] += consensus(sw_)[0]
                else:
                    region['end'] =  start_pos + pos_
                    break
            regions.append(region)
    return regions


def clip_end_consensuses(weights, clip_end_weights, clip_end_cov, clip_cov_decay_threshold=0.5):
    positions = list(range(len(weights)))
    masked_positions = positions[:50] + positions[-50:]
    reversed_weights = list(sorted(zip(positions, clip_end_cov, clip_end_weights, weights),
                                   key=lambda x: x[0], reverse=True))
    regions = []
    region_active = False
    for i, (pos, ec, ew, w) in enumerate(reversed_weights):
        if ec/sum(w.values()) > 0.5 and pos not in masked_positions and not region_active:
            region = {'end': pos+1}  # Start with end since we're iterating in reverse
            region_active = True
            rev_clip_cns = None
            for pos_, ec_, ew_, w_ in reversed_weights[len(positions)-pos:]:
                if ec_ > sum(w_.values()) * clip_cov_decay_threshold:
                    if not rev_clip_cns:  # Add first base to account for lag in clip coverage
                        rev_clip_cns = consensus(clip_end_weights[pos_+1])[0]
                    rev_clip_cns += consensus(ew_)[0]
                else:
                    region['start'] =  pos_
                    region['clip_cns'] = rev_clip_cns[::-1]
                    break
            regions.append(region)
    return regions


# find_gaps() > find_concordant_alternatives()
# For destruction
def discordant_sites(weights, clip_starts, clip_ends, min_depth, terminal_mask=50):
    masked_positions = list(range(0, terminal_mask)) + list(range(len(ten.weights)-terminal_mask, len(ten.weights)))
    gaps = []
    gap = namedtuple('gap', ['start', 'end'])
    aligned_depth = [sum({nt:w[nt] for nt in list('ACGT')}.values()) for w in weights]
    gap_open = False
    for i, (aln_dep, clip_s, clip_e) in enumerate(zip(aligned_depth, clip_starts, clip_ends)):
        threshold_freq = int(max(aln_dep*0.5, min_depth))
        if clip_s > threshold_freq and not gap_open and i:
            # print(clip_s, threshold_freq)
            gap_start = i
            gap_open = True
        elif gap_open and clip_e >= threshold_freq and i:
            # print(clip_s, threshold_freq)
            gap_end = i
            gaps.append(gap(gap_start, gap_end))
            gap_open = False
    # print(gaps)
    return gaps


def consensus(weight):
    '''
    Returns tuple of consensus base, weight and flag indicating a tie for consensus
    Removed namedtuples for performance
    '''
    base, frequency = max(weight.items(), key=lambda x:x[1]) if sum(weight.values()) else ('N', 0)
    weight_sans_consensus = {k:d for k, d in weight.items() if k != base}
    tie = True if frequency and frequency in weight_sans_consensus.values() else False
    aligned_depth = sum(weight.values())
    proportion = round(frequency/aligned_depth, 2) if aligned_depth else 0
    return(base, frequency, proportion, tie)


def s_overhang_consensus(clip_start_weights, start_pos, min_depth, max_len=500):
    '''
    Returns consensus sequence (string) of clipped reads at specified position
    start_pos is the first position described by the CIGAR-S
    '''
    consensus_overhang = ''
    for pos in range(start_pos, start_pos+max_len):
        pos_consensus = consensus(clip_start_weights[pos])
        if pos_consensus[1] >= min_depth:
            consensus_overhang += pos_consensus[0]
        else:
            break
    return consensus_overhang


def e_overhang_consensus(clip_end_weights, start_pos, min_depth, max_len=500):
    '''
    Returns consensus sequence (string) of clipped reads at specified position
    end_pos is the last position described by the CIGAR-S
    '''
    rev_consensus_overhang = ''
    for pos in range(start_pos, start_pos-max_len, -1):
        pos_consensus = consensus(clip_end_weights[pos])
        # print(clip_end_weights[pos], pos, pos_consensus[1], min_depth)
        if pos_consensus[1] >= min_depth:
            rev_consensus_overhang += pos_consensus[0]
        else:
            break
    consensus_overhang = rev_consensus_overhang[::-1]
    return consensus_overhang


def s_flanking_seq(start_pos, weights, min_depth, k):
    '''
    Returns consensus sequence (string) flanking LHS of soft-clipped gaps
    '''
    flank_seq = ''
    for pos in range(start_pos-k, start_pos):
        pos_consensus = consensus(weights[pos])
        if pos_consensus[1] >= min_depth:
            flank_seq += pos_consensus[0]
    return flank_seq


def e_flanking_seq(end_pos, weights, min_depth, k):
    '''
    Returns consensus sequence (string) flanking RHS of soft-clipped gaps
    '''
    flank_seq = ''
    for pos in range(end_pos+1, end_pos+k+1):
        pos_consensus = consensus(weights[pos])
        if pos_consensus[1] >= min_depth:
            flank_seq += pos_consensus[0]

    return flank_seq


def lcs(s1, s2):
    '''
    Returns longest common substring
    https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Longest_common_substring
    '''
    m = [[0] * (1 + len(s2)) for i in range(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in range(1, 1 + len(s1)):
        for y in range(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return s1[x_longest - longest: x_longest]


def close_by_lcs(l_seq, r_seq):
    '''
    Returns sequence built from merged overlapping left- and right-clipped consensus sequences
    Potentially very slow; used as last resort
    '''
    _lcs = lcs(l_seq, r_seq)
    l_trim = l_seq[:l_seq.find(_lcs)]
    r_trim = r_seq[r_seq.find(_lcs)+len(_lcs):]
    return l_trim + _lcs + r_trim


def reconcile_gaps(gaps, weights, clip_start_weights, clip_end_weights, min_depth, closure_k):
    '''
    Returns dict of consensus insertions between to gap coordinates from find_gaps()
    Dict keys are gap start positions (left of gap)
    '''
    gap_consensuses = {}
    for gap in gaps:
        s_overhang_seq = s_overhang_consensus(clip_start_weights, gap.start, min_depth)
        e_overhang_seq = e_overhang_consensus(clip_end_weights, gap.end, min_depth)
        s_flank_seq = s_flanking_seq(gap.start, weights, min_depth, closure_k)
        e_flank_seq = e_flanking_seq(gap.end, weights, min_depth, closure_k)
        # print(s_overhang_seq)
        # print(e_overhang_seq)
        # print(s_flank_seq)
        # print(e_flank_seq)
        if e_flank_seq in s_overhang_seq:  # Close gap using right-clipped read consensus
            # print('e_flank_seq in s_overhang_seq')
            i = s_overhang_seq.find(e_flank_seq)  # str.find() returns -1 in absence of match
            gap_consensus = s_overhang_seq[:i]
        elif s_flank_seq in e_overhang_seq:  # Close gap using left-clipped read consensus
            # print('s_flank_seq in e_overhang_seq')
            i = e_overhang_seq.find(s_flank_seq)
            gap_consensus = e_overhang_seq[i:]
        elif len(lcs(s_overhang_seq, e_overhang_seq)) >= closure_k:
            gap_consensus = close_by_lcs(s_overhang_seq, e_overhang_seq)
            # print('len(lcs(s_overhang_seq, e_overhang_seq)) >= closure_k')
        else:
            print('Unable to close gap', file=sys.stderr)
            continue
        gap_consensuses[gap.start] = gap_consensus.lower()
    return gap_consensuses


def consensus_sequence(weights, clip_start_weights, clip_end_weights, insertions, deletions,
                       trim_ends, min_depth, uppercase):
    consensus_seq = ''
    changes = [None] * len(weights)
    for pos, weight in tqdm.tqdm(enumerate(weights), total=len(weights), desc='building consensus'):
        ins_freq = sum(insertions[pos].values()) if insertions[pos] else 0
        del_freq = deletions[pos]
        aligned_depth = sum({nt: weight[nt] for nt in list('ACGT')}.values())
        try:
            aligned_depth_next = sum({nt: weights[pos+1][nt] for nt in list('ACGT')}.values())
        except IndexError:
            aligned_depth_next = 0
        threshold_freq = aligned_depth * 0.5
        indel_threshold_freq = min(threshold_freq, aligned_depth_next * 0.5)
        if del_freq > threshold_freq:
            changes[pos] = 'D'
        elif aligned_depth < min_depth:
            consensus_seq += 'N'
            changes[pos] = 'N'
        else:
            if ins_freq > indel_threshold_freq:
                insertion = consensus(insertions[pos])
                consensus_seq += insertion[0].lower() if not insertion[3] else 'N'
                changes[pos] = 'I'
            pos_consensus = consensus(weight)
            consensus_seq += pos_consensus[0] if not pos_consensus[3] else 'N'
    if trim_ends:
        consensus_seq = consensus_seq.strip('N')
    if uppercase:
        consensus_seq = consensus_seq.upper()
    return consensus_seq, changes


def consensus_seqrecord(consensus, ref_id):
    return SeqRecord(Seq(consensus), id=ref_id + '_cns', description='')


def build_report(weights, changes, gaps, gap_consensuses, bam_path, fix_gaps, trim_ends,
                 min_depth, closure_k, uppercase):
    aligned_depth = [sum({nt:w[nt] for nt in list('ACGT')}.values()) for w in weights]
    ambiguous_sites = []
    insertion_sites = []
    deletion_sites = []
    gaps_fmt = ['-'.join([str(g+1) for g in gap]) for gap in gaps] if gaps else []
    gap_consensuses_fmt = ', '.join(gap_consensuses.values()) if gap_consensuses else ''
    for pos, change in enumerate(changes):
        if change == 'N':
            ambiguous_sites.append(str(pos))
        elif change == 'I':
            insertion_sites.append(str(pos+1))
        elif change == 'D':
            deletion_sites.append(str(pos))
    report = '========================= REPORT ===========================\n'
    report += 'options:\n'
    report += '- bam_path: {}\n'.format(bam_path)
    report += '- fix_gaps: {}\n'.format(fix_gaps)
    report += '- trim_ends: {}\n'.format(trim_ends)
    report += '- uppercase: {}\n'.format(uppercase)
    report += '- min_depth: {}\n'.format(min_depth)
    report += '- closure_k: {}\n'.format(closure_k)
    report += 'min,max observed depth[50:-50]: {},{}\n'.format(min(aligned_depth[50:-50]), max(aligned_depth))
    report += 'ambiguous sites: {}\n'.format(', '.join(ambiguous_sites))
    report += 'insertion sites: {}\n'.format(', '.join(insertion_sites))
    report += 'deletion sites: {}\n'.format(', '.join(deletion_sites))
    report += 'soft-clipped gaps: {}\n'.format(', '.join(gaps_fmt))
    report += 'inserted sequences: {}\n'.format(gap_consensuses_fmt)
    report += '============================================================\n'
    return report


def bam_to_consensus(bam_path, fix_gaps=False, trim_ends=False, min_depth=2,
                     closure_k=7, uppercase=False, reporting=True):
    refs_consensuses = []
    refs_changes = {}
    # print(list(parse_bam(bam_path).keys()))
    for ref_id, aln in parse_bam(bam_path).items():
        if fix_gaps:
            gaps = find_gaps(aln.weights, aln.clip_starts, aln.clip_ends, min_depth)
            gap_consensuses = reconcile_gaps(gaps, aln.weights, aln.clip_start_weights,
                                             aln.clip_end_weights, min_depth, closure_k)
        else:
            gaps, gap_consensuses = None, None
        # print(ref_id, file=sys.stderr)
        consensus, changes = consensus_sequence(aln.weights, aln.clip_start_weights,
                                                aln.clip_end_weights, aln.insertions, aln.deletions,
                                                trim_ends, min_depth, uppercase)
        report = build_report(aln.weights, changes, gaps, gap_consensuses, bam_path, fix_gaps,
                              trim_ends, min_depth, closure_k, uppercase)
    refs_consensuses.append(consensus_seqrecord(consensus, ref_id))
    refs_changes[ref_id] = changes
    result = namedtuple('result', ['consensuses', 'refs_changes', 'report'])
    return result(refs_consensuses, refs_changes, report)


def weights(bam_path: 'path to SAM/BAM file',
            relative: 'output relative nucleotide frequencies'=False,
            no_confidence: 'skip confidence calculation'=False):
    '''
    Returns DataFrame of per-site nucleotide frequencies, depth, consensus, lower bound of the
    99% consensus confidence interval, and Shannon entropy
    '''
    def binomial_ci(count, nobs, alpha=0.01):
        '''Returns lower, upper bounds of the Jeffrey binomial proportion confidence interval'''
        lower_ci, upper_ci = scipy.stats.beta.interval(1-alpha, count+0.5, nobs-count+0.5)
        return lower_ci, upper_ci
    
    refs_alns = parse_bam(bam_path)
    weights_fmt = []
    for ref, aln in refs_alns.items():
        weights_fmt.extend([dict(w, ref=ref, pos=i) for i, w in enumerate(aln.weights, start=1)])
    
    weights_df = pd.DataFrame(weights_fmt, columns=['ref','pos','A','C','G','T','N'])
    weights_df['depth'] = weights_df[['A','C','G','T','N']].sum(axis=1)
    consensus_depths_df = weights_df[['A','C','G','T','N']].max(axis=1)
    weights_df['consensus'] = consensus_depths_df.divide(weights_df.depth)
    
    rel_weights_df = pd.DataFrame()
    for nt in ['A','C','G','T','N']:
        rel_weights_df[[nt]] = weights_df[[nt]].divide(weights_df.depth, axis=0)
    
    weights_df['shannon'] = [scipy.stats.entropy(x) for x in rel_weights_df[['A','C','G','T']].as_matrix()]
    
    if not no_confidence:
        weights_df['lower_ci'] = [binomial_ci(c, t)[0] for c, t, in zip(consensus_depths_df, weights_df['depth'])]
    
    if relative:
        for nt in ['A','C','G','T','N']:
            weights_df[[nt]] = rel_weights_df[[nt]]
    
    return weights_df.round(dict(consensus=3, lower_ci=3, shannon=3))


def features(bam_path: 'path to SAM/BAM file'):
    '''
    Returns DataFrame of relative per-site nucleotide frequencies, insertions, deletions and entropy
    '''
    
    refs_alns = parse_bam(bam_path)
    weights_fmt = []
    for ref, aln in refs_alns.items():
        weights_fmt.extend([dict(w, ref=ref, pos=i) for i, w in enumerate(aln.weights, start=1)])
    for pos, weight in enumerate(weights_fmt):
        weight['i'] = sum(aln.insertions[pos].values())
        weight['d'] = aln.deletions[pos]

    # Think about which columns should sum to 1
    weights_df = pd.DataFrame(weights_fmt, columns=['ref','pos','A','C','G','T','N','i','d'])
    weights_df['depth'] = weights_df[['A','C','G','T','N', 'd']].sum(axis=1)
    consensus_depths_df = weights_df[['A','C','G','T','N']].max(axis=1)
    weights_df['consensus'] = consensus_depths_df.divide(weights_df.depth)
    
    for nt in ['A','C','G','T','N','i','d']:
        weights_df[[nt]] = weights_df[[nt]].divide(weights_df.depth, axis=0)
    
    weights_df['shannon'] = [scipy.stats.entropy(x) for x in weights_df[['A','C','G','T','i','d']].as_matrix()]

    return weights_df.round(3)


def variants(bam_path: 'path to SAM/BAM file',
             abs_threshold: 'absolute frequency (0-∞) threshold above which to call variants'=1,
             rel_threshold: 'relative frequency (0.0-1.0) threshold above which to call variants'=0.01,
             only_variants: 'exclude invariant sites from output'=False,
             absolute: 'report absolute variant frequencies'=False):
    '''
    Returns DataFrame of single nucleotide variant frequencies exceeding specified frequency
    thresholds from an aligned BAM
    '''
    weights_df = weights(bam_path, no_confidence=True)
    weights = weights_df[['A','C','G','T']].to_dict('records')
    variant_sites = []
    
    for i, weight in enumerate(weights):
        depth = sum(weight.values())
        consensus = consensus(weight)
        alt_weight = {nt:w for nt, w in weight.items() if nt != consensus[0]}
        alt_weight_rel = {nt:w/depth for nt, w in alt_weight.items() if depth}
        alt_depths = alt_weight.values()
        max_alt_weight = max(alt_weight, key=alt_weight.get)
        max_alt_depth = max(alt_depths)
        alts_above_thresholds = {nt:w for nt, w in alt_weight.items()  # Get weights >= abs & rel
                                 if depth and w >= abs_threshold and w/depth >= rel_threshold}
        if absolute:
            variant_sites.append(alts_above_thresholds)
        else:
            variant_sites.append({nt:round(w/depth, 3) for nt, w in alts_above_thresholds.items() if depth})

    variants_df = pd.DataFrame(variant_sites, columns=['A','C','G','T'])
    variants_df = pd.concat([weights_df.ref,
                             weights_df.pos,
                             variants_df,
                             weights_df.depth,
                             weights_df.consensus,
                             weights_df.shannon], axis=1)
    if only_variants:
        variants_df = variants_df[variants_df['A'].notnull()
                                  | variants_df['C'].notnull()
                                  | variants_df['G'].notnull()
                                  | variants_df['T'].notnull()]
    return variants_df


def parse_samtools_depth(*args):
    ids_depths = {}
    for arg in args:
        id = arg
        depths_df = pd.read_table(arg, names=['contig', 'position', 'depth'])
        ids_depths[id] = depths_df.depth.tolist()
    return ids_depths


def plotly_samtools_depth(ids_depths):
    n_positions = len(ids_depths[max(ids_depths, key=lambda x: len(set(ids_depths[x])))])
    traces = []
    for id, depths in sorted(ids_depths.items()):
        traces.append(
            go.Scattergl(
                x=list(range(1, n_positions)),
                y=depths,
                mode='lines',
                name=id,
                text=id))
    layout = go.Layout(
        title='Depth of coverage',
        xaxis=dict(
            title='Position',
            gridcolor='rgb(255, 255, 255)',
            gridwidth=2),
        yaxis=dict(
            title='Depth',
            gridcolor='rgb(255, 255, 255)',
            gridwidth=2,
            type='log'),
        paper_bgcolor='rgb(243, 243, 243)',
        plot_bgcolor='rgb(243, 243, 243)')

    fig = go.Figure(data=traces, layout=layout)
    py.plot(fig, filename='depths.html')


def parse_variants(*args):
    ids_data = {}
    for arg in args:
        id = arg
        df = pd.read_table(arg, sep='\t')
        df['max_alt'] = df[['A', 'C', 'G', 'T']].max(axis=1)
        ids_data[id] = df.to_dict('series')
    return ids_data


def plotly_variants(ids_data):
    import plotly.offline as py
    import plotly.graph_objs as go
    traces = []
    for id, data in sorted(ids_data.items()):
        traces.append(
            go.Scattergl(
                x=data['pos'],
                y=data['max_alt'],
                mode='markers',
                name=id,
                text=id))
    layout = go.Layout(
        title='Variants',
        xaxis=dict(
            title='Position',
            gridcolor='rgb(255, 255, 255)',
            gridwidth=2),
        yaxis=dict(
            title='Abundance',
            gridcolor='rgb(255, 255, 255)',
            gridwidth=2,
            type='linear'),
        paper_bgcolor='rgb(243, 243, 243)',
        plot_bgcolor='rgb(243, 243, 243)')

    fig = go.Figure(data=traces, layout=layout)
    py.plot(fig, filename='variants.html')


# def samtools_depth(bam_path):
#     '''Return list of samtools reported depths (all positions, coverage limit 1m)'''
#     cmd = subprocess.run('samtools depth -a -m 1000000 {}'.format(bam_path),
#                          shell=True,
#                          check=True,
#                          stdout=subprocess.PIPE,
#                          universal_newlines=True)
#     depths_fh = io.StringIO(cmd.stdout)
#     df = pd.read_table(depths_fh, sep='\t', names=['contig', 'pos', 'depth'])
#     return df.depth.tolist()


def plotly_clips(bam_path):
    import plotly.offline as py
    import plotly.graph_objs as go
    aln = list(parse_bam(bam_path).items())[0][1]
    # depth = samtools_depth(bam_path)
    aligned_depth = [sum(weight.values()) for weight in aln.weights]
    # print(True if depth == aligned_depth else False)
    ins = [sum(i.values()) for i in aln.insertions]
    x_axis = list(range(len(aligned_depth)))
    traces = [
        # go.Scattergl(
        #     x = x_axis,
        #     y = depth,
        #     mode = 'lines',
        #     name = 'Samtools depth'),
        go.Scattergl(
            x = x_axis,
            y = aligned_depth,
            mode = 'lines',
            name = 'Aligned depth'),
        go.Scattergl(
            x = x_axis,
            y = aln.clip_starts,
            mode = 'markers',
            name = 'Soft clip starts'),
        go.Scattergl(
            x = x_axis,
            y = aln.clip_ends,
            mode = 'markers',
            name = 'Soft clip ends'),
        go.Scattergl(
            x = x_axis,
            y = aln.clip_start_cov,
            mode = 'lines',
            name = 'Soft clip start coverage'),
        go.Scattergl(
            x = x_axis,
            y = aln.clip_end_cov,
            mode = 'lines',
            name = 'Soft clip end coverage'),
        go.Scattergl(
            x = x_axis,
            y = ins,
            mode = 'lines',
            name = 'Insertions'),
        go.Scattergl(
            x = x_axis,
            y = aln.deletions,
            mode = 'lines',
            name = 'Deletions')]
    layout = go.Layout(
        xaxis=dict(
            type='linear',
            autorange=True),
        yaxis=dict(
            type='linear',
            autorange=True))
    fig = go.Figure(data=traces, layout=layout)
    out_fn = os.path.splitext(os.path.split(bam_path)[1])[0]
    py.plot(fig, filename=out_fn + '.clips.html')



if __name__ == '__main__':
    kindel.cli.main()

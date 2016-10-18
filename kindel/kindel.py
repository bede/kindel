# TODO
# - reporting:
#   - mix/max coverage
#   - consensus_sequence()… namedtuple causing slowdown?

import os
import sys
import collections

import argh
import tqdm
import simplesam

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_alignment(bam_path):
    '''
    Iterate over records, returning lists of base frequency dicts, indels and soft clpping info
    weights: base frequencies after CIGAR reconciliation
    insertions, deletions: per base frequencies
    gap_starts, gap_ends: gap positions in a LTR direction
    clip_s_weights, clip_e_weights: base frequencies from left and right clipped regions
    '''
    with open(bam_path, 'r') as bam_fh:
        records = simplesam.Reader(bam_fh)
        first_sq = list(records.header['@SQ'].values())[0] if '@SQ' in records.header else None
        ref_name = os.path.splitext(os.path.basename(bam_path))[0]
        ref_len = int(next(iter(first_sq)).replace('LN:','')) if first_sq else 100000
        weights = [{'A':0,'T':0,'G':0,'C':0,'N':0} for p in range(ref_len)]
        clip_s_weights = [{'A':0,'T':0,'G':0,'C':0,'N':0} for p in range(ref_len)]
        clip_e_weights = [{'A':0,'T':0,'G':0,'C':0,'N':0} for p in range(ref_len)]
        insertions = [collections.defaultdict(int) for p in range(ref_len)]
        deletions = [0] * ref_len
        clip_starts = [0] * (ref_len + 1)
        clip_ends = [0] * (ref_len + 1)
        for record in tqdm.tqdm(records, desc='loading sequences'): # Progress bar
            q_pos = 0 
            r_pos = record.pos-1 # Zero indexed genome coordinates
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
                    if i == 0: # Count right-of-gap / l-clipped start positions (e.g. start of ref)
                        clip_ends[r_pos] += 1
                        for gap_i in range(length):
                            q_nt = record.seq[gap_i].upper()
                            rel_r_pos = r_pos - length + gap_i
                            if rel_r_pos >= 0:
                                clip_e_weights[rel_r_pos][q_nt] += 1
                        q_pos += length
                    else: # Count left-of-gap / r-clipped start position (e.g. end of ref)
                        clip_starts[r_pos] += 1 
                        for pos in range(length):
                            q_nt = record.seq[q_pos].upper()
                            if r_pos < ref_len:
                                clip_s_weights[r_pos][q_nt] += 1
                                r_pos += 1
                                q_pos += 1
    return (ref_name, weights, insertions, deletions, clip_starts, clip_ends, clip_s_weights,
           clip_e_weights)


def find_gaps(weights, clip_starts, clip_ends, threshold_weight, min_depth):
    '''
    Return list of inclusive soft-clipped consensus gap coordinates as tuples. 
    '''
    gaps = []
    coverage = [sum({nt:w[nt] for nt in list('ACGT')}.values()) for w in weights]
    gap_open = False
    for i, (cov, clip_s, clip_e) in enumerate(zip(coverage, clip_starts, clip_ends)):
        threshold_weight_freq = max(cov*threshold_weight, min_depth)
        if clip_s > threshold_weight_freq and i:
            gap_start = i
            gap_open = True
        elif gap_open and clip_e > threshold_weight_freq and i:
            gap_end = i
            gaps.append((gap_start, gap_end))
            gap_open = False
    return gaps


def consensus(weight):
    '''
    Returns namedtuple of consensus base, weight and flag indicating a tie for consensus
    '''
    consensus_base, consensus_weight = max(weight.items(), key=lambda x:x[1])
    weight_sans_consensus = {k:d for k, d in weight.items() if k != consensus_base}
    tie = True if consensus_weight in weight_sans_consensus.values() else False
    pos_consensus = collections.namedtuple('pos_consensus', ['base', 'weight', 'tie'])
    return pos_consensus(consensus_base, consensus_weight, tie)


def s_overhang_consensus(clip_s_weights, start_pos, min_depth, max_len=500):
    '''
    Returns consensus sequence (string) of clipped reads at specified position
    start_pos is the first position described by the CIGAR-S
    '''
    consensus_overhang = ''
    for pos in range(start_pos, start_pos+max_len):
        pos_consensus = consensus(clip_s_weights[pos])
        if pos_consensus.weight >= min_depth:
            consensus_overhang += pos_consensus.base
        else:
            break
    return consensus_overhang


def e_overhang_consensus(clip_e_weights, start_pos, min_depth, max_len=500):
    '''
    Returns consensus sequence (string) of clipped reads at specified position
    end_pos is the last position described by the CIGAR-S
    '''
    rev_consensus_overhang = ''
    for pos in range(start_pos, start_pos-max_len, -1):
        pos_consensus = consensus(clip_e_weights[pos])
        if pos_consensus.weight >= min_depth:
            rev_consensus_overhang += pos_consensus.base
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
        if pos_consensus.weight >= min_depth:
            flank_seq += pos_consensus.base
    return flank_seq


def e_flanking_seq(end_pos, weights, min_depth, k):
    '''
    Returns consensus sequence (string) flanking RHS of soft-clipped gaps
    '''
    flank_seq = ''
    for pos in range(end_pos+1, end_pos+k+1):
        pos_consensus = consensus(weights[pos])
        if pos_consensus.weight >= min_depth:
            flank_seq += pos_consensus.base

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
    '''
    _lcs = lcs(l_seq, r_seq)
    l_trim = l_seq[:l_seq.find(_lcs)]
    r_trim = r_seq[r_seq.find(_lcs)+len(_lcs):]
    return l_trim + _lcs + r_trim


def reconcile_gaps(gaps, weights, clip_s_weights, clip_e_weights, min_depth, closure_k, uppercase):
    '''
    Returns dict of consensus insertions between to gap coordinates from find_gaps()
    Dict keys are gap start positions (left of gap)
    '''
    gap_consensuses = {}
    for gap in gaps:
        s_overhang_seq = s_overhang_consensus(clip_s_weights, gap[0], min_depth)
        e_overhang_seq = e_overhang_consensus(clip_e_weights, gap[1], min_depth)
        s_flank_seq = s_flanking_seq(gap[0], weights, min_depth, closure_k)
        e_flank_seq = e_flanking_seq(gap[1], weights, min_depth, closure_k)
        if e_flank_seq in s_overhang_seq: # Close gap using right-clipped read consensus
            i = s_overhang_seq.find(e_flank_seq) # str.find() returns -1 in absence of match
            gap_consensus = s_overhang_seq[:i]
        elif s_flank_seq in e_overhang_seq: # Close gap using left-clipped read consensus
            i = e_overhang_seq.find(s_flank_seq)
            gap_consensus = e_overhang_seq[i:]
        elif len(lcs(s_overhang_seq, e_overhang_seq)) >= closure_k:
            gap_consensus = close_by_lcs(s_overhang_seq, e_overhang_seq)
        else:
            print('Failed to close gap') # Stub... Needs tests
            break
        if uppercase:
            gap_consensuses[gap[0]] = gap_consensus
        else:
            gap_consensuses[gap[0]] = gap_consensus.lower()
    return gap_consensuses


def consensus_sequence(weights, clip_s_weights, clip_e_weights, insertions, deletions, gaps,
                       gap_consensuses, fix_gaps, trim_ends, threshold_weight, min_depth):
    consensus_seq = ''
    changes = [None] * len(weights)
    gap_starts = [g[0] for g in gaps] if fix_gaps else []
    skip_pos = False

    for pos, weight in tqdm.tqdm(enumerate(weights), total=len(weights), desc='building consensus'):
        ins_freq = sum(insertions[pos].values()) if insertions[pos] else 0
        del_freq = deletions[pos]
        coverage = sum({nt: weight[nt] for nt in list('ACGT')}.values())
        threshold_weight_freq = coverage * threshold_weight
        if pos in gap_starts and pos in gap_consensuses:
            consensus_seq += gap_consensuses[pos]
            gap_i = gap_starts.index(pos)
            skip_pos = gaps[gap_i][1]-pos
        elif skip_pos:
            skip_pos -= 1
            continue
        elif del_freq > threshold_weight_freq:
            changes[pos] = 'D'
        elif coverage < min_depth:
            consensus_seq += 'N'
            changes[pos] = 'N'
        else:
            pos_consensus = consensus(weight)
            consensus_seq += pos_consensus.base if not pos_consensus.tie else 'N'
        if ins_freq > threshold_weight_freq:
            insertion = consensus(insertions[pos])
            consensus_seq += insertion.base if not insertion.tie else 'N'
            changes[pos] = 'I'
    if trim_ends:
        consensus_seq = consensus_seq.strip('N')
    return consensus_seq, changes


def consensus_seqrecord(consensus, ref_name):
    return SeqRecord(Seq(consensus), id=ref_name + '_cns', description='')


def build_report(weights, changes, gaps, gap_consensuses, bam_path, fix_gaps, trim_ends,
                 threshold_weight, min_depth, closure_k, uppercase):
    coverage = [sum({nt:w[nt] for nt in list('ACGT')}.values()) for w in weights]
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
    report += '- threshold_weight: {}\n'.format(threshold_weight)
    report += '- min_depth: {}\n'.format(min_depth)
    report += '- closure_k: {}\n'.format(closure_k)
    report += 'min,max observed depth: {},{}\n'.format(min(coverage), max(coverage))
    report += 'ambiguous sites: {}\n'.format(', '.join(ambiguous_sites))
    report += 'insertion sites: {}\n'.format(', '.join(insertion_sites))
    report += 'deletion sites: {}\n'.format(', '.join(deletion_sites))
    report += 'soft-clipped gaps: {}\n'.format(', '.join(gaps_fmt))
    report += 'inserted sequences: {}\n'.format(gap_consensuses_fmt)
    report += '============================================================\n'
    return report


def bam_to_consensus(bam_path,
                     fix_gaps=False,
                     trim_ends=False,
                     threshold_weight=0.5,
                     min_depth=2,
                     closure_k=7,
                     uppercase=False):
    
    (ref_name, weights, insertions, deletions, clip_starts, clip_ends, clip_s_weights,
     clip_e_weights) = parse_alignment(bam_path)

    if fix_gaps:
        gaps = find_gaps(weights, clip_starts, clip_ends, threshold_weight, min_depth)
        gap_consensuses = reconcile_gaps(gaps, weights, clip_s_weights, clip_e_weights, min_depth,
                                         closure_k, uppercase)
    else:
        gaps, gap_consensuses = None, None
    
    consensus, changes = consensus_sequence(weights, clip_s_weights, clip_e_weights, insertions,
                                            deletions, gaps, gap_consensuses, fix_gaps, trim_ends,
                                            threshold_weight, min_depth)
    
    consensus_record = consensus_seqrecord(consensus, ref_name)
    
    report = build_report(weights, changes, gaps, gap_consensuses, bam_path, fix_gaps, trim_ends,
                          threshold_weight, min_depth, closure_k, uppercase)
    
    result = collections.namedtuple('result', ['record', 'report'])

    return result(consensus_record, report)


def bam_to_consensus_fasta(bam_path: 'path to SAM/BAM file',
                           fix_gaps: 'attempt to reconcile reference at soft-clip boundaries'=False,
                           trim_ends: 'trim ambiguous nucleotides (Ns) from sequence ends'=False,
                           threshold_weight: 'consensus threshold weight'=0.5,
                           min_depth: 'substitute Ns at coverage depths beneath this value'=2,
                           closure_k: 'match length required to close soft-clipped gaps'=7,
                           uppercase: 'close gaps using uppercase alphabet'=False):
    result = bam_to_consensus(bam_path, fix_gaps, trim_ends, threshold_weight, min_depth, closure_k,
                              uppercase)
    print(result.report, file=sys.stderr)
    return result.record.format('fasta')


if __name__ == '__main__':
    argh.dispatch_command(bam_to_consensus_fasta)

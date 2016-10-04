import os
import sys
import collections

import argh
import tqdm
import simplesam

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_records(bam_path):
    with open(bam_path, 'r') as bam_fh:
        records = simplesam.Reader(bam_fh)
        first_sq = list(records.header['@SQ'].values())[0] if '@SQ' in records.header else None
        ref_name = os.path.splitext(os.path.basename(bam_path))[0]
        ref_len = int(next(iter(first_sq)).replace('LN:','')) if first_sq else 100000
        weights = [{'A':0,'T':0,'G':0,'C':0,'N':0} for p in range(ref_len)]
        l_clip_weights = [{'A':0,'T':0,'G':0,'C':0,'N':0} for p in range(ref_len)]
        r_clip_weights = [{'A':0,'T':0,'G':0,'C':0,'N':0} for p in range(ref_len)]
        insertions = [collections.defaultdict(int) for p in range(ref_len)]
        deletions = [0] * ref_len
        clip_starts = [[0] * (ref_len+1), [0] * (ref_len+1)] # Genome len list for l/r-clipped seqs
        for record in tqdm.tqdm(records):
            q_pos = 0 
            r_pos = record.pos-1 # Zero indexed genome coordinates
            for i, cigarette in enumerate(record.cigars):
                length, operation = cigarette
                if operation == 'M':
                    for pos in range(length):
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
                    if i == 0:
                        clip_starts[0][r_pos] += 1 # Count left-clipped start position
                        for pos in reversed(range(length)): # Count backwards
                            q_nt = record.seq[pos].upper()
                            if r_pos >= 0:
                                l_clip_weights[r_pos][q_nt] += 1
                                r_pos -= 1
                                q_pos -= 1
                        # print(record.qname, record.seq[length-3]+record.seq[length-2]+record.seq[length-1]+record.seq[length]+record.seq[length+1]+record.seq[length+2])

                    else:
                        clip_starts[1][r_pos] += 1 # Count right-clipped start position
                        for pos in range(length): # Count forwards
                            q_nt = record.seq[q_pos].upper()
                            if r_pos < ref_len:
                                r_clip_weights[r_pos][q_nt] += 1
                                r_pos += 1
                                q_pos += 1
                    # q_pos += length
        # for i, w in enumerate(l_clip_weights):
        #     print(max(w, key=lambda k: w[k]), end='')

    return ref_name, weights, insertions, deletions, clip_starts, r_clip_weights


def find_gaps(weights, clip_starts, threshold_weight, min_depth):
    # Returns list of soft-clipped alignment gaps as tuples in format [(start_pos, end_pos)]
    gaps = []
    weights_acgt = [{nt: weights[i][nt] for nt in list('ACGT')} for i in range(len(weights))]
    coverage = [sum(weight.values()) for weight in weights_acgt]
    for i, (c, l, r) in enumerate(zip(coverage, clip_starts[0], clip_starts[1])):
        gap_open = False
        threshold_weight_freq = max(c * threshold_weight, min_depth)
        if r > threshold_weight_freq and i:
            gap_start = i
            gap_open = True
        if gap_open and l > threshold_weight_freq and i:
            gap_end = i
            gaps.append((gap_start, gap_end))
            gap_open = False
    return gaps

def l_overhang_consensus(l_clip_weights, start_pos, min_depth):
    # Returns consensus sequence (string) of clipped reads at specified position
    consensus_overhang = ''
    max_overhang_len = 500 # Arbitrary figure greater than read length as safety net
    for pos in range(start_pos, start_pos-max_overhang_len, -1):
        heaviest_base, heaviest_weight = max(l_clip_weights[pos].items(), key=lambda x:x[1])
        if heaviest_weight >= min_depth:
            consensus_overhang += heaviest_base
        else:
            break
    return consensus_overhang

def r_overhang_consensus(r_clip_weights, start_pos, min_depth):
    # Returns consensus sequence (string) of clipped reads at specified position
    consensus_overhang = ''
    max_overhang_len = 500 # Arbitrary figure greater than read length as safety net
    for pos in range(start_pos, start_pos+max_overhang_len):
        heaviest_base, heaviest_weight = max(r_clip_weights[pos].items(), key=lambda x:x[1])
        if heaviest_weight >= min_depth:
            consensus_overhang += heaviest_base
        else:
            break
    return consensus_overhang

def r_flanking_seq(start_pos, weights, min_depth, k):
    # Returns consensus sequence (string) flanking RHS of soft-clipped gaps
    r_flank_seq = ''
    for pos in range(start_pos, start_pos+k):
        heaviest_base, heaviest_weight = max(weights[pos].items(), key=lambda x:x[1])
        if heaviest_weight >= min_depth:
            r_flank_seq += heaviest_base
    return r_flank_seq

def reconcile_gaps(gaps, weights, r_clip_weights, min_depth, bridge_k):
    # Returns list of consensus strings corresponding to gap coordinates generated by find_gaps()
    gap_consensuses = []
    for gap in gaps:
        r_overhang = r_overhang_consensus(r_clip_weights, gap[0], min_depth)
        r_flank_seq = r_flanking_seq(gap[1], weights, min_depth, bridge_k)
        index = r_overhang.find(r_flank_seq) # str.find() returns -1 in absence of match
        gap_consensuses.append(r_overhang[:index].lower() if index >= 0 else '')
        # print(r_overhang, file=sys.stderr)
        # print(r_flank_seq, file=sys.stderr)
        # print(index, file=sys.stderr)
    return gap_consensuses


def consensus_sequence(weights, r_clip_weights, insertions, deletions, gaps, gap_consensuses,
                       fix_gaps, trim_ends, threshold_weight, min_depth):
    consensus = ''
    changes = [None] * len(weights)
    gap_starts = [g[0] for g in gaps] if fix_gaps else []
    skip_pos = False
    for pos, weight in enumerate(weights):
        ins_freq = sum(insertions[pos].values()) if insertions[pos] else 0
        del_freq = deletions[pos]
        coverage = sum({nt: weight[nt] for nt in list('ACGT')}.values())
        threshold_weight_freq = coverage * threshold_weight
        if pos in gap_starts and gap_consensuses[gap_starts.index(pos)]:
            gap_i = gap_starts.index(pos)
            consensus += gap_consensuses[gap_i]
            skip_pos = gaps[gap_i][1]-gaps[gap_i][0] 
        elif skip_pos:
            skip_pos -= 1
            continue
        elif del_freq > threshold_weight_freq:
            changes[pos] = 'D'
        elif coverage < min_depth:
            consensus += 'N'
            changes[pos] = 'N'
        else:
            consensus += max(weight, key=lambda k: weight[k])
        if ins_freq > threshold_weight_freq:
            top_ins, top_ins_freq = max(insertions[pos].items(), key=lambda x:x[1])
            consensus += top_ins
            changes[pos] = 'I'
    if trim_ends:
        consensus = consensus.strip('N')
    return consensus, changes


def consensus_seqrecord(consensus, ref_name):
    return SeqRecord(Seq(consensus), id=ref_name + '_cns', description='')


def build_report(weights, changes, gaps, gap_consensuses, bam_path, fix_gaps, trim_ends,
                 threshold_weight, min_depth, bridge_k):
    # weights_acgt = [{nt: weights[i][nt] for nt in list('ACGT')} for i in range(len(weights))]
    # coverage = [sum(weight.values()) for weight in weights_acgt]

    coverage = [sum({nt:w[nt] for nt in list('ACGT')}.values()) for w in weights]


    for i, cov in enumerate(coverage):
        print(str(i), cov)
    # print(coverage)
    ambiguous_sites = []
    insertion_sites = []
    deletion_sites = []
    gaps_fmt = ['-'.join([str(g) for g in gap]) for gap in gaps] if gaps else []
    gap_consensuses_fmt = ', '.join(gap_consensuses) if gap_consensuses else ''
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
    report += '- threshold_weight: {}\n'.format(threshold_weight)
    report += '- min_depth: {}\n'.format(min_depth)
    report += '- bridge_k: {}\n'.format(bridge_k)
    report += 'min,max observed depth: {},{}\n'.format(min(coverage), max(coverage))
    report += 'ambiguous sites: {}\n'.format(', '.join(ambiguous_sites))
    report += 'insertion sites: {}\n'.format(', '.join(insertion_sites))
    report += 'deletion sites: {}\n'.format(', '.join(deletion_sites))
    report += 'soft-clipped gaps: {}\n'.format(', '.join(gaps_fmt))
    report += 'inserted sequences: {}\n'.format(gap_consensuses_fmt)
    report += '============================================================\n'
    return report


def bam_to_consensus_seqrecord(bam_path,
                               fix_gaps=False,
                               trim_ends=False,
                               threshold_weight=0.5,
                               min_depth=2,
                               bridge_k=7):
    ref_name, weights, insertions, deletions, clip_starts, r_clip_weights = parse_records(bam_path)
    if fix_gaps:
        gaps = find_gaps(weights, clip_starts, threshold_weight, min_depth)
        gap_consensuses = reconcile_gaps(gaps, weights, r_clip_weights, min_depth, bridge_k)
    else: gaps, gap_consensuses = None, None
    consensus, changes = consensus_sequence(weights, r_clip_weights, insertions, deletions, gaps,
                                            gap_consensuses, fix_gaps, trim_ends, threshold_weight,
                                            min_depth)
    consensus_record = consensus_seqrecord(consensus, ref_name)
    report_fmt = build_report(weights, changes, gaps, gap_consensuses, bam_path, fix_gaps,
                              trim_ends, threshold_weight, min_depth, bridge_k)

    return consensus_record, report_fmt


def bam_to_consensus_fasta(bam_path: 'path to SAM/BAM file',
                           fix_gaps: 'attempt to reconcile reference at soft-clip boundaries'=False,
                           trim_ends: 'trim ambiguous nucleotides (Ns) from sequence ends'=False,
                           threshold_weight: 'consensus threshold weight'=0.5,
                           min_depth: 'substitute Ns at coverage depths beneath this value'=2,
                           bridge_k: 'match length required to bridge soft-clipped gaps'=7):
    consensus_record, report_fmt = bam_to_consensus_seqrecord(bam_path,
                                                              fix_gaps,
                                                              trim_ends,
                                                              threshold_weight,
                                                              min_depth,
                                                              bridge_k)
    print(report_fmt, file=sys.stderr)
    return consensus_record.format('fasta')


if __name__ == '__main__':
    argh.dispatch_command(bam_to_consensus_fasta)

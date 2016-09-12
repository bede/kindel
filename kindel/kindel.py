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
        ref_name = list(list(records.header.values())[0].keys())[0].replace('SN:','') if first_sq else 'aln'
        ref_len = int(next(iter(first_sq)).replace('LN:','')) if first_sq else 100000
        weights = [{'A':0,'T':0,'G':0,'C':0,'N':0} for p in range(ref_len)]
        clip_weights = [{'A':0,'T':0,'G':0,'C':0,'N':0} for p in range(ref_len)]
        insertions = [collections.defaultdict(int) for p in range(ref_len)]
        deletions = [0] * ref_len
        clip_starts = [[0] * (ref_len+1), [0] * (ref_len+1)] # Genome length lists for left- and right-clipped seqs
        for record in tqdm.tqdm(records):
            q_pos = 0 
            r_pos = record.pos-1 # Use zero-based coordinates
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
                    else:
                        # print(clip_starts)
                        # print('\n', r_pos, length, record.qname)
                        clip_starts[1][r_pos] += 1 # Count right-clipped start position
                        for pos in range(length):
                            # print(q_pos, r_pos)
                            q_nt = record.seq[q_pos].upper()
                            if r_pos < ref_len:
                                clip_weights[r_pos][q_nt] += 1
                                r_pos += 1
                                q_pos += 1
                    q_pos += length
    return ref_name, weights, insertions, deletions, clip_starts, clip_weights


def clipping_boundaries(weights, clip_starts, threshold_weight, min_depth):
    clip_boundaries = {'l': [], 'r': []}
    coverage = [sum(weight.values()) for weight in weights]
    for i, (c, l, r) in enumerate(zip(coverage, clip_starts[0], clip_starts[1])):
        threshold_weight_freq = max(c * threshold_weight, min_depth)
        if l > threshold_weight_freq and i:
            clip_boundaries['l'].append(i)
        if r > threshold_weight_freq and i:
            clip_boundaries['r'].append(i)
    return clip_boundaries


def consensus_sequence(weights, insertions, deletions, clip_boundaries, threshold_weight, min_depth):
    consensus = ''
    changes = [None] * len(weights)
    for pos, weight in enumerate(weights):


        if pos in clip_boundaries[1]:
            pass # Do reconciliation based on clip_weights 


        ins_freq = sum(insertions[pos].values()) if insertions[pos] else 0
        del_freq = deletions[pos]
        coverage = sum(weight.values())
        threshold_weight_freq = coverage * threshold_weight
        if del_freq > threshold_weight_freq:
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
    return consensus, changes


def consensus_seqrecord(consensus, ref_name):
    return SeqRecord(Seq(consensus), id=ref_name + '_cns', description='')


def report(weights, changes, threshold_weight, min_depth):
    coverage = [sum(weight.values()) for weight in weights] 
    ambiguous_sites = []
    insertion_sites = []
    deletion_sites = []
    for pos, change in enumerate(changes):
        if change == 'N':
            ambiguous_sites.append(str(pos))
        elif change == 'I':
            ambiguous_sites.append(str(pos+1))
        elif change == 'D':
            ambiguous_sites.append(str(pos))
    report = '========================= REPORT ===========================\n'
    report += 'consensus weight: {}\n'.format(threshold_weight)
    report += 'minimum depth: {}\n'.format(min_depth)
    report += 'min,max observed depth: {},{}\n'.format(min(coverage), max(coverage))
    report += 'ambiguous sites: {}\n'.format(', '.join(ambiguous_sites))
    report += 'insertion sites: {}\n'.format(', '.join(insertion_sites))
    report += 'deletion sites: {}\n'.format(', '.join(deletion_sites))
    report += '============================================================\n'
    return report


def bam_to_consensus_seqrecord(bam_path, threshold_weight=0.5, min_depth=1):
    ref_name, weights, insertions, deletions, clip_starts, clip_weights = parse_records(bam_path)
    print('\n'.join([str(x) for x in clip_weights]))
    clip_boundaries = clipping_boundaries(weights, clip_starts,
                                          threshold_weight, min_depth)
    print('clip_boundaries')
    print(clip_boundaries)
    consensus, changes = consensus_sequence(weights, insertions, deletions, clip_boundaries, threshold_weight, min_depth)
    consensus_record = consensus_seqrecord(consensus, ref_name)
    return consensus_record


def bam_to_consensus_fasta(bam_path: 'path to SAM/BAM file',
                           threshold_weight: 'consensus threshold weight'=0.5,
                           min_depth: 'substitute Ns at coverage depths beneath this value'=1):
    consensus_fasta = bam_to_consensus_seqrecord(bam_path, threshold_weight, min_depth).format('fasta')
    return consensus_fasta


if __name__ == '__main__':
    argh.dispatch_command(bam_to_consensus_fasta)

#!/usr/bin/env python3

import sys

import tqdm
import simplesam

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from collections import defaultdict


def parse_records(bam_path):
    with open(bam_path, 'r') as bam_fh:
        records = simplesam.Reader(bam_fh)
        first_sq = list(records.header['@SQ'].values())[0] if '@SQ' in records.header else None
        ref_name = list(list(records.header.values())[0].keys())[0].replace('SN:','') if first_sq else 'aln'
        ref_len = int(next(iter(first_sq)).replace('LN:',''))+3 if first_sq else 100000
        weights = [{'A':0,'T':0,'G':0,'C':0,'N':0} for p in range(ref_len)]
        insertions = [defaultdict(int) for p in range(ref_len)]
        deletions = [0] * ref_len
        for i, record in tqdm.tqdm(enumerate(records)):
            q_pos = 0 
            r_pos = record.pos-1
            first_iteration = True
            for cigarette in record.cigars:
                length, operation = cigarette
                if operation == 'M':
                    for pos in range(length):
                        q_nt = record.seq[q_pos]
                        weights[r_pos][q_nt] += 1
                        r_pos += 1
                        q_pos += 1
                elif operation == 'I':
                    nts = record.seq[q_pos:q_pos+length]
                    insertions[r_pos][nts] += 1
                    q_pos += length
                elif operation == 'D':
                    deletions[r_pos] += 1
                    r_pos += length
                elif operation == 'S':
                    if first_iteration:
                        r_pos += length
                        q_pos += length
                    else:
                        q_pos += length
                first_iteration = False
    return weights, insertions, deletions, ref_name


def consensus_sequence(weights, insertions, deletions, threshold_weight, min_depth):
    consensus = ''
    changes = [None] * len(weights)
    for pos, weight in enumerate(weights):
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
    weights, insertions, deletions, ref_name = parse_records(bam_path)
    consensus, changes = consensus_sequence(weights, insertions, deletions, threshold_weight, min_depth)
    consensus_record = consensus_seqrecord(consensus, ref_name)
    return consensus_record


def bam_to_consensus_fasta(bam_path, fasta_path=sys.stdout, threshold_weight=0.5, min_depth=1):
    weights, insertions, deletions, ref_name = parse_records(bam_path)
    consensus, changes = consensus_sequence(weights, insertions, deletions, threshold_weight, min_depth)
    consensus_record = consensus_seqrecord(consensus, ref_name)
    SeqIO.write(consensus_record, fasta_path, format='fasta')


if __name__ == '__main__':

    bam_path = sys.argv[1]
    threshold_weight = 0.5
    min_depth = 1
    weights, insertions, deletions, ref_name = parse_records(bam_path)
    consensus, changes = consensus_sequence(weights, insertions, deletions, threshold_weight, min_depth)
    consensus_record = consensus_seqrecord(consensus, ref_name)
    SeqIO.write(consensus_record, sys.stdout, format='fasta')
    print(report(weights, changes, threshold_weight, min_depth), file=sys.stderr, end='')

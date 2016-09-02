#!/usr/bin/env python3

import sys
import simplesam

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from collections import defaultdict


def parse_records(sam_path):
    with open(sam_path, 'r') as sam_fh:
        records = simplesam.Reader(sam_fh)
        first_sq = list(records.header['@SQ'].values())[0] if '@SQ' in records.header else None
        ref_name = list(list(records.header.values())[0].keys())[0].replace('SN:','') if first_sq else 'aln'
        ref_len = int(next(iter(first_sq)).replace('LN:',''))+3 if first_sq else 100000
        weights = [{'A':0,'T':0,'G':0,'C':0,'N':0} for p in range(ref_len)]
        insertions = [defaultdict(int) for p in range(ref_len)]
        deletions = [0] * ref_len
        for i, record in enumerate(records):
            if i % 100000 == 0 and i:
                print('…Processed', str(i), 'records', file=sys.stderr) 
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


def consensus_sequence(weights, insertions, deletions, prop_threshold, min_coverage):
    consensus = ''
    changes = [None] * len(weights)
    for pos, weight in enumerate(weights):
        ins_freq = sum(insertions[pos].values()) if insertions[pos] else 0
        del_freq = deletions[pos]
        coverage = sum(weight.values())
        freq_threshold = coverage * prop_threshold
        if del_freq > freq_threshold:
            changes[pos] = 'D'
        elif coverage < min_coverage:
            consensus += 'N'
            changes[pos] = 'N'
        else:
            consensus += max(weight, key=lambda k: weight[k])
        if ins_freq > freq_threshold:
            top_ins, top_ins_freq = max(insertions[pos].items(), key=lambda x:x[1])
            consensus += top_ins
            changes[pos] = 'I'

    return consensus, changes

def consensus_seqrecord(consensus, ref_name):
    return SeqRecord(Seq(consensus), id=ref_name + '_cns', description='')

def report(changes, prop_threshold, min_coverage):
    ambiguous_sites = []
    insertion_sites = []
    deletion_sites = []
    for pos, change in enumerate(changes):
        if change == 'N':
            ambiguous_sites.append(str(pos))
        if change == 'I':
            ambiguous_sites.append(str(pos+1))
        if change == 'D':
            ambiguous_sites.append(str(pos))
    report = '========================= REPORT ===========================\n'
    report += 'consensus threshold: {}\n'.format(prop_threshold)
    report += 'minimum coverage: {}\n'.format(min_coverage)
    report += 'ambiguous sites: {}\n'.format(', '.join(ambiguous_sites))
    report += 'insertion sites: {}\n'.format(', '.join(insertion_sites))
    report += 'deletion sites: {}\n'.format(', '.join(deletion_sites))
    report += '============================================================\n'

    return report


if __name__ == '__main__':

    sam_path = sys.argv[1]
    prop_threshold = 0.5
    min_coverage = 1
    weights, insertions, deletions, ref_name = parse_records(sam_path)
    consensus, changes = consensus_sequence(weights, insertions, deletions, prop_threshold, min_coverage)    
    consensus_record = consensus_seqrecord(consensus, ref_name)
    SeqIO.write(consensus_record, sys.stdout, format='fasta')
    print(report(changes, prop_threshold, min_coverage))

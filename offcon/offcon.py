#!/usr/bin/env python3

import sys
import simplesam

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from collections import defaultdict


def parse_records(sam_path, ref_max_len=100000):
    
    with open(sam_path, 'r') as sam_fh:
        records = simplesam.Reader(sam_fh)

    weights = [{'A':0,'T':0,'G':0,'C':0,'N':0} for p in range(ref_max_len)]
    insertions = [defaultdict(int) for p in range(ref_max_len)]
    deletions = [0] * ref_max_len
    for i, record in enumerate(records):
        if i % 100000 == 0 and i:
            print('Processed', str(i), 'records...', file=sys.stderr) 
        q_pos = 0 
        r_pos = record.pos
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
                insertions[r_pos-1][nts] += 1
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
    return weights, insertions, deletions


def coverage(weights):
    coverage = [sum(w.values()) for w in weights]
    return coverage


def weights_prop(weights):
    cov = [sum(w.values()) for w in weights]
    weights_prop = [{nt:round(wt/cov[p], 3) if cov[p] else 0 for nt, wt in wts.items() if cov[p]}
                    for p, wts in enumerate(weights)]
    return weights_prop


def consensus_sequence(weights, insertions, deletions, min_coverage=1):
    consensus = ''
    ins_total, del_total, n_total = 0, 0, 0
    for pos, weight in enumerate(weights):
        ins_freq = sum(insertions[pos].values()) if insertions[pos] else 0
        del_freq = deletions[pos]
        coverage = sum(weight.values())
        consensus_threshold = coverage * 0.5
        if del_freq > consensus_threshold:
            del_total += 1
        elif coverage < min_coverage:
            consensus += 'N'
            n_total += 1
        else:
            consensus += max(weight, key=lambda k: weight[k])
        if ins_freq > consensus_threshold:
            top_ins, top_ins_freq = max(insertions[pos].items(), key=lambda x:x[1])
            consensus += top_ins
            ins_total += 1

    return consensus.strip('N')


if __name__ == '__main__':

    weights, insertions, deletions = parse_records(records)
    consensus = consensus_sequence(weights, insertions, deletions)
    
    consensus_record = SeqRecord(Seq(cns), id='cns', description='')
    SeqIO.write(cns_record, sys.stdout, format='fasta')

    # print(coverage(weights))
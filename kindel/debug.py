import sys

import kindel

from Bio import SeqIO


def show_weights(sample_ids):
    sample_weights = [SeqIO.read('tests/bam/master_cns/' + id + '.bam.cns.fa', 'fasta') for id in sample_ids]

    alignments = []
    for sample_id in sample_ids:
        alignments.append(kindel.parse_alignment('tests/bam/' + sample_id + '.bam'))

    for id, alignment in zip(sample_ids, alignments):
        print(id)
        for i, w in enumerate(alignment.weights, start=1):
            coverage = sum({nt:w[nt] for nt in list('ACGT')}.values())
            consensus = kindel.consensus(w)
            print(i,
                  coverage,
                  consensus.base,
                  consensus.frequency,
                  consensus.prop,
                  'TIE' if consensus.tie else '',
                  'DIVERGENT' if consensus.prop and consensus.prop <= 0.75 else '')


if __name__ == '__main__':

    sample_ids = ['HCV_AVU_AB_1.12345.R12.sub']
             # 'HCV_AVU_AB_3.4.R12.sub',
             # 'HCV_AVU_AB_7.1.R12.sub',
             # 'HCV_AVU_AB_9.1.R12.sub',
             # 'HCV_AVU_AB_11.1.R12.sub']

    show_weights(sample_ids)
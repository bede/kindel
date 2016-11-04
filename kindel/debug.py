import sys

import kindel

from Bio import SeqIO


def show_weights(sample_ids):
    # sample_weights = [SeqIO.read('tests/bam/test_cns/' + id + '.bam.cns.fa', 'fasta') for id in sample_ids]

    # print(sample_weights)

    alignments = []
    for sample_id in sample_ids:
        alignments.append(kindel.parse_alignment('tests/bam/' + sample_id + '.bam'))

    for id, alignment in zip(sample_ids, alignments):
        print(id)
        for i, w in enumerate(alignment.weights):
            coverage = sum({nt:w[nt] for nt in list('ACGT')}.values())
            consensus = kindel.consensus(w)
            print(i+1,
                  coverage,
                  consensus[0],
                  consensus[1],
                  consensus[2],
                  str(alignment.clip_starts[i]) + 'clip' + str(alignment.clip_ends[i]),
                  'TIE' if consensus[3] else '',
                  'DIVERGENT' if consensus[2] and consensus[2] <= 0.75 else '',
                  w)


if __name__ == '__main__':

    sample_ids = ['HCV_AVU_AB_1.12345.R12.sub']
             # 'HCV_AVU_AB_3.4.R12.sub',
             # 'HCV_AVU_AB_7.1.R12.sub',
             # 'HCV_AVU_AB_9.1.R12.sub',
             # 'HCV_AVU_AB_11.1.R12.sub']

    show_weights(sample_ids)
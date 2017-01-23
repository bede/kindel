import sys
import functools
import concurrent.futures

from kindel import kindel

from Bio import SeqIO

def run(function, iterable):
    with concurrent.futures.ProcessPoolExecutor() as x:
        return {iterable:r for iterable, r in zip(iterable, x.map(function, iterable))}


sample_ids = ['HCV_AVU_AB_1.12345.R12.sub',
             'HCV_AVU_AB_3.4.R12.sub',
             'HCV_AVU_AB_7.1.R12.sub',
             'HCV_AVU_AB_9.1.R12.sub',
             'HCV_AVU_AB_11.1.R12.sub']
# sample_bam_fns = [f'tests/bwa/{id}.bam' for id in sample_ids]
# sample_master_cns_fns = [f'tests/bwa/master_cns/{id}.bam.cns.fa' for id in sample_ids]


### UNIT ###

def test_consensus():
    pos_weight = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'N': 5}
    assert kindel.consensus(pos_weight)[0] == 'N'
    assert kindel.consensus(pos_weight)[1] == 5
    assert kindel.consensus(pos_weight)[2] == 0.33
    assert kindel.consensus(pos_weight)[3] is False
    pos_weight_tie = {'A': 5, 'C': 5, 'G': 3, 'T': 4, 'N': 1}
    assert kindel.consensus(pos_weight_tie)[2]


### FUNCTIONAL ###
# cd /Users/bede/Research/Tools/kindel/kindel/tests/bwa
# for file in *.bam; do echo $file; python3 /Users/Bede/Research/Tools/kindel/kindel/kindel.py --fix-gaps --trim-ends $file > master_cns/$file.cns.fa; done

def test_basic_bwa():
    bam_path = 'tests/bwa/' + sample_ids[0] + '.bam'
    assert kindel.bam_to_consensus(bam_path)

def test_outputs_fix_gaps_trim_ends_bwa(): # Slow AF
    test_records = [kindel.bam_to_consensus('tests/bwa/' + id + '.bam',
                                            fix_gaps=True,
                                            trim_ends=True)[0] for id in sample_ids]
    master_records = [SeqIO.read('tests/bwa/master_cns/' + id + '.bam.cns.fa', 'fasta') for id in sample_ids]
    for t, m in zip(test_records, master_records):
        if str(t.seq) != str(m.seq):
            print('test len:' + str(len(t.seq)), 'master len: ' + str(len(m.seq)))
            lcs = kindel.lcs(str(t.seq), str(m.seq))
            print(lcs)
            print('lcs len: ' + str(len(lcs)))
            # SeqIO.write(t, sys.stdout, 'fasta')
            # SeqIO.write(m, sys.stdout, 'fasta')

# def test_outputs_fix_gaps_trim_ends_bwa_parallel():
#     partial_bam_to_consensus = functools.partial(kindel.bam_to_consensus,
#                                                  fix_gaps=True,
#                                                  trim_ends=True,
#                                                  reporting=False)
#     print(sample_bam_fns)
#     test_records = run(partial_bam_to_consensus, sample_bam_fns)
#     print('=======')
#     print(test_records)
#     print('=======')
#     master_records = {fn:SeqIO.read(fn, 'fasta') for fn in sample_master_cns_fns}
#     # print(test_records, master_records)
#     for (tk, tv), (mk, mv), in zip(test_records.items(), master_records.items()):
#         print('>>>>>>>>>>>')
#         print(tk, tv, mk, mv, sep='\n')
#         print('>>>>>>>>>>>')
#         if str(t.seq) != str(m.seq):
#             print('test len:' + str(len(t.seq)), 'master len: ' + str(len(m.seq)))
#             lcs = kindel.lcs(str(t.seq), str(m.seq))
#             print(lcs)
#             print('lcs len: ' + str(len(lcs)))
            # SeqIO.write(t, sys.stdout, 'fasta')
            # SeqIO.write(m, sys.stdout, 'fasta')


def test_basic_multiple_contigs_bbmap():
    pass
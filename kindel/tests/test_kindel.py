import sys

from kindel import kindel

from Bio import SeqIO


sample_ids = ['HCV_AVU_AB_1.12345.R12.sub',
             'HCV_AVU_AB_3.4.R12.sub',
             'HCV_AVU_AB_7.1.R12.sub',
             'HCV_AVU_AB_9.1.R12.sub',
             'HCV_AVU_AB_11.1.R12.sub']

### UNIT ###

def test_consensus():
    pos_weight = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'N': 5}
    assert kindel.consensus(pos_weight).base == 'N'
    assert kindel.consensus(pos_weight).frequency == 5
    assert kindel.consensus(pos_weight).prop == 0.33
    assert kindel.consensus(pos_weight).tie is False
    pos_weight_tie = {'A': 5, 'C': 5, 'G': 3, 'T': 4, 'N': 1}
    assert kindel.consensus(pos_weight_tie)[2]


### FUNCTIONAL ###
# cd /Users/bede/Research/Tools/kindel/kindel/tests/bam
# for file in *.bam; do echo $file; python3 /Users/Bede/Research/Tools/kindel/kindel/kindel.py --fix-gaps --trim-ends $file > master_cns/$file.cns.fa; done

def test_basic():
    bam_path = 'tests/bam/' + sample_ids[0] + '.bam'
    assert kindel.bam_to_consensus(bam_path)

def test_outputs_fix_gaps_trim_ends():
    test_records = [kindel.bam_to_consensus('tests/bam/' + id + '.bam',
                                            fix_gaps=True,
                                            trim_ends=True)[0] for id in sample_ids]
    master_records = [SeqIO.read('tests/bam/master_cns/' + id + '.bam.cns.fa', 'fasta') for id in sample_ids]
    for t, m in zip(test_records, master_records):
        if str(t.seq) != str(m.seq):
            print('test len:' + str(len(t.seq)), 'master len: ' + str(len(m.seq)))
            lcs = kindel.lcs(str(t.seq), str(m.seq))
            print(lcs)
            print('lcs len: ' + str(len(lcs)))
            SeqIO.write(t, sys.stdout, 'fasta')
            SeqIO.write(m, sys.stdout, 'fasta')

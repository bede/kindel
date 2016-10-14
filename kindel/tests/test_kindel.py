from kindel import kindel

test_bams = ['HCV_AVU_AB_1.12345.R12.sub.bam',
             'HCV_AVU_AB_3.4.R12.sub.bam',
             'HCV_AVU_AB_7.1.R12.sub.bam',
             'HCV_AVU_AB_9.1.R12.sub.bam',
             'HCV_AVU_AB_11.1.R12.sub.bam']

test_bam_records = [kindel.bam_to_consensus_seqrecord('tests/' + b) for b in test_bams]

def test_bam():
    assert test_bam_records

def test_consensus():
    weight = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'N': 5}
    assert kindel.consensus(weight)[0] == 'N'
    assert kindel.consensus(weight)[1] == 5
    assert kindel.consensus(weight)[2] is False
    weight_tie = {'A': 5, 'C': 5, 'G': 3, 'T': 4, 'N': 1}
    assert kindel.consensus(weight_tie)[2] is True
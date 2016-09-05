from kindel import kindel

bam_consensus_record = kindel.bam_to_consensus_seqrecord('tests/hcv_k21c100.EU155341.bam')
sam_consensus_record = kindel.bam_to_consensus_seqrecord('tests/hcv_k21c100.EU155341.sam')

def test_bam():
    assert bam_consensus_record

# def test_sam():
#     assert sam_consensus_record

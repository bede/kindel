import os
import sys
import functools
import concurrent.futures

from Bio import SeqIO

from kindel import kindel


bwa_path = 'tests/data_bwa_mem/'
seg_path = 'tests/data_sehemehl/'

bwa_fns = [bwa_path + fn for fn in os.listdir(bwa_path) if fn.endswith('.bam')]
seg_fns = [seg_path + fn for fn in os.listdir(seg_path) if fn.endswith('.bam')]

test_aln = list(kindel.parse_bam(bwa_path + '1.1.sub_test.bam').values())[0]


# UNIT

def test_consensus():
    pos_weight = {'A': 1, 'C': 2, 'G': 3, 'T': 4, 'N': 5}
    assert kindel.consensus(pos_weight)[0] == 'N'
    assert kindel.consensus(pos_weight)[1] == 5
    assert kindel.consensus(pos_weight)[2] == 0.33
    assert kindel.consensus(pos_weight)[3] is False
    pos_weight_tie = {'A': 5, 'C': 5, 'G': 3, 'T': 4, 'N': 1}
    assert kindel.consensus(pos_weight_tie)[2]


def test_merge_by_lcs():
    one = ('AACTGCCGCTAGGGGCGCGTTCGGGCTCGCCAACATCTTCAGTCCGGG',
           'GCCGCTAGGGGCGCGTTCGGGCTCGCCAACATCTTCAGTCCGGGCGCTAAGCAGAACA')
    two = ('AACTGCCGCTAGGGGCGCGTTCGGGCTCGCCAACATCTTCAGTCCGGGCGCTAAGCAGAACATC',
           'GCAGATACCTACACCACCGGGGGAACTGCCGCTAGGGGCGCGTTCGGGCTCGCCAACATCTTCAGTCCGGGCGCTAAGCAGAACA')
    short = ('AT', 'CG')
    assert kindel.merge_by_lcs(*one) == 'AACTGCCGCTAGGGGCGCGTTCGGGCTCGCCAACATCTTCAGTCCGGGCGCTAAGCAGAACA'
    assert kindel.merge_by_lcs(*two) == 'AACTGCCGCTAGGGGCGCGTTCGGGCTCGCCAACATCTTCAGTCCGGGCGCTAAGCAGAACA'
    assert kindel.merge_by_lcs(*short) == None


# FUNCTIONAL

def test_parse_bam():
    assert test_aln.ref_id == 'ENA|EU155341|EU155341.2'
    assert len(test_aln.weights) == 9306 


def test_cdrp_consensuses():
    cdrps = kindel.cdrp_consensuses(test_aln.weights, test_aln.clip_start_weights,
                                    test_aln.clip_end_weights, test_aln.clip_start_depth,
                                    test_aln.clip_end_depth, 0.1, 10)
    print(cdrps)
    assert cdrps[0][0].seq == 'AACTGCCGCTAGGGGCGCGTTCGGGCTCGCCAACATCTTCAGTCCGGGCGCTAAGCAGAACATCCAGCTGATCAACA'
    assert cdrps[0][1].seq == 'AGCGTCGATGCAGATACCTACACCACCGGGGGAACTGCCGCTAGGGGCGCGTTCGGGCTCGCCAACATCTTCAGTCCGGGCGCTAAGCAGAACA'


def test_bam_to_consensus_bwa():
    for fn in bwa_fns:
        assert kindel.bam_to_consensus(fn)


def test_bam_to_consensus_realign_bwa():
    for fn in bwa_fns:
        assert kindel.bam_to_consensus(fn, realign=True)

def test_weights():
    kindel.weights(bwa_fns[0], relative=True) 


# CLI




# SAMPLE-SPECIFIC FUNCTIONAL REGRESSION
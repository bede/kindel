#!/usr/bin/env python3

from sod import sod

read_len = 250

records = sod.simulate_reads('tests/hxb2.fasta', n_reads=1000, read_len=read_len,
                             sub_rate=0.1, ins_rate=0.1, del_rate=0.1,
                             fastq=False, uppercase=False)

def test_simulate():
    assert records

def test_read_lens():
    for record in records:
        assert len(record.seq) == read_len

def test_fastq():
    records_mixedcase = sod.simulate_reads('tests/hxb2.fasta', n_reads=1000,
                                           read_len=read_len, sub_rate=0.1,
                                           ins_rate=0.1, del_rate=0.1,
                                           fastq=True, uppercase=False)
    records_uppercase = sod.simulate_reads('tests/hxb2.fasta', n_reads=1000,
                                           read_len=read_len, sub_rate=0.1,
                                           ins_rate=0.1, del_rate=0.1,
                                           fastq=True, uppercase=True)
    quality_mixedcase = records_mixedcase[0].letter_annotations['phred_quality']
    quality_uppercase = records_uppercase[0].letter_annotations['phred_quality']




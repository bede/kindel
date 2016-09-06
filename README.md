# Kindel – indel-aware offline consensus calling for DNA alignments
[![PyPI version](https://badge.fury.io/py/kindel.svg)](https://badge.fury.io/py/kindel)

Kindel is a simple consensus caller which takes a SAM/BAM file as input and generates a *single chromosome* consensus sequence in fasta format. Unlike similar tools it correctly reconciles small indels described in CIGAR fields. At this time no local realignment is performed and so iterative consensus construction and remapping may be necessary to converge on an optimal consensus sequence.

Existing consensus calling approaches are complicated and often involve a variant calling step. While an [elegant and sophisticated streaming approach](https://github.com/karel-brinda/ococo) was recently released, it does not reconcile indels. Kindel generates similar results to the consensus calling tool in [Geneious](http://www.geneious.com/) and doesn't cost $1000.

## Installation
```
pip3 install kindel
```
Dependencies should automatically installed, except for Samtools which is needed for BAM input.

## Usage
### Command line
```
$ kindel alignment.bam --min-depth 1 --threshold-weight 0.5
```
The consensus fasta is sent to `stdout` and a report to `stderr`
```
$ kindel -h
usage: kindel [-h] [-t THRESHOLD_WEIGHT] [-m MIN_DEPTH] bam-path

positional arguments:
  bam-path              path to SAM/BAM file

optional arguments:
  -h, --help            show this help message and exit
  -t THRESHOLD_WEIGHT, --threshold-weight THRESHOLD_WEIGHT
                        consensus threshold weight (default: 0.5)
  -m MIN_DEPTH, --min-depth MIN_DEPTH
                        substitute Ns at coverage depths beneath this value
                        (default: 1)
```

### Python3
```
from kindel import kindel

kindel.bam_to_consensus_seqrecord(bam_path, threshold_weight=0.5, min_depth=1) # returns BioPython SeqRecord
kindel.bam_to_consensus_fasta(bam_path, threshold_weight=0.5, min_depth=1) # returns fasta string
```

## Issues
Please let me know if you run into problems by opening a GitHub issue, [tweeting](https://twitter.com/beconstant) or mailing me via `b at bede dawt im`.
- Intended for use with viral genomes – requires a single chromosome reference sequence
- SAM/BAM files must contain an SQ header line containing the reference sequence length.
- Does not perform local realignment around soft clipped regions – this is a feature for the next release
- Slow (10-20k records/s)

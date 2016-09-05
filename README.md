# Kindel – indel-aware offline consensus calling for DNA alignments
[![PyPI version](https://badge.fury.io/py/kindel.svg)](https://badge.fury.io/py/kindel)

Kindel is a no frills majority rule consensus sequence caller which accepts a SAM/BAM file as input and generates a consensus sequence in fasta format. Unlike other such tools, it correctly reconciles insertion/deletions in the generated reference sequence. It is intended for use with small, single 'chromosome' genomes such as viruses. SAM/BAM files must contain an SQ header line.

Existing open source consensus calling approaches are complicated and typically involve a variant calling step. While an [elegant streaming approach](https://github.com/karel-brinda/ococo) was recently released, it does not reconcile indels. Kindel generates similar results to the consensus caller in the proprietary software [Geneious](http://www.geneious.com/) (and doesn't cost $1000).

## Installation
```
pip3 install kindel
```
Dependencies: `simplesam`, `biopython`, `argh`, `tqdm` and `samtools` (for BAM input)

## Usage
### Command line
```
kindel alignment.bam --min-depth 1 --threshold-weight 0.5
```
The consensus fasta is sent to `stdout` and a report to `stderr`

### Python3
```
from kindel import kindel

kindel.bam_to_consensus_seqrecord(bam_path, threshold_weight=0.5, min_depth=1) # returns BioPython SeqRecord
kindel.bam_to_consensus_fasta(bam_path, fasta_path=sys.stdout, threshold_weight=0.5, min_depth=1)
```

## Issues
Open a GitHub issue, [tweet me](https://twitter.com/beconstant) or mail me via `b at bede dawt im`  
- **Crashes on lots of BAM files at the moment… To be remedied ASAP**
- Intended for use with viral genomes – requires single chromosome reference sequence
- Slow (10-20k records/s)

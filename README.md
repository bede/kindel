# Kindel – indel-aware offline consensus calling for DNA alignments
Kindel is a no frills majority rule consensus sequence caller which accepts a SAM/BAM file as input and generates a consensus sequence in fasta format. Unlike other such tools, it correctly reconciles insertion/deletions in the generated reference sequence. This Python implementation is intended for use with small, **single contig** genomes such as viruses, and processes approximately 15,000 Illumina reads/s in my testing. SAM/BAM files must contain an SQ header line.

Existing open source consensus calling approaches are complicated and typically involve a variant calling step. While an [elegant streaming approach](https://github.com/karel-brinda/ococo) was recently released, it does not reconcile indels. Kindel generates similar results to the consensus caller in the proprietary software [Geneious](http://www.geneious.com/) (and doesn't cost $1000).

## Installation
```
pip3 install kindel
```
Dependencies: `simplesam`, `biopython`, `argh` and `samtools` (for BAM input)

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

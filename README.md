# Kindel – indel-aware offline consensus calling for nucleotide alignments
[![PyPI version](https://badge.fury.io/py/kindel.svg)](https://badge.fury.io/py/kindel)  
[![Build Status](https://travis-ci.org/bede/kindel.svg?branch=master)](https://travis-ci.org/bede/kindel)  


A straightforward consensus caller which accepts a headed SAM/BAM input and generates a majority rule consensus in fasta format,  reconciling small indels described in CIGAR fields so as to maximise read-reference concordance. Furthermore, Kindel permits recovery of consensus sequence across highly divergent regions (such as  those encoding viral envelope proteins) where regions of reads cannot be aligned. With `--realign`, Kindel identifies regions of the reference where soft clipped coverage depth exceeds coverage\*0.5 and attempts to assemble a patched consensus using unaligned sequence context.

Existing consensus calling approaches are complicated and often involve a variant calling step. While an [elegant streaming approach](https://github.com/karel-brinda/ococo) was recently released, it cannot reconcile indels.


## Installation
```python
# Requires Python 3.5+
pip install kindel
```
Dependencies should automatically installed, except for Samtools which is needed for BAM input.


## Usage
### Command line
```
$ kindel consensus --realign alignment.bam
```
The consensus fasta is sent to `stdout` and a report to `stderr`
```
$ kindel -h
usage: kindel [-h] {consensus,weights,variants} ...

positional arguments:
  {consensus,weights,variants}
    consensus           Infer consensus sequence(s) from alignment in SAM/BAM
                        format
    weights             Returns DataFrame of per-site nucleotide frequencies
                        and coverage
    variants            Output variants exceeding specified absolute and
                        relative frequency thresholds

optional arguments:
  -h, --help            show this help message and exit
```
---
```
$  kindel consensus -h
usage: kindel consensus [-h] [-r] [--min-depth MIN_DEPTH]
                        [--min-overlap MIN_OVERLAP] [-c CLIP_DECAY_THRESHOLD]
                        [-t] [-u]
                        bam-path

Infer consensus sequence(s) from alignment in SAM/BAM format

positional arguments:
  bam-path              path to SAM/BAM file

optional arguments:
  -h, --help            show this help message and exit
  -r, --realign         attempt to reconstruct reference around soft-clip
                        boundaries (default: False)
  --min-depth MIN_DEPTH
                        substitute Ns at coverage depths beneath this value
                        (default: 2)
  --min-overlap MIN_OVERLAP
                        match length required to close soft-clipped gaps
                        (default: 7)
  -c CLIP_DECAY_THRESHOLD, --clip-decay-threshold CLIP_DECAY_THRESHOLD
                        read depth fraction at which to cease clip extension
                        (default: 0.1)
  -t, --trim-ends       trim ambiguous nucleotides (Ns) from sequence ends
                        (default: False)
  -u, --uppercase       close gaps using uppercase alphabet (default: False)
```
---
```
$ kindel weights -h
usage: kindel weights [-h] [-r] [-n] bam-path

Returns DataFrame of per-site nucleotide frequencies and coverage

positional arguments:
  bam-path             path to SAM/BAM file

optional arguments:
  -h, --help           show this help message and exit
  -r, --relative       output relative nucleotide frequencies (default: False)
  -n, --no-confidence  skip confidence calculation (default: False)

```
---
```
$ kindel variants -h
usage: kindel variants [-h] [-a ABS_THRESHOLD] [-r REL_THRESHOLD] [-o]
                       bam-path

Output variants exceeding specified absolute and relative frequency thresholds

positional arguments:
  bam-path              path to SAM/BAM file

optional arguments:
  -h, --help            show this help message and exit
  -a ABS_THRESHOLD, --abs-threshold ABS_THRESHOLD
                        absolute frequency (0-∞) threshold above which to call
                        variants (default: 1)
  -r REL_THRESHOLD, --rel-threshold REL_THRESHOLD
                        relative frequency (0.0-1.0) threshold above which to
                        call variants (default: 0.01)
  -o, --only-variants   exclude invariant sites from output (default: False)
```


### Python API
```python
from kindel import kindel

kindel.bam_to_consensus(bam_path, realign=False, min_depth=2, min_overlap=7,
                        clip_decay_threshold=0.1, trim_ends=False, uppercase=False)
```

## Issues
Please let me know if you run into problems by opening a GitHub issue, [tweeting](https://twitter.com/beconstant) or mailing me via `b at bede dawt im`.
- Conceived for use with highly diverse Hepatitis C populations – untested with anything larger
- SAM/BAM files must contain an SQ header line containing the reference sequence length.
- Able to close gaps of up to 2x read length given adequate depth of coverage
- Sometimes requires multiple runs to converge on an optimal consensus
- Slow (10-20k records/s)

## Features
- [x] Reconciliation of CIGAR described insertions and deletions
- [x] Gap closure (`--realign`) using soft-clipped alignment context
- [ ] Support mutiple reference sequences (needs testing)
- [x] Support SAMs from multiple aligners – tested BWA MEM, Segemehl
- [x] Frequency based variant calling with `kindel variants` (no VCF output)
- [x] Plotting of clip frequencies
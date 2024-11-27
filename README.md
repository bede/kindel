# Kindel: indel-aware consensus from aligned BAM

[![JOSS status](http://joss.theoj.org/papers/117efd1fc35bb2011311f73d3fa0b545/status.svg)](http://joss.theoj.org/papers/117efd1fc35bb2011311f73d3fa0b545)  
[![PyPI version](https://badge.fury.io/py/kindel.svg)](https://badge.fury.io/py/kindel)  
[![Build Status](https://travis-ci.org/bede/kindel.svg?branch=master)](https://travis-ci.org/bede/kindel)  



Kindel reconciles substitutions and CIGAR-described indels to to produce a majority consensus from a SAM/BAM file. Using the `--realign` option, Kindel can optionally recover consensus across short alignment gaps using soft-clipped sequence information. Where Kindel finds 'clip-dominant' regions of an alignment, in realignment mode it attempts to reassemble the consensus sequence using unaligned sequence context. Primarily intended for use with small alignments of e.g. virus genomes, it has been tested with BAMs created by aligners BWA and Minimap2. If you encounter problems, please open an issue. Please also cite the [JOSS article](http://joss.theoj.org/papers/117efd1fc35bb2011311f73d3fa0b545) if you find this useful.



### Core functionality

![clip-dominant region](kindelflow.png)



### Reassembly of clip-dominant regions (CDRs) with `--realign`

![clip-dominant region](cdrs.png)


## Features

- Consensus of aligned substititutions, insertions and deletions
- Gap closure (`--realign`) using overlapping soft-clipped alignment context
- Tested with Illumina alignments from BWA, Minimap2 and Segemehl 
- Support for BAMs with multiple reference contigs, chromosomes
- Crude frequency-based variant calling with `kindel variants` (no VCF output)



## Limitations

- Intended for use with small alignments of e.g. virus genomes. Expect slow performance with megabase genomes.
- SAM/BAM files must contain an SQ header line with reference sequence(s) length
- Realignment mode (`--realign`) is able to close gaps of up to 2x read length given ample depth of coverage



## Installation

```shell
# Requires Python 3.9+ and Samtools
pip install kindel
```
For a complete installation using a conda-compatible package manager:

```
conda create -y -n kindel python=3.13 samtools
conda activate kindel
pip install kindel
```

For a local development install:

```
pip install --editable '.[dev]'  # pip install kindel '.[dev]'
```



## Usage

Also see [`usage.ipynb`](usage.ipynb)



### Command line

```bash
$ kindel consensus alignment.bam > cns.fa
```
Generate a consensus sequence from an aligned BAM, saving the consensus sequence to `cns.fa`



```bash
$ kindel consensus --realign alignment.bam > cns.fa
```

Generate a consensus sequence from an aligned BAM with realignment mode enabled, allowing closure of small gaps in the consensus sequence



```bash
$ kindel plot alignment.bam
```

Generate an interactive plot showing aligned depth alongside insertion, deletion and soft clipping frequency across the genome



```bash
$ kindel -h
usage: kindel [-h] {consensus,weights,features,variants,plot} ...

positional arguments:
  {consensus,weights,features,variants,plot,version}
    consensus           Infer consensus sequence(s) from alignment in SAM/BAM
                        format
    weights             Returns table of per-site nucleotide frequencies and
                        coverage
    features            Returns table of per-site nucleotide frequencies and
                        coverage including indels
    variants            Output variants exceeding specified absolute and
                        relative frequency thresholds
    plot                Plot sitewise soft clipping frequency across reference
                        and genome
    version             Show version

optional arguments:
  -h, --help            show this help message and exit

```
---
```bash
$  kindel consensus -h
usage: kindel consensus [-h] [-r] [--min-depth MIN_DEPTH]
                        [--min-overlap MIN_OVERLAP] [-c CLIP_DECAY_THRESHOLD]
                        [--mask-ends MASK_ENDS] [-t] [-u]
                        bam_path

Infer consensus sequence(s) from alignment in SAM/BAM format

positional arguments:
  bam_path              path to SAM/BAM file

optional arguments:
  -h, --help            show this help message and exit
  -r, --realign         attempt to reconstruct reference around soft-clip
                        boundaries (default: False)
  --min-depth MIN_DEPTH
                        substitute Ns at coverage depths beneath this value
                        (default: 1)
  --min-overlap MIN_OVERLAP
                        match length required to close soft-clipped gaps
                        (default: 7)
  -c CLIP_DECAY_THRESHOLD, --clip-decay-threshold CLIP_DECAY_THRESHOLD
                        read depth fraction at which to cease clip extension
                        (default: 0.1)
  --mask-ends MASK_ENDS
                        ignore clip dominant positions within n positions of
                        termini (default: 50)
  -t, --trim-ends       trim ambiguous nucleotides (Ns) from sequence ends
                        (default: False)
  -u, --uppercase       close gaps using uppercase alphabet (default: False)
```



### Python API

```python
from kindel import kindel

kindel.bam_to_consensus(bam_path, realign=False, min_depth=2, min_overlap=7,
                        clip_decay_threshold=0.1, trim_ends=False, uppercase=False)
```



## Issues

Please let me know if you run into problems by opening a GitHub issue, tweeting [@beconstant](https://twitter.com/beconstant) or mailing me via `b at bede dawt im`. Ideally send me your BAM, or a subsample of it!



## Contributing

If you would like to contribute to this project, please open an issue or contact the author directly using the details above. Please note that this project is released with a [Contributor Code of Conduct](https://github.com/statsmaths/kerasR/blob/master/CONDUCT.md), and by participating in this project you agree to abide by its terms.

Before opening a pull request, please:

- Ensure tests pass by executing `pytest` inside the package directory  (requires `pytest` package)
- Increment the version number inside `__init__.py` according to [SemVer](http://semver.org/)
- Update documentation and/or tests if possible
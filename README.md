# Kindel: indel-aware consensus

[![JOSS status](http://joss.theoj.org/papers/117efd1fc35bb2011311f73d3fa0b545/status.svg)](http://joss.theoj.org/papers/117efd1fc35bb2011311f73d3fa0b545)  
[![PyPI version](https://badge.fury.io/py/kindel.svg)](https://badge.fury.io/py/kindel)  
[![Build Status](https://travis-ci.org/bede/kindel.svg?branch=master)](https://travis-ci.org/bede/kindel)  



Kindel reconciles substitutions and CIGAR-described indels to to produce a majority consensus from a SAM or BAM file. Kindel can optionally further recover consensus across unaligned regions (such as those frequently seen in RNA virus populations) using soft-clipped sequence information. With **`--realign`**, Kindel identifies regions of the reference that are clip-dominant (>depth\*0.5) and attempts to assemble a patched consensus using unaligned sequence context. Intended for use with small virus genomes and Illumina data, and tested against aligners BWA, Minimap2 and Segemehl. Developed after growing frustrated with dubious and overly complicated existing methods. I may improve performance for bacterial genomes in a future update. Update: [Pilon](https://github.com/broadinstitute/pilon) also now solves similar problems rather well. **If your BAM generates errors, please open an issue and I'll try to fix it. Please also cite the [JOSS article](http://joss.theoj.org/papers/117efd1fc35bb2011311f73d3fa0b545) if you find this useful.**





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





### Todo

- [ ] Customisable clipping threshold weight
- [ ] Improve performance with bacterial genomes and long reads
- [ ] Optionally fill gaps using a supplied reference
- [ ] Parallelisation





## Limitations

- Slow (10-20k records/s), & will probably explode with bacterial genomes
- SAM/BAM files must contain an SQ header line with reference sequence(s) length
- Able to close gaps of up to 2x read length given adequate depth of coverage
- May require multiple runs to converge on a consensus





## Installation

```python
# Requires Python 3.6+
pip install kindel
```
Dependencies should automatically installed, except for Samtools which is needed for BAM input.



## Usage

Also see [`usage.ipynb`](usage.ipynb)



### Command line

```bash
$ kindel consensus alignment.bam > cns.fa
```
The consensus fasta is sent to `stdout` and a report to `stderr`
```bash
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
```bash
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
```bash
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
```bash
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

Please let me know if you run into problems by opening a GitHub issue, tweeting [@beconstant](https://twitter.com/beconstant) or mailing me via `b at bede dawt im`. Ideally send me your BAM, or a subsample of it!



## Contributing

If you would like to contribute to this project, please open an issue or contact the author directly using the details above. Please note that this project is released with a [Contributor Code of Conduct](https://github.com/statsmaths/kerasR/blob/master/CONDUCT.md), and by participating in this project you agree to abide by its terms.

Before issuing a pull request, please:

- Ensure tests pass by executing `pytest` inside the package directory  (requires `pytest` package)
- Increment the version number inside `__init__.py` according to [SemVer](http://semver.org/)
- Update documentation and/or tests if possible
import sys
import argh
import pandas as pd

from kindel import kindel


def consensus(bam_path: 'path to SAM/BAM file',
              fix_gaps: 'attempt to reconcile reference at soft-clip boundaries'=False,
              trim_ends: 'trim ambiguous nucleotides (Ns) from sequence ends'=False,
              threshold: 'consensus threshold weight'=0.5,
              min_depth: 'substitute Ns at coverage depths beneath this value'=2,
              closure_k: 'match length required to close soft-clipped gaps'=7,
              uppercase: 'close gaps using uppercase alphabet'=False):
    '''Infer consensus sequence(s) from alignment in SAM/BAM format'''
    result = kindel.bam_to_consensus(bam_path,
                                     fix_gaps,
                                     trim_ends,
                                     threshold,
                                     min_depth,
                                     closure_k)
    print(result.report, file=sys.stderr)
    return SeqIO.write(result.consensuses, sys.stdout,'fasta')


def weights(bam_path: 'path to SAM/BAM file',
            relative: 'output relative nucleotide frequencies'=False):
    '''Returns DataFrame of per-site nucleotide frequencies and coverage'''
    weights_df = kindel.weights(bam_path, relative)
    return weights_df.to_csv(sys.stdout, sep='\t', index=False)


def variants(bam_path: 'path to SAM/BAM file',
             abs_threshold: 'absolute frequency (0-∞) threshold above which to call variants'=1,
             rel_threshold: 'relative frequency (0.0-1.0) threshold above which to call variants'=0.01,
             consensus_threshold: 'consensus threshold weight'=0.5):
    '''Output variants exceeding specified absolute and relative frequency thresholds'''
    return kindel.variants(bam_path, abs_threshold, rel_threshold, consensus_threshold)

def main():
    parser = argh.ArghParser()
    parser.add_commands([consensus, weights, variants])
    parser.dispatch()


if __name__ == '__main__':
    main()

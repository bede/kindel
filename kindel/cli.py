import sys

import argh

from kindel import kindel


def main(bam_path, fasta_path=sys.stdout, threshold_weight=0.5, min_depth=1):
    kindel.bam_to_consensus_fasta(bam_path,
                                  fasta_path,
                                  threshold_weight,
                                  min_depth)

argh.dispatch_command(main)
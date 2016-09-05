import argh

from kindel import kindel


def main():
    argh.dispatch_command(kindel.bam_to_consensus_fasta)


if __name__ == '__main__':
    main()
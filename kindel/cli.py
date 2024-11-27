import sys
import argh

from kindel import kindel

from kindel import __version__


def consensus(
    bam_path: "path to SAM/BAM file",
    realign: "attempt to reconstruct reference around soft-clip boundaries" = False,
    min_depth: "substitute Ns at coverage depths beneath this value" = 1,
    min_overlap: "match length required to close soft-clipped gaps" = 7,
    clip_decay_threshold: "read depth fraction at which to cease clip extension" = 0.1,
    mask_ends: "ignore clip dominant positions within n positions of termini" = 50,
    trim_ends: "trim ambiguous nucleotides (Ns) from sequence ends" = False,
    uppercase: "close gaps using uppercase alphabet" = False,
):
    """Infer consensus sequence(s) from alignment in SAM/BAM format"""
    result = kindel.bam_to_consensus(
        bam_path,
        realign,
        min_depth,
        min_overlap,
        clip_decay_threshold,
        mask_ends,
        trim_ends,
        uppercase,
    )
    print("\n".join([r for r in result.refs_reports.values()]), file=sys.stderr)
    for consensus_record in result.consensuses:
        print(f">{consensus_record.name}")
        print(consensus_record.sequence)


def weights(
    bam_path: "path to SAM/BAM file",
    relative: "output relative nucleotide frequencies" = False,
    no_confidence: "skip confidence calculation" = False,
):
    """Returns table of per-site nucleotide frequencies and coverage"""
    weights_df = kindel.weights(bam_path, relative, no_confidence)
    weights_df.to_csv(sys.stdout, sep="\t", index=False)


def features(bam_path: "path to SAM/BAM file"):
    """Returns table of per-site nucleotide frequencies and coverage including indels"""
    weights_df = kindel.features(bam_path)
    weights_df.to_csv(sys.stdout, sep="\t", index=False)


def variants(
    bam_path: "path to SAM/BAM file",
    abs_threshold: "absolute frequency (0-∞) threshold above which to call variants" = 1,
    rel_threshold: "relative frequency (0.0-1.0) threshold above which to call variants" = 0.01,
    only_variants: "exclude invariant sites from output" = False,
    absolute: "report absolute variant frequencies" = False,
):
    """Output variants exceeding specified absolute and relative frequency thresholds"""
    variants_df = kindel.variants(
        bam_path, abs_threshold, rel_threshold, only_variants, absolute
    )
    variants_df.to_csv(sys.stdout, sep="\t", index=False, na_rep=0)


def plot(bam_path: "path to SAM/BAM file"):
    """Plot sitewise soft clipping frequency across reference and genome"""
    return kindel.plotly_clips(bam_path)


def version():
    """Show version"""
    return f"kindel {__version__}"


def main():
    parser = argh.ArghParser()
    parser.add_commands([consensus, weights, features, variants, plot, version])
    parser.dispatch()


if __name__ == "__main__":
    main()

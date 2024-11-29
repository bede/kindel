import os
import tqdm
import simplesam
import scipy.stats

import logging
from collections import OrderedDict, defaultdict, namedtuple

import dnaio

import numpy as np
import pandas as pd

from . import cli

# logging.basicConfig(level=logging.DEBUG)

Region = namedtuple("Region", ["start", "end", "seq", "direction"])


def parse_records(ref_id, ref_len, records):
    """
    Iterate over records, returning namedtuple of base frequencies, indels and soft clipping info
    weights: lists of dicts of base frequencies after CIGAR reconciliation
    insertions, deletions: list of dicts of base frequencies
    clip_starts, clip_ends: lists of clipping start and end frequencies in an LTR direction
    clip_start_weights, clip_end_weights: base frequencies from left- and right-clipped regions
    """
    weights = [{"A": 0, "T": 0, "G": 0, "C": 0, "N": 0} for p in range(ref_len)]
    clip_start_weights = [
        {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0} for p in range(ref_len)
    ]
    clip_end_weights = [
        {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0} for p in range(ref_len)
    ]
    clip_starts = [0] * (ref_len + 1)
    clip_ends = [0] * (ref_len + 1)
    insertions = [defaultdict(int) for p in range(ref_len + 1)]
    deletions = [0] * (ref_len + 1)
    for record in tqdm.tqdm(records, desc="loading sequences"):  # Progress bar
        q_pos = 0
        r_pos = record.pos - 1  # Zero indexed genome coordinates
        if (
            not record.mapped or len(record.seq) <= 1
        ):  # Skips unmapped, reads with sequence '*'
            continue
        for i, cigarette in enumerate(record.cigars):  # StopIteration -> RuntimeError
            length, operation = cigarette
            if operation in {"M", "=", "X"}:  # Catch SAM 1.{3,4} matches and subs
                for _ in range(length):
                    q_nt = record.seq[q_pos].upper()
                    weights[r_pos][q_nt] += 1
                    r_pos += 1
                    q_pos += 1
            elif operation == "I":
                nts = record.seq[q_pos : q_pos + length].upper()
                insertions[r_pos][nts] += 1
                q_pos += length
            elif operation == "D":
                for del_i in range(length):
                    deletions[r_pos + del_i] += 1
                r_pos += length
            elif operation == "S":
                if (
                    i == 0
                ):  # Count right-of-gap / l-clipped start positions (e.g. start of ref)
                    clip_ends[r_pos] += 1
                    for gap_i in range(length):
                        q_nt = record.seq[gap_i].upper()
                        rel_r_pos = r_pos - length + gap_i
                        if rel_r_pos >= 0:
                            clip_end_weights[rel_r_pos][q_nt] += 1
                    q_pos += length
                else:  # Count left-of-gap / r-clipped start position (e.g. end of ref)
                    clip_starts[r_pos] += 1
                    for pos in range(length):
                        q_nt = record.seq[q_pos].upper()
                        if r_pos < ref_len:
                            clip_start_weights[r_pos][q_nt] += 1
                            r_pos += 1
                            q_pos += 1

    aligned_depth = [sum(w.values()) for w in weights]
    weights_consensus_seq = "".join([consensus(w)[0] for w in weights])
    discordant_depth = [
        sum({nt: w[nt] for nt in [k for k in w.keys() if k != cns_nt]}.values())
        for w, cns_nt in zip(weights, weights_consensus_seq)
    ]
    consensus_depth = np.array(aligned_depth) - np.array(discordant_depth)
    clip_start_depth = [
        sum({nt: w[nt] for nt in list("ACGT")}.values()) for w in clip_start_weights
    ]
    clip_end_depth = [
        sum({nt: w[nt] for nt in list("ACGT")}.values()) for w in clip_end_weights
    ]
    clip_depth = list(map(lambda x, y: x + y, clip_start_depth, clip_end_depth))
    alignment = namedtuple(
        "alignment",
        [
            "ref_id",
            "weights",
            "insertions",
            "deletions",
            "clip_starts",
            "clip_ends",
            "clip_start_weights",
            "clip_end_weights",
            "clip_start_depth",
            "clip_end_depth",
            "clip_depth",
            "consensus_depth",
        ],
    )

    return alignment(
        ref_id,
        weights,
        insertions,
        deletions,
        clip_starts,
        clip_ends,
        clip_start_weights,
        clip_end_weights,
        clip_start_depth,
        clip_end_depth,
        clip_depth,
        consensus_depth,
    )


def parse_bam(bam_path):
    """
    Returns alignment information for each reference sequence as an OrderedDict
    """
    alignments = OrderedDict()
    with open(bam_path, "r") as bam_fh:
        bam = simplesam.Reader(bam_fh)
        refs_lens = {
            name.replace("SN:", ""): int(length[0].replace("LN:", ""))
            for name, length in bam.header["@SQ"].items()
        }

        refs_records = defaultdict(list)
        for r in bam:
            for id in refs_lens:
                refs_records[r.rname].append(r)

        if "*" in refs_records:
            del refs_records["*"]

    # assert len(refs_records) <= 1, 'Detected primary mappings to more than one reference'
    # Use samtools view to extract single contig primary mappings
    # Otherwise would make a useful enhancement

    for ref_id, records in refs_records.items():
        alignments[ref_id] = parse_records(ref_id, refs_lens[ref_id], records)

    return alignments


def cdr_start_consensuses(
    weights,
    deletions,
    clip_start_weights,
    clip_start_depth,
    clip_decay_threshold,
    mask_ends,
):
    """
    Returns list of Region instances for right clipped consensuses of clip-dominant region
    """
    positions = list(range(len(weights)))
    masked_positions = positions[:mask_ends] + positions[-mask_ends:]
    regions = []
    logging.debug("cdr_start_consensuses()")
    for (
        pos,
        csd,
        w,
        dp,
    ) in zip(positions, clip_start_depth, weights, deletions):
        cdr_positions = []
        for r in regions:
            range_list = range(r.start, r.end)
            for s in range_list:
                cdr_positions.append(s)
        if (
            csd / (sum(w.values()) + dp + 1) > 0.5
            and pos not in masked_positions + cdr_positions
        ):
            logging.debug("Starting extension")
            start_pos = pos
            clip_consensus = ""
            for pos_, (csd_, sw_, w_, dd_) in enumerate(
                zip(
                    clip_start_depth[pos:],
                    clip_start_weights[pos:],
                    weights[pos:],
                    deletions[pos:],
                )
            ):
                logging.debug(
                    f"start_pos: {start_pos}, pos: {pos_}, csd: {csd_}, sum (weights): {sum(w_.values())}, sum (deletions): {dd_}"
                )
                if csd_ > sum(w_.values(), dd_) * clip_decay_threshold:
                    clip_consensus += consensus(sw_)[0]
                else:
                    end_pos = start_pos + pos_
                    logging.debug("Stopping extension")
                    break
            regions.append(Region(start_pos, end_pos, clip_consensus, "→"))

    for region in regions:
        logging.debug(region)
    return regions


def cdr_end_consensuses(
    weights,
    deletions,
    clip_end_weights,
    clip_end_depth,
    clip_decay_threshold,
    mask_ends,
):
    """
    Returns list of Region instances for left clipped consensuses of clip-dominant region
    """
    positions = list(range(len(weights)))
    masked_positions = positions[:mask_ends] + positions[-mask_ends:]
    reversed_weights = list(
        sorted(
            zip(positions, clip_end_depth, clip_end_weights, weights, deletions),
            key=lambda x: x[0],
            reverse=True,
        )
    )
    regions = []
    logging.debug("cdr_end_consensuses()")
    for pos, ced, _, w, dp in reversed_weights:
        cdr_positions = []
        for r in regions:
            for s in range(r.start, r.end):
                cdr_positions.append(s)
        if (
            ced / (sum(w.values()) + dp + 1) > 0.5
            and pos not in masked_positions + cdr_positions
        ):
            logging.debug("Starting extension")
            end_pos = pos + 1  # Start with end since we're iterating in reverse
            rev_clip_consensus = None
            for pos_, ced_, cew_, w_, dd_ in reversed_weights[len(positions) - pos :]:
                logging.debug(
                    f"end_pos: {end_pos}, pos: {pos_}, ced: {ced_}, sum (weights): {sum(w_.values())}, sum (deletions): {dd_}"
                )
                if ced_ > sum(w_.values(), dd_) * clip_decay_threshold:
                    if (
                        not rev_clip_consensus
                    ):  # Add first base to account for lag in clip coverage
                        rev_clip_consensus = consensus(clip_end_weights[pos_ + 1])[0]
                    rev_clip_consensus += consensus(cew_)[0]
                else:
                    start_pos = pos_
                    clip_consensus = rev_clip_consensus[::-1]
                    logging.debug("Stopping extension")
                    break
            regions.append(Region(start_pos, end_pos, clip_consensus, "←"))

    for region in regions:
        logging.debug(region)
    return regions


def cdrp_consensuses(
    weights,
    deletions,
    clip_start_weights,
    clip_end_weights,
    clip_start_depth,
    clip_end_depth,
    clip_decay_threshold,
    mask_ends,
):
    """
    Returns list of 2-tuples of L&R clipped consensus sequences around clip-dominant regions
    Pairs overlapping right- (→) and left- (←) clipped sequences around CDRs
    """
    combined_cdrs = cdr_start_consensuses(
        weights,
        deletions,
        clip_start_weights,
        clip_start_depth,
        clip_decay_threshold,
        mask_ends,
    ) + cdr_end_consensuses(
        weights,
        deletions,
        clip_end_weights,
        clip_end_depth,
        clip_decay_threshold,
        mask_ends,
    )
    paired_cdrs = []
    fwd_cdrs = [r for r in combined_cdrs if r.direction == "→"]
    rev_cdrs = [r for r in combined_cdrs if r.direction == "←"]
    for fwd_cdr in fwd_cdrs:
        fwd_cdr_range = range(fwd_cdr.start, fwd_cdr.end)
        for rev_cdr in rev_cdrs:
            rev_cdr_range = range(rev_cdr.start, rev_cdr.end)
            if set(fwd_cdr_range).intersection(rev_cdr_range):
                paired_cdrs.append((fwd_cdr, rev_cdr))
                break

    logging.debug("cdrp_consensuses()")

    return paired_cdrs


def merge_by_lcs(s1, s2, min_overlap):
    """Returns superstring of s1 and s2 about an exact overlap of len > min_overlap"""

    def lcs(s1, s2):
        m = [[0] * (1 + len(s2)) for i in range(1 + len(s1))]
        longest, x_longest = 0, 0
        for x in range(1, 1 + len(s1)):
            for y in range(1, 1 + len(s2)):
                if s1[x - 1] == s2[y - 1]:
                    m[x][y] = m[x - 1][y - 1] + 1
                    if m[x][y] > longest:
                        longest = m[x][y]
                        x_longest = x
                else:
                    m[x][y] = 0
        return s1[x_longest - longest : x_longest]

    logging.debug(f"merge_by_lcs(): s1: {s1}, s2: {s2}, min_overlap: {min_overlap}")
    longest_common_subsequence = lcs(s1, s2)
    if len(longest_common_subsequence) < min_overlap:
        return None  # Failed
    left_part = s1.split(longest_common_subsequence, 1)[0]
    right_part = s2.split(longest_common_subsequence, 1)[1]
    merged_sequence = left_part + longest_common_subsequence + right_part
    return merged_sequence


def merge_cdrps(cdrps, min_overlap):
    """Returns merged clip-dominant region pairs as Region instances"""
    merged_cdrps = []
    for cdrp in cdrps:
        fwd_cdr, rev_cdr = cdrp
        merged_seq = merge_by_lcs(
            fwd_cdr.seq, rev_cdr.seq, min_overlap
        )  # Fails as None
        if not merged_seq:
            logging.warning(
                f"No overlap found for CDRP spanning positions {fwd_cdr.start}-{rev_cdr.end} (min_overlap = {min_overlap})"
            )
        merged_cdrps.append(Region(fwd_cdr.start, rev_cdr.end, merged_seq, None))

    logging.debug("Merged CDRPs")
    logging.debug(merged_cdrps)
    return merged_cdrps


def consensus(weight):
    """
    Returns tuple of consensus base, weight and flag indicating a tie for consensus
    """
    base, frequency = (
        max(weight.items(), key=lambda x: x[1]) if sum(weight.values()) else ("N", 0)
    )
    weight_sans_consensus = {k: d for k, d in weight.items() if k != base}
    tie = True if frequency and frequency in weight_sans_consensus.values() else False
    aligned_depth = sum(weight.values())
    proportion = round(frequency / aligned_depth, 2) if aligned_depth else 0

    return (base, frequency, proportion, tie)


def s_overhang_consensus(clip_start_weights, start_pos, min_depth, max_len=1000):
    """
    Returns consensus sequence (string) of clipped reads at specified position
    start_pos is the first position described by the CIGAR-S
    """
    consensus_overhang = ""
    for pos in range(start_pos, start_pos + max_len):
        pos_consensus = consensus(clip_start_weights[pos])
        if pos_consensus[1] >= min_depth:
            consensus_overhang += pos_consensus[0]
        else:
            break

    return consensus_overhang


def e_overhang_consensus(clip_end_weights, start_pos, min_depth, max_len=1000):
    """
    Returns consensus sequence (string) of clipped reads at specified position
    end_pos is the last position described by the CIGAR-S
    """
    rev_consensus_overhang = ""
    for pos in range(start_pos, start_pos - max_len, -1):
        pos_consensus = consensus(clip_end_weights[pos])
        if pos_consensus[1] >= min_depth:
            rev_consensus_overhang += pos_consensus[0]
        else:
            break
    consensus_overhang = rev_consensus_overhang[::-1]

    return consensus_overhang


def s_flanking_seq(start_pos, weights, min_depth, k):
    """
    Returns consensus sequence (string) flanking LHS of soft-clipped gaps
    """
    flank_seq = ""
    for pos in range(start_pos - k, start_pos):
        pos_consensus = consensus(weights[pos])
        if pos_consensus[1] >= min_depth:
            flank_seq += pos_consensus[0]
    return flank_seq


def e_flanking_seq(end_pos, weights, min_depth, k):
    """
    Returns consensus sequence (string) flanking RHS of soft-clipped gaps
    """
    flank_seq = ""
    for pos in range(end_pos + 1, end_pos + k + 1):
        pos_consensus = consensus(weights[pos])
        if pos_consensus[1] >= min_depth:
            flank_seq += pos_consensus[0]

    return flank_seq


def consensus_sequence(
    weights, insertions, deletions, cdr_patches, trim_ends, min_depth, uppercase
):
    consensus_seq = ""
    changes = [None] * len(weights)
    skip_positions = 0
    for pos, weight in tqdm.tqdm(
        enumerate(weights), total=len(weights), desc="building consensus"
    ):
        if skip_positions:
            skip_positions -= 1
            continue
        if cdr_patches and any(r.start == pos and r.seq for r in cdr_patches):
            cdr_patch = next(r for r in cdr_patches if r.start == pos)
            consensus_seq += cdr_patch.seq.lower()
            cdr_span = cdr_patch.end - cdr_patch.start
            skip_positions += cdr_span - 1
            continue
        ins_freq = sum(insertions[pos].values()) if insertions[pos] else 0
        del_freq = deletions[pos]
        aligned_depth = sum({nt: weight[nt] for nt in list("ACGT")}.values())
        try:
            aligned_depth_next = sum(
                {nt: weights[pos + 1][nt] for nt in list("ACGT")}.values()
            )
        except IndexError:
            aligned_depth_next = 0
        threshold_freq = aligned_depth * 0.5
        indel_threshold_freq = min(threshold_freq, aligned_depth_next * 0.5)
        if del_freq > threshold_freq:
            changes[pos] = "D"
        elif aligned_depth < min_depth:
            consensus_seq += "N"
            changes[pos] = "N"
        else:
            if ins_freq > indel_threshold_freq:
                insertion = consensus(insertions[pos])
                consensus_seq += insertion[0].lower() if not insertion[3] else "N"
                changes[pos] = "I"
            pos_consensus = consensus(weight)
            consensus_seq += pos_consensus[0] if not pos_consensus[3] else "N"
    if trim_ends:
        consensus_seq = consensus_seq.strip("N")
    if uppercase:
        consensus_seq = consensus_seq.upper()

    return consensus_seq, changes


def consensus_seqrecord(consensus, ref_id):
    return dnaio.Sequence(name=f"{ref_id}_cns", sequence=consensus, qualities=None)


def build_report(
    ref_id,
    weights,
    changes,
    cdr_patches,
    bam_path,
    realign,
    min_depth,
    min_overlap,
    clip_decay_threshold,
    trim_ends,
    uppercase,
):
    aligned_depth = [sum({nt: w[nt] for nt in list("ACGT")}.values()) for w in weights]
    ambiguous_sites = []
    insertion_sites = []
    deletion_sites = []
    cdr_patches_fmt = (
        ["{}-{}: {}".format(r.start, r.end, r.seq) for r in cdr_patches]
        if cdr_patches
        else ""
    )
    for pos, change in enumerate(changes):
        if change == "N":
            ambiguous_sites.append(str(pos))
        elif change == "I":
            insertion_sites.append(str(pos + 1))
        elif change == "D":
            deletion_sites.append(str(pos))
    report = "========================= REPORT ===========================\n"
    report += "reference: {}\n".format(ref_id)
    report += "options:\n"
    report += "- bam_path: {}\n".format(bam_path)
    report += "- min_depth: {}\n".format(min_depth)
    report += "- realign: {}\n".format(realign)
    report += "    - min_overlap: {}\n".format(min_overlap)
    report += "    - clip_decay_threshold: {}\n".format(clip_decay_threshold)
    report += "- trim_ends: {}\n".format(trim_ends)
    report += "- uppercase: {}\n".format(uppercase)
    report += "observations:\n"
    report += "- min, max observed depth: {}, {}\n".format(
        min(aligned_depth), max(aligned_depth)
    )
    report += "- ambiguous sites: {}\n".format(", ".join(ambiguous_sites))
    report += "- insertion sites: {}\n".format(", ".join(insertion_sites))
    report += "- deletion sites: {}\n".format(", ".join(deletion_sites))
    report += "- clip-dominant regions: {}\n".format(", ".join(cdr_patches_fmt))

    return report


def bam_to_consensus(
    bam_path,
    realign=False,
    min_depth=1,
    min_overlap=7,
    clip_decay_threshold=0.1,
    mask_ends=50,
    trim_ends=False,
    uppercase=False,
):
    consensuses = []
    refs_changes = {}
    refs_reports = {}
    for ref_id, aln in parse_bam(bam_path).items():
        if realign:
            cdrps = cdrp_consensuses(
                aln.weights,
                aln.deletions,
                aln.clip_start_weights,
                aln.clip_end_weights,
                aln.clip_start_depth,
                aln.clip_end_depth,
                clip_decay_threshold,
                mask_ends,
            )
            cdr_patches = merge_cdrps(cdrps, min_overlap)

        else:
            cdr_patches = None
        # logging.debug(aln.weights, aln.clip_start_weights, aln.clip_end_weights, aln.insertions,  aln.deletions, cdr_patches, trim_ends, min_depth, uppercase)
        consensus, changes = consensus_sequence(
            aln.weights,
            aln.insertions,
            aln.deletions,
            cdr_patches,
            trim_ends,
            min_depth,
            uppercase,
        )
        report = build_report(
            ref_id,
            aln.weights,
            changes,
            cdr_patches,
            bam_path,
            realign,
            min_depth,
            min_overlap,
            clip_decay_threshold,
            trim_ends,
            uppercase,
        )
        consensuses.append(consensus_seqrecord(consensus, ref_id))
        refs_reports[ref_id] = report
        refs_changes[ref_id] = changes
    result = namedtuple("result", ["consensuses", "refs_changes", "refs_reports"])

    return result(consensuses, refs_changes, refs_reports)


def weights(
    bam_path: "path to SAM/BAM file",
    relative: "output relative nucleotide frequencies" = False,
    confidence: "calculate confidence interval" = True,
    confidence_alpha: "confidence interval alpha" = 0.01,
):
    """
    Returns DataFrame of per-site nucleotide frequencies, depth, consensus, clip start
    clip end, consensus, confidence intervals, and entropy
    """

    def binomial_ci(count, nobs, alpha=0.01):
        """Returns lower, upper bounds of the Jeffrey binomial proportion confidence interval"""
        lower_ci, upper_ci = scipy.stats.beta.interval(
            1 - alpha, count + 0.5, nobs - count + 0.5
        )
        return lower_ci, upper_ci

    refs_alns = parse_bam(bam_path)
    weights_fmt = []
    for chrom, aln in refs_alns.items():
        for i, w in enumerate(aln.weights, start=1):
            weight_dict = dict(w, chrom=chrom, pos=i)
            weight_dict["insertions"] = sum(aln.insertions[i - 1].values())
            weight_dict["deletions"] = aln.deletions[i - 1]
            weight_dict["clip_starts"] = aln.clip_starts[i - 1]
            weight_dict["clip_ends"] = aln.clip_ends[i - 1]
            weights_fmt.append(weight_dict)

    weights_df = pd.DataFrame(
        weights_fmt,
        columns=[
            "chrom",
            "pos",
            "A",
            "C",
            "G",
            "T",
            "N",
            "insertions",
            "deletions",
            "clip_starts",
            "clip_ends",
        ],
    )
    weights_df["depth"] = weights_df[["A", "C", "G", "T", "N", "deletions"]].sum(axis=1)
    consensus_depths_df = weights_df[["A", "C", "G", "T", "N", "deletions"]].max(axis=1)
    weights_df["consensus"] = consensus_depths_df.divide(weights_df.depth)

    rel_weights_df = pd.DataFrame()
    for nt in ["A", "C", "G", "T", "N", "deletions"]:
        rel_weights_df[[nt]] = weights_df[[nt]].divide(weights_df.depth, axis=0)
        rel_weights_df = rel_weights_df.round(
            {"A": 4, "C": 4, "G": 4, "T": 4, "N": 4, "deletions": 4}
        )

    weights_df["shannon"] = [
        scipy.stats.entropy(x) for x in rel_weights_df[["A", "C", "G", "T"]].values
    ]

    if confidence:
        conf_ints = [
            binomial_ci(c, t, confidence_alpha)
            for c, t in zip(consensus_depths_df, weights_df["depth"])
        ]
        weights_df["lower_ci"] = [ci[0] for ci in conf_ints]
        weights_df["upper_ci"] = [ci[1] for ci in conf_ints]

    if relative:
        for nt in ["A", "C", "G", "T", "N"]:
            weights_df[[nt]] = rel_weights_df[[nt]]

    return weights_df.round(dict(consensus=3, lower_ci=3, upper_ci=3, shannon=3))


def features(bam_path: "path to SAM/BAM file"):
    """
    Returns DataFrame of relative per-site nucleotide frequencies, insertions, deletions
    and entropy
    """
    refs_alns = parse_bam(bam_path)
    weights_fmt = []
    for ref, aln in refs_alns.items():
        weights_fmt.extend(
            [dict(w, chrom=ref, pos=i) for i, w in enumerate(aln.weights, start=1)]
        )
    for pos, weight in enumerate(weights_fmt):
        weight["i"] = sum(aln.insertions[pos].values())
        weight["d"] = aln.deletions[pos]

    # Think about which columns should sum to 1
    weights_df = pd.DataFrame(
        weights_fmt, columns=["chrom", "pos", "A", "C", "G", "T", "N", "i", "d"]
    )
    weights_df["depth"] = weights_df[["A", "C", "G", "T", "N", "d"]].sum(axis=1)
    consensus_depths_df = weights_df[["A", "C", "G", "T", "N"]].max(axis=1)
    weights_df["consensus"] = consensus_depths_df.divide(weights_df.depth)

    for nt in ["A", "C", "G", "T", "N", "i", "d"]:
        weights_df[[nt]] = weights_df[[nt]].divide(weights_df.depth, axis=0)

    weights_df["shannon"] = [
        scipy.stats.entropy(x)
        for x in weights_df[["A", "C", "G", "T", "i", "d"]].values
    ]

    return weights_df.round(3)


def plotly_clips(bam_path):
    import plotly.offline as py
    import plotly.graph_objs as go

    aln = list(parse_bam(bam_path).items())[0][1]
    aligned_depth = [sum(weight.values()) for weight in aln.weights]
    ins = [sum(i.values()) for i in aln.insertions]
    x_axis = list(range(len(aligned_depth)))
    traces = [
        go.Scattergl(x=x_axis, y=aligned_depth, mode="lines", name="Aligned depth"),
        go.Scattergl(
            x=x_axis, y=aln.clip_depth, mode="lines", name="Soft clip total depth"
        ),
        go.Scattergl(
            x=x_axis, y=aln.clip_start_depth, mode="lines", name="Soft clip start depth"
        ),
        go.Scattergl(
            x=x_axis, y=aln.clip_end_depth, mode="lines", name="Soft clip end depth"
        ),
        go.Scattergl(
            x=x_axis, y=aln.clip_starts, mode="markers", name="Soft clip starts"
        ),
        go.Scattergl(x=x_axis, y=aln.clip_ends, mode="markers", name="Soft clip ends"),
        go.Scattergl(x=x_axis, y=ins, mode="markers", name="Insertions"),
        go.Scattergl(x=x_axis, y=aln.deletions, mode="markers", name="Deletions"),
    ]
    layout = go.Layout(
        xaxis=dict(type="linear", autorange=True),
        yaxis=dict(type="linear", autorange=True),
    )
    fig = go.Figure(data=traces, layout=layout)
    out_fn = os.path.splitext(os.path.split(bam_path)[1])[0]
    py.plot(fig, filename=out_fn + ".plot.html")


if __name__ == "__main__":
    cli.main()

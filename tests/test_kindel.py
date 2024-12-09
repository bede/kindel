import subprocess
from pathlib import Path
import dnaio
from kindel import kindel


bwa_path = Path("tests/data_bwa_mem")
seg_path = Path("tests/data_segemehl")
mm2_path = Path("tests/data_minimap2")
ext_path = Path("tests/data_ext")


bwa_fns = {fn.name: fn for fn in bwa_path.iterdir() if fn.suffix == ".bam"}
seg_fns = {fn.name: fn for fn in seg_path.iterdir() if fn.suffix == ".bam"}
mm2_fns = {fn.name: fn for fn in mm2_path.iterdir() if fn.suffix == ".bam"}
ext_fns = {fn.name: fn for fn in ext_path.iterdir() if fn.suffix == ".bam"}

test_aln = list(kindel.parse_bam(bwa_path / "1.1.sub_test.bam").values())[0]


# UNIT


def test_consensus():
    pos_weight = {"A": 1, "C": 2, "G": 3, "T": 4, "N": 5}
    assert kindel.consensus(pos_weight)[0] == "N"
    assert kindel.consensus(pos_weight)[1] == 5
    assert kindel.consensus(pos_weight)[2] == 0.33
    assert kindel.consensus(pos_weight)[3] is False
    pos_weight_tie = {"A": 5, "C": 5, "G": 3, "T": 4, "N": 1}
    assert kindel.consensus(pos_weight_tie)[2]


def test_merge_by_lcs():
    one = (
        "AACTGCCGCTAGGGGCGCGTTCGGGCTCGCCAACATCTTCAGTCCGGG",
        "GCCGCTAGGGGCGCGTTCGGGCTCGCCAACATCTTCAGTCCGGGCGCTAAGCAGAACA",
    )
    two = (
        "AACTGCCGCTAGGGGCGCGTTCGGGCTCGCCAACATCTTCAGTCCGGGCGCTAAGCAGAACATC",
        "GCAGATACCTACACCACCGGGGGAACTGCCGCTAGGGGCGCGTTCGGGCTCGCCAACATCTTCAGTCCGGGCGCTAAGCAGAACA",
    )
    short = ("AT", "CG")
    assert (
        kindel.merge_by_lcs(*one, min_overlap=7)
        == "AACTGCCGCTAGGGGCGCGTTCGGGCTCGCCAACATCTTCAGTCCGGGCGCTAAGCAGAACA"
    )
    assert (
        kindel.merge_by_lcs(*two, min_overlap=7)
        == "AACTGCCGCTAGGGGCGCGTTCGGGCTCGCCAACATCTTCAGTCCGGGCGCTAAGCAGAACA"
    )
    assert kindel.merge_by_lcs(*short, min_overlap=7) is None


def test_version():
    subprocess.run("kindel version", shell=True, check=True)


# FUNCTIONAL


def test_parse_bam():
    assert test_aln.ref_id == "ENA|EU155341|EU155341.2"
    assert len(test_aln.weights) == 9306


def test_cdrp_consensuses():
    cdrps = kindel.cdrp_consensuses(
        test_aln.weights,
        test_aln.deletions,
        test_aln.clip_start_weights,
        test_aln.clip_end_weights,
        test_aln.clip_start_depth,
        test_aln.clip_end_depth,
        0.1,
        10,
    )
    print(cdrps)
    assert (
        cdrps[0][0].seq
        == "AACTGCCGCTAGGGGCGCGTTCGGGCTCGCCAACATCTTCAGTCCGGGCGCTAAGCAGAACATCCAGCTGATCAACA"
    )
    assert (
        cdrps[0][1].seq
        == "AGCGTCGATGCAGATACCTACACCACCGGGGGAACTGCCGCTAGGGGCGCGTTCGGGCTCGCCAACATCTTCAGTCCGGGCGCTAAGCAGAACA"
    )


def test_consensus_bwa(tmp_path):
    for fn, path in bwa_fns.items():
        print(f"Processing {fn}")
        with dnaio.open(path.with_suffix(".fa"), mode="r") as reader:
            expected_seq = next(iter(reader)).sequence  # Convert reader to an iterator
        subprocess.run(
            f"kindel consensus {path} > {tmp_path / fn}.fa", shell=True, check=True
        )
        with dnaio.open(tmp_path / f"{fn}.fa", mode="r") as reader:
            observed_seq = next(iter(reader)).sequence  # Convert reader to an iterator
        assert observed_seq.upper() == expected_seq.upper()


def test_consensus_bwa_realign(tmp_path):
    for fn, path in bwa_fns.items():
        print(f"Processing {fn}")
        with dnaio.open(path.with_suffix(".realign.fa"), mode="r") as reader:
            expected_seq = next(iter(reader)).sequence  # Convert reader to an iterator
        subprocess.run(
            f"kindel consensus -r {path} > {tmp_path / fn}.realign.fa",
            shell=True,
            check=True,
        )
        with dnaio.open(tmp_path / f"{fn}.realign.fa", mode="r") as reader:
            observed_seq = next(iter(reader)).sequence  # Convert reader to an iterator
        print(observed_seq)
        assert observed_seq.upper() == expected_seq.upper()


def test_consensus_mm2(tmp_path):
    for fn, path in mm2_fns.items():
        print(f"Processing {fn}")
        with dnaio.open(path.with_suffix(".fa"), mode="r") as reader:
            expected_records = {
                record.name: record.sequence for record in iter(reader)
            }  # Convert reader to iterator
        subprocess.run(
            f"kindel consensus {path} > {tmp_path / fn}.fa", shell=True, check=True
        )
        with dnaio.open(tmp_path / f"{fn}.fa", mode="r") as reader:
            observed_records = {
                record.name: record.sequence for record in iter(reader)
            }  # Convert reader to iterator
        for r_id in expected_records:
            assert observed_records[r_id].upper() == expected_records[r_id].upper()


def test_consensus_mm2_realign(tmp_path):
    for fn, path in mm2_fns.items():
        print(f"Processing {fn}")
        with dnaio.open(path.with_suffix(".realign.fa"), mode="r") as reader:
            expected_records = {
                record.name: record.sequence for record in iter(reader)
            }  # Convert reader to iterator
        subprocess.run(
            f"kindel consensus {path} > {tmp_path / fn}.realign.fa",
            shell=True,
            check=True,
        )
        with dnaio.open(tmp_path / f"{fn}.realign.fa", mode="r") as reader:
            observed_records = {
                record.name: record.sequence for record in iter(reader)
            }  # Convert reader to iterator
        for r_id in expected_records:
            assert observed_records[r_id].upper() == expected_records[r_id].upper()


def test_consensus_ext(tmp_path):
    for fn, path in ext_fns.items():
        print(f"Processing {fn}")
        with dnaio.open(path.with_suffix(".fa"), mode="r") as reader:
            expected_records = {
                record.name: record.sequence for record in iter(reader)
            }  # Convert reader to iterator
        subprocess.run(
            f"kindel consensus {path} > {tmp_path / fn}.fa", shell=True, check=True
        )
        with dnaio.open(tmp_path / f"{fn}.fa", mode="r") as reader:
            observed_records = {
                record.name: record.sequence for record in iter(reader)
            }  # Convert reader to iterator
        for r_id in expected_records:
            assert observed_records[r_id].upper() == expected_records[r_id].upper()


def test_consensus_ext_realign(tmp_path):
    for fn, path in ext_fns.items():
        print(f"Processing {fn}")
        with dnaio.open(path.with_suffix(".realign.fa"), mode="r") as reader:
            expected_records = {
                record.name: record.sequence for record in iter(reader)
            }  # Convert reader to iterator
        subprocess.run(
            f"kindel consensus {path} > {tmp_path / fn}.realign.fa",
            shell=True,
            check=True,
        )
        with dnaio.open(tmp_path / f"{fn}.realign.fa", mode="r") as reader:
            observed_records = {
                record.name: record.sequence for record in iter(reader)
            }  # Convert reader to iterator
        for r_id in expected_records:
            assert observed_records[r_id].upper() == expected_records[r_id].upper()


def test_plot():
    subprocess.run(f"kindel plot {bwa_fns['1.1.sub_test.bam']}", shell=True, check=True)
    plot_file = Path("1.1.sub_test.plot.html")
    if plot_file.exists():
        plot_file.unlink()


def test_weights():
    kindel.weights(bwa_fns["1.1.sub_test.bam"])


def test_weights_relative():
    kindel.weights(bwa_fns["1.1.sub_test.bam"], relative=True)


def test_features():
    kindel.features(bwa_fns["1.1.sub_test.bam"])

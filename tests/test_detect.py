import csv
import gzip
import json
from pathlib import Path

import pytest

from eber_scope.detect import correct_barcode, find_anchor_window, reverse_complement, run_detection


REFERENCE = (
    "AGGACCTACGCTGCCCTAGAGGTTTTGCTAGGGAGGAGACGTGTGTGGCTGTAGCCACCCGTCCCGGGTACAAGTCCCGG"
    "GTGGTGAGGACGGTGTCTGTGGTTGTCTTCCCAGACTCTGCTTTCTGCCGTCTTCGGTCAAGTACCAGCTGGTGGTCCGC"
    "ATGTTTT"
)
LEFT = "CTTCCGATCT"
RIGHT = "TTTCTTATATGGG"
BARCODE = "ACGTACGTACGTACGT"


def write_fastq(path, records):
    with Path(path).open("w") as handle:
        for read_id, sequence in records:
            handle.write(f"@{read_id}\n{sequence}\n+\n{'I' * len(sequence)}\n")


@pytest.mark.parametrize("umi_length", [10, 12])
def test_anchor_layouts(umi_length):
    sequence = "GG" + LEFT + BARCODE + "A" * umi_length + RIGHT + "CCC"
    window = find_anchor_window(sequence, LEFT, RIGHT, 16 + umi_length, 1)
    assert sequence[window["insert_start"] : window["insert_end"]] == BARCODE + "A" * umi_length


def test_anchor_mismatch_limit():
    insert = BARCODE + "A" * 12
    one_mismatch = "A" + LEFT[1:] + insert + RIGHT
    two_mismatches = "AA" + LEFT[2:] + insert + RIGHT
    assert find_anchor_window(one_mismatch, LEFT, RIGHT, len(insert), 1) is not None
    assert find_anchor_window(two_mismatches, LEFT, RIGHT, len(insert), 1) is None


def test_barcode_correction_statuses():
    whitelist = {"AAAAAAAAAAAAAAAA", "CAAAAAAAAAAAAAAA"}
    neighbor_map = {"GAAAAAAAAAAAAAAA": whitelist, "AAAAAAAAAAAAAAAT": {"AAAAAAAAAAAAAAAA"}}
    assert correct_barcode("AAAAAAAAAAAAAAAA", whitelist, neighbor_map) == ("AAAAAAAAAAAAAAAA", "exact")
    assert correct_barcode("AAAAAAAAAAAAAAAT", whitelist, neighbor_map) == (
        "AAAAAAAAAAAAAAAA",
        "corrected_one_mismatch",
    )
    assert correct_barcode("GAAAAAAAAAAAAAAA", whitelist, neighbor_map)[1] == "ambiguous_one_mismatch"
    assert correct_barcode("TTTTTTTTTTTTTTTT", whitelist, neighbor_map)[1] == "not_in_whitelist"


@pytest.mark.parametrize(
    ("chemistry", "umi"),
    [("5p-v2", "ACGTACGTAC"), ("5p-v3", "ACGTACGTACGT")],
)
def test_detect_counts_read_pairs_and_unique_umis(tmp_path, chemistry, umi):
    reference_path = tmp_path / "reference.fa"
    reference_path.write_text(f">EBER1\n{REFERENCE}\n")
    whitelist_path = tmp_path / "whitelist.txt"
    whitelist_path.write_text(BARCODE + "-1\n")
    r1_sequence = LEFT + BARCODE + umi + RIGHT + REFERENCE[:70]
    r2_sequence = REFERENCE[70:150]
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    write_fastq(r1_path, [("pair1/1", r1_sequence), ("pair2/1", r1_sequence)])
    write_fastq(r2_path, [("pair1/2", r2_sequence), ("pair2/2", r2_sequence)])

    _, cell_path, qc_path = run_detection(
        r1_path,
        r2_path,
        tmp_path / "result",
        "synthetic_sample",
        reference_path,
        whitelist_path,
        chemistry=chemistry,
    )
    with cell_path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["n_qualifying_read_pairs"] == "2"
    assert rows[0]["n_unique_umis"] == "1"
    assert rows[0]["targeted_positive"] == "1"
    qc = json.loads(qc_path.read_text())
    assert qc["counts"]["qualifying_reference_hits"] == 2


def test_detect_rejects_mismatched_read_names(tmp_path):
    reference_path = tmp_path / "reference.fa"
    reference_path.write_text(f">EBER1\n{REFERENCE}\n")
    insert = BARCODE + "A" * 12
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    write_fastq(r1_path, [("pair1/1", LEFT + insert + RIGHT + REFERENCE[:70])])
    write_fastq(r2_path, [("other/2", REFERENCE[70:150])])
    with pytest.raises(ValueError, match="read-name mismatch"):
        run_detection(r1_path, r2_path, tmp_path / "result", "sample", reference_path)


def test_detect_rejects_mismatched_record_counts(tmp_path):
    reference_path = tmp_path / "reference.fa"
    reference_path.write_text(f">EBER1\n{REFERENCE}\n")
    insert = BARCODE + "A" * 12
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    write_fastq(r1_path, [("pair1/1", LEFT + insert + RIGHT + REFERENCE[:70])])
    write_fastq(r2_path, [("pair1/2", REFERENCE[70:150]), ("pair2/2", REFERENCE[70:150])])
    with pytest.raises(ValueError, match="different numbers"):
        run_detection(r1_path, r2_path, tmp_path / "result", "sample", reference_path)


def test_reverse_complement_r2_is_detected(tmp_path):
    reference_path = tmp_path / "reference.fa"
    reference_path.write_text(f">EBER1\n{REFERENCE}\n")
    insert = BARCODE + "A" * 12
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    write_fastq(r1_path, [("pair1/1", LEFT + insert + RIGHT + REFERENCE[:70])])
    write_fastq(r2_path, [("pair1/2", reverse_complement(REFERENCE[70:150]))])
    _, cell_path, _ = run_detection(r1_path, r2_path, tmp_path / "result", "sample", reference_path)
    with cell_path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1


def test_alignment_threshold_is_inclusive(monkeypatch, tmp_path):
    class FixedThirtyBaseAligner:
        def __init__(self, reference):
            self.reference = reference

        def __call__(self, target):
            return {"target_begin": 0, "target_end_optimal": 30}

    monkeypatch.setattr("eber_scope.detect.StripedSmithWaterman", FixedThirtyBaseAligner)
    reference_path = tmp_path / "reference.fa"
    reference_path.write_text(f">EBER1\n{'ACGT' * 25}\n")
    insert = BARCODE + "A" * 12
    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    write_fastq(r1_path, [("pair1/1", LEFT + insert + RIGHT + "A")])
    write_fastq(r2_path, [("pair1/2", "A" * 30)])

    _, passing_path, _ = run_detection(
        r1_path,
        r2_path,
        tmp_path / "passing",
        "sample",
        reference_path,
        threshold=0.60,
    )
    _, failing_path, _ = run_detection(
        r1_path,
        r2_path,
        tmp_path / "failing",
        "sample",
        reference_path,
        threshold=0.600001,
    )
    with passing_path.open() as handle:
        assert len(list(csv.DictReader(handle, delimiter="\t"))) == 1
    with failing_path.open() as handle:
        assert len(list(csv.DictReader(handle, delimiter="\t"))) == 0

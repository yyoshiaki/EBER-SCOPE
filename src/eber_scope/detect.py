import csv
import gzip
import itertools
import json
from collections import defaultdict
from pathlib import Path

from skbio.alignment import StripedSmithWaterman


CHEMISTRIES = {
    "5p-v2": (16, 10),
    "5p-v3": (16, 12),
}
DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
HIT_COLUMNS = [
    "sample",
    "reference_id",
    "read_id",
    "barcode_raw",
    "barcode",
    "barcode_status",
    "umi",
    "anchor_left_mismatches",
    "anchor_right_mismatches",
    "alignment_span_fraction_r1",
    "alignment_span_fraction_r2",
    "summed_alignment_span_fraction",
    "cdna_length_r1",
    "read_length_r2",
    "qualifying_hit",
]


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def read_fasta(path):
    records = {}
    current_id = None
    sequence = []
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = "".join(sequence).upper()
                current_id = line[1:].split()[0]
                sequence = []
            else:
                if current_id is None:
                    raise ValueError("FASTA sequence encountered before a header")
                sequence.append(line)
    if current_id is not None:
        records[current_id] = "".join(sequence).upper()
    if not records:
        raise ValueError("No FASTA records found")
    return records


def iter_fastq(path):
    with open_text(path) as handle:
        record_number = 0
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().rstrip("\r\n")
            plus = handle.readline()
            quality = handle.readline().rstrip("\r\n")
            record_number += 1
            if not sequence or not plus or not quality:
                raise ValueError(f"Truncated FASTQ record {record_number} in {path}")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ record {record_number} in {path}")
            if len(sequence) != len(quality):
                raise ValueError(f"Sequence/quality length mismatch at record {record_number} in {path}")
            read_id = header[1:].split()[0]
            if read_id.endswith("/1") or read_id.endswith("/2"):
                read_id = read_id[:-2]
            yield read_id, sequence.upper()


def hamming_distance(left, right):
    if len(left) != len(right):
        raise ValueError("Hamming distance requires equal-length sequences")
    return sum(a != b for a, b in zip(left, right))


def find_anchor_window(sequence, left_anchor, right_anchor, insert_length, max_mismatches):
    total_length = len(left_anchor) + insert_length + len(right_anchor)
    for start in range(len(sequence) - total_length + 1):
        left_observed = sequence[start : start + len(left_anchor)]
        left_mismatches = hamming_distance(left_observed, left_anchor)
        if left_mismatches > max_mismatches:
            continue
        right_start = start + len(left_anchor) + insert_length
        right_observed = sequence[right_start : right_start + len(right_anchor)]
        right_mismatches = hamming_distance(right_observed, right_anchor)
        if right_mismatches <= max_mismatches:
            return {
                "stop": right_start + len(right_anchor),
                "insert_start": start + len(left_anchor),
                "insert_end": right_start,
                "left_mismatches": left_mismatches,
                "right_mismatches": right_mismatches,
            }
    return None


def load_whitelist(path, barcode_length):
    if path is None:
        return set(), {}
    whitelist = set()
    with open_text(path) as handle:
        for line in handle:
            barcode = line.strip()
            if barcode.endswith("-1"):
                barcode = barcode[:-2]
            if len(barcode) == barcode_length:
                whitelist.add(barcode.upper())
    if not whitelist:
        raise ValueError("No barcodes of the expected length were found in the whitelist")
    neighbor_map = defaultdict(set)
    for barcode in whitelist:
        for index, base in enumerate(barcode):
            for alternative in "ACGT":
                if alternative != base:
                    neighbor = barcode[:index] + alternative + barcode[index + 1 :]
                    neighbor_map[neighbor].add(barcode)
    return whitelist, neighbor_map


def correct_barcode(barcode, whitelist, neighbor_map):
    if not whitelist:
        return barcode, "not_checked"
    if barcode in whitelist:
        return barcode, "exact"
    matches = neighbor_map.get(barcode, set())
    if len(matches) == 1:
        return next(iter(matches)), "corrected_one_mismatch"
    if len(matches) > 1:
        return barcode, "ambiguous_one_mismatch"
    return barcode, "not_in_whitelist"


def reverse_complement(sequence):
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def alignment_span(alignment):
    return max(0, alignment["target_end_optimal"] - alignment["target_begin"])


def run_detection(
    r1_path,
    r2_path,
    output_prefix,
    sample,
    reference_path,
    barcode_whitelist=None,
    chemistry="5p-v3",
    left_anchor="CTTCCGATCT",
    right_anchor="TTTCTTATATGGG",
    anchor_mismatches=1,
    threshold=0.60,
    max_reads=None,
):
    barcode_length, umi_length = CHEMISTRIES[chemistry]
    references = read_fasta(reference_path)
    aligners = {name: StripedSmithWaterman(sequence) for name, sequence in references.items()}
    whitelist, neighbor_map = load_whitelist(barcode_whitelist, barcode_length)
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    hit_path = Path(f"{output_prefix}.hit_read_pairs.tsv.gz")
    cell_path = Path(f"{output_prefix}.cell_summary.tsv")
    qc_path = Path(f"{output_prefix}.qc.json")

    counts = {
        "read_pairs_examined": 0,
        "read_pairs_with_anchors": 0,
        "qualifying_reference_hits": 0,
        "qualifying_hits_with_usable_barcode": 0,
    }
    barcode_status_counts = defaultdict(int)
    cell_hits = defaultdict(lambda: {"read_pairs": 0, "umis": set(), "max_fraction": 0.0})
    r1_records = iter_fastq(r1_path)
    r2_records = iter_fastq(r2_path)
    if max_reads is not None:
        r1_records = itertools.islice(r1_records, max_reads)
        r2_records = itertools.islice(r2_records, max_reads)

    with gzip.open(hit_path, "wt", newline="") as hit_handle:
        writer = csv.DictWriter(hit_handle, fieldnames=HIT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for pair_number, pair in enumerate(itertools.zip_longest(r1_records, r2_records), start=1):
            r1_record, r2_record = pair
            if r1_record is None or r2_record is None:
                raise ValueError("R1 and R2 contain different numbers of records")
            r1_id, r1_sequence = r1_record
            r2_id, r2_sequence = r2_record
            if r1_id != r2_id:
                raise ValueError(f"R1/R2 read-name mismatch at pair {pair_number}: {r1_id} != {r2_id}")
            counts["read_pairs_examined"] += 1
            window = find_anchor_window(
                r1_sequence,
                left_anchor,
                right_anchor,
                barcode_length + umi_length,
                anchor_mismatches,
            )
            if window is None:
                continue
            counts["read_pairs_with_anchors"] += 1
            insert = r1_sequence[window["insert_start"] : window["insert_end"]]
            barcode_raw = insert[:barcode_length]
            umi = insert[barcode_length:]
            barcode, barcode_status = correct_barcode(barcode_raw, whitelist, neighbor_map)
            barcode_status_counts[barcode_status] += 1
            r2_reverse = reverse_complement(r2_sequence)

            for reference_id, aligner in aligners.items():
                r1_span = alignment_span(aligner(r1_sequence))
                r2_span = max(alignment_span(aligner(r2_sequence)), alignment_span(aligner(r2_reverse)))
                reference_length = len(references[reference_id])
                summed_fraction = min(reference_length, r1_span + r2_span) / reference_length
                if summed_fraction < threshold:
                    continue
                counts["qualifying_reference_hits"] += 1
                usable_barcode = barcode_status in {"not_checked", "exact", "corrected_one_mismatch"}
                if usable_barcode:
                    counts["qualifying_hits_with_usable_barcode"] += 1
                    cell = cell_hits[(reference_id, barcode)]
                    cell["read_pairs"] += 1
                    cell["umis"].add(umi)
                    cell["max_fraction"] = max(cell["max_fraction"], summed_fraction)
                writer.writerow(
                    {
                        "sample": sample,
                        "reference_id": reference_id,
                        "read_id": r1_id,
                        "barcode_raw": barcode_raw,
                        "barcode": barcode,
                        "barcode_status": barcode_status,
                        "umi": umi,
                        "anchor_left_mismatches": window["left_mismatches"],
                        "anchor_right_mismatches": window["right_mismatches"],
                        "alignment_span_fraction_r1": f"{r1_span / reference_length:.6f}",
                        "alignment_span_fraction_r2": f"{r2_span / reference_length:.6f}",
                        "summed_alignment_span_fraction": f"{summed_fraction:.6f}",
                        "cdna_length_r1": len(r1_sequence) - window["stop"],
                        "read_length_r2": len(r2_sequence),
                        "qualifying_hit": 1,
                    }
                )

    cell_columns = [
        "sample",
        "reference_id",
        "barcode",
        "n_qualifying_read_pairs",
        "n_unique_umis",
        "max_summed_alignment_span_fraction",
        "targeted_positive",
    ]
    with cell_path.open("w", newline="") as cell_handle:
        writer = csv.DictWriter(cell_handle, fieldnames=cell_columns, delimiter="\t")
        writer.writeheader()
        for (reference_id, barcode), values in sorted(cell_hits.items()):
            writer.writerow(
                {
                    "sample": sample,
                    "reference_id": reference_id,
                    "barcode": barcode,
                    "n_qualifying_read_pairs": values["read_pairs"],
                    "n_unique_umis": len(values["umis"]),
                    "max_summed_alignment_span_fraction": f"{values['max_fraction']:.6f}",
                    "targeted_positive": 1,
                }
            )

    qc = {
        "tool": "EBER-SCOPE",
        "version": "0.1.0",
        "sample": sample,
        "inputs": {
            "r1": Path(r1_path).name,
            "r2": Path(r2_path).name,
            "reference": Path(reference_path).name,
            "barcode_whitelist": Path(barcode_whitelist).name if barcode_whitelist else None,
        },
        "parameters": {
            "chemistry": chemistry,
            "barcode_length": barcode_length,
            "umi_length": umi_length,
            "left_anchor": left_anchor,
            "right_anchor": right_anchor,
            "anchor_mismatches": anchor_mismatches,
            "summed_alignment_span_fraction_threshold": threshold,
            "cell_positive_minimum_qualifying_read_pairs": 1,
            "max_reads": max_reads,
        },
        "counts": counts,
        "barcode_status_counts": dict(sorted(barcode_status_counts.items())),
        "targeted_positive_cells": len(cell_hits),
    }
    qc_path.write_text(json.dumps(qc, indent=2) + "\n")
    return hit_path, cell_path, qc_path

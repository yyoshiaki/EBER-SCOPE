import pandas as pd

from eber_scope.merge_calls import merge_calls


def test_merge_calls_builds_union_and_sources(tmp_path):
    existing = pd.DataFrame(
        {
            "sample": ["sample1"] * 4,
            "barcode": ["AAAA-1", "CCCC-1", "GGGG-1", "TTTT-1"],
            "existing_positive": [True, False, True, False],
        }
    )
    targeted = pd.DataFrame(
        {
            "sample": ["sample1", "sample1"],
            "reference_id": ["EBER1", "EBER1"],
            "barcode": ["CCCC", "GGGG"],
            "n_qualifying_read_pairs": [1, 2],
            "n_unique_umis": [1, 2],
            "max_summed_alignment_span_fraction": [0.7, 0.8],
            "targeted_positive": [1, 1],
        }
    )
    existing_path = tmp_path / "existing.tsv"
    targeted_path = tmp_path / "targeted.tsv"
    output_path = tmp_path / "combined.tsv"
    existing.to_csv(existing_path, sep="\t", index=False)
    targeted.to_csv(targeted_path, sep="\t", index=False)
    combined = merge_calls(existing_path, targeted_path, output_path)
    assert combined["detection_source"].tolist() == ["existing_only", "targeted_only", "both", "negative"]
    assert combined["union_positive"].tolist() == [True, True, True, False]

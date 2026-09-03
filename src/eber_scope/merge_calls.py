from pathlib import Path

import pandas as pd


def parse_boolean(series):
    normalized = series.astype(str).str.strip().str.lower()
    valid = {"true", "false", "1", "0", "yes", "no", "y", "n"}
    unexpected = sorted(set(normalized) - valid)
    if unexpected:
        raise ValueError(f"Unrecognized boolean values: {unexpected}")
    return normalized.isin({"true", "1", "yes", "y"})


def normalize_barcode(series):
    return series.astype(str).str.replace(r"-1$", "", regex=True)


def merge_calls(existing_path, targeted_path, output_path):
    existing = pd.read_csv(existing_path, sep="\t", dtype={"sample": str, "barcode": str})
    targeted = pd.read_csv(targeted_path, sep="\t", dtype={"sample": str, "barcode": str})
    required_existing = {"sample", "barcode", "existing_positive"}
    required_targeted = {
        "sample",
        "barcode",
        "n_qualifying_read_pairs",
        "n_unique_umis",
        "max_summed_alignment_span_fraction",
        "targeted_positive",
    }
    if missing := sorted(required_existing - set(existing.columns)):
        raise ValueError(f"Existing-call TSV is missing columns: {missing}")
    if missing := sorted(required_targeted - set(targeted.columns)):
        raise ValueError(f"Targeted-call TSV is missing columns: {missing}")

    existing = existing.copy()
    targeted = targeted.copy()
    existing["barcode"] = normalize_barcode(existing["barcode"])
    targeted["barcode"] = normalize_barcode(targeted["barcode"])
    existing["existing_positive"] = parse_boolean(existing["existing_positive"])
    targeted["targeted_positive"] = parse_boolean(targeted["targeted_positive"])
    if existing.duplicated(["sample", "barcode"]).any():
        raise ValueError("Existing-call TSV has duplicate sample/barcode rows")
    if targeted.duplicated(["sample", "barcode"]).any():
        raise ValueError("Targeted-call TSV has duplicate sample/barcode rows")

    targeted_columns = [
        "sample",
        "barcode",
        "n_qualifying_read_pairs",
        "n_unique_umis",
        "max_summed_alignment_span_fraction",
        "targeted_positive",
    ]
    combined = existing.merge(targeted[targeted_columns], on=["sample", "barcode"], how="left")
    combined["n_qualifying_read_pairs"] = combined["n_qualifying_read_pairs"].fillna(0).astype(int)
    combined["n_unique_umis"] = combined["n_unique_umis"].fillna(0).astype(int)
    combined["max_summed_alignment_span_fraction"] = combined[
        "max_summed_alignment_span_fraction"
    ].fillna(0.0)
    combined["targeted_positive"] = combined["targeted_positive"].fillna(False).astype(bool)
    combined["union_positive"] = combined["existing_positive"] | combined["targeted_positive"]
    combined["detection_source"] = "negative"
    combined.loc[combined["existing_positive"] & ~combined["targeted_positive"], "detection_source"] = "existing_only"
    combined.loc[~combined["existing_positive"] & combined["targeted_positive"], "detection_source"] = "targeted_only"
    combined.loc[combined["existing_positive"] & combined["targeted_positive"], "detection_source"] = "both"
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_path, sep="\t", index=False)
    return combined

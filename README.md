# EBER-SCOPE

EBER-SCOPE detects EBER1-positive cell barcodes in targeted 10x 5′ paired-end
FASTQs and can combine those calls with an existing single-cell EBV call set.
Version 0.1 is intentionally EBER1-specific and reproduces the decision rule
used in the associated study: a cell is targeted-positive when it has at least
one qualifying read pair.

## Install

```bash
conda env create -f environment.yml
conda activate eber-scope
pip install -e .
```

## Detect EBER1

```bash
eber-scope detect \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --sample sample_id \
  --barcode-whitelist barcodes.txt.gz \
  --chemistry 5p-v3 \
  --output-prefix results/sample
```

The bundled EBER1 reference and the study thresholds are defaults. Use
`--chemistry 5p-v2` for the 16-bp cell barcode plus 10-bp UMI layout; `5p-v3`
uses a 16-bp barcode plus 12-bp UMI. The v3/GEM-X configuration was used in
the production analysis. The v2 preset is covered by synthetic tests but was
not empirically validated in that analysis.

Outputs:

- `*.hit_read_pairs.tsv.gz`: qualifying read pairs, for audit/QC.
- `*.cell_summary.tsv`: targeted-positive cell barcodes with both read-pair
  counts and unique-UMI counts.
- `*.qc.json`: input and barcode-processing summary.

`n_qualifying_read_pairs` is the paper-defining quantity. Unique UMIs are
reported as QC only; they are not used to change the positivity call.

## Combine with existing calls

Prepare a TSV with `sample`, `barcode`, and `existing_positive`, then run:

```bash
eber-scope merge \
  --existing existing_calls.tsv \
  --targeted results/sample.cell_summary.tsv \
  --output combined_calls.tsv
```

The output reports `existing_only`, `targeted_only`, `both`, or `negative`,
plus the union call. Targeted calls are left-joined to the existing cell table;
barcodes absent from that table are intentionally not added.

## Method notes

- The default anchors allow one mismatch per anchor.
- EBER1 support is the sum of local-alignment spans from full R1 and the best
  R2 orientation, capped at the reference length; the default threshold is
  0.60. This matches the production calculation and is not unique-base
  coverage or an independent identity threshold.
- Whitelist barcodes are accepted exactly or corrected only when a one-base
  neighbor maps uniquely. Ambiguous/unmatched barcodes remain in read-level QC
  but are excluded from cell calls.
- R1/R2 record counts and read names are checked before analysis.

Do not commit human FASTQs, cell-level outputs, sample manifests, or protected
metadata to this repository. Tests use synthetic reads only.

## Validation

The packaged workflow passed 11 synthetic regression tests and an equivalence
check on 2,000,000 read pairs from a representative production 10x 5′ v3
library. The legacy and packaged workflows identified the same 17 qualifying
read pairs and 15 usable positive cell barcodes, with exact agreement across
read-level fields and cell-level summaries. The privacy-safe aggregate report
is provided in `validation/260903_production_subset_equivalence.json`; no FASTQ,
cell barcode, UMI, sample identifier, or internal path is included.

## Citation

When using EBER-SCOPE for this workflow, please cite:

> Yasumizu Y, Kim N, Rivier CA, et al. A Genetically Driven Immunologic
> Mechanism Underlying the Link between EBV and Multiple Sclerosis. medRxiv.
> 2026. https://doi.org/10.64898/2025.12.11.25342083

Software citation metadata are provided in `CITATION.cff`.

## License

EBER-SCOPE is released under the MIT License.

# Reproducibility check

Status: **PASS_LOCAL_RELEASE_CANDIDATE**

## Frozen-input integrity

- Frozen manifest: `manuscript_freeze/RESULTS_FREEZE_MANIFEST.tsv` (565 rows).
- Local frozen files: 559/559 present and SHA256-matched; 0 missing; 0 mismatched.
- Server-key frozen files: 6/6 were rechecked against read-only verified local cache copies; 0 missing; 0 mismatched.
- No frozen analytical result, model, ROI, module, gene set, parameter, figure, or evidence grade was recalculated or modified.

## Package inventory checks

- 816 accession/file records include official URL, download date or provenance note, byte size/checksum where available, and raw-data exclusion status.
- 1,208 patient/donor-sample-library mapping rows retain replicate and pair structure.
- 13 dataset-independence rows retain one-cohort counting rules and overlap limitations.
- The canonical M02 set contains exactly 36 unique genes; the 35-gene INPP5D sensitivity and frozen cohort-common mapped sets remain separately labeled.
- 159 source scripts were copied from the frozen source set; `.gitkeep` and compiled `.pyc` caches were excluded.
- 39 small source-data/audit files are included; no raw expression matrix or large analysis object is included.
- Stage 10F and 10G skip/closure records are present.
- Nine frozen negative, skipped, or NOT_ESTIMABLE rows are explicitly indexed.

## Clean-environment minimum reproduction

A fresh copy of `release_candidate/` was made under the system temporary directory and checked with:

```text
python reproduce_minimal.py
```

Result: `PASS` (exit code 0). The verifier uses only the Python standard library and checks every entry in `SHA256SUMS`, required files, canonical module cardinality, dataset independence, script mapping, branch-closure records, and the explicit boundary-result index. This is a provenance reproduction only; it does not rerun bioinformatics analyses or inspect raw data.

## Reproduction limits

Public raw/processed source data must be reacquired from `DATA_ACCESSIONS.tsv`. Some historical bulk cohorts remain sample-scoped because patient identities cannot be recovered. Figshare primary 36/36 M02 and patient-level spatial inference remain NOT_ESTIMABLE. These limitations are preserved rather than repaired or recoded.

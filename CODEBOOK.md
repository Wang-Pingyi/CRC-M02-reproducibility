# Codebook

## Core tables

- `DATA_ACCESSIONS.tsv`: official accession/DOI, asset, URL, download date, byte size, checksum, license note, and verification status. Raw files are not included.
- `PATIENT_SAMPLE_MAP.tsv`: public/de-identified patient or donor, sample, library, biological sample, replicate and pair mapping. Technical replicates and nested spatial units are explicitly identified.
- `DATASET_INDEPENDENCE.tsv`: one-row-per-dataset reuse, overlap risk, permitted counting rule and Stage 10I audit status.
- `GENE_ID_MAPPING.tsv`: canonical M02 identities and cohort-specific feature mapping. Missing genes remain missing; no proxy rescue is used.
- `COHORT_COMMON_GENESETS.tsv`: canonical 36-gene module, frozen 35-gene INPP5D sensitivity, and cohort-common mapped sets with weights and provenance.
- `SCRIPT_OUTPUT_MAP.tsv`: released source code, original and sanitized hashes, role, direct inputs where detectable, and stage-scoped outputs/figures.

## Semantics

- `patient/donor` is the biological inference unit. Sample, library, cell, array, ROI and spot are not silently promoted to independent patients.
- `sample-scoped` means patient identity or independence could not be verified.
- `NOT_ESTIMABLE` means the prespecified estimand or coverage gate could not be evaluated; it is neither positive nor negative.
- `SKIPPED_NOT_ESTIMABLE` and `SKIPPED_GATE_NOT_MET` are governance outcomes, not biological null results.
- `NA` means unavailable, inapplicable, or not recoverable from public records; it is never encoded as zero.
- `primary` and `sensitivity` roles are frozen and must not be exchanged based on result direction.

## Source-data layout

- `source_data/manuscript_figures/`: Stage 11B per-figure TSVs.
- `source_data/tables/main/` and `source_data/tables/supplementary/`: frozen table inputs and audit material.
- `source_data/stage10fg/`: formal Stage 10F/10G skip and branch-closure records.
- `source_data/NEGATIVE_AND_NOT_ESTIMABLE.tsv`: index of explicit negative, skipped and non-estimable frozen rows.

## Security and portability

All text copied into this package is scanned for credentials, private-key material, `.env` files, private host/user identifiers and absolute private paths. Portable placeholders are documented in `README.md`. Checksums cover every package file except `SHA256SUMS` itself.

# CRC M02 cross-dataset reproducibility package

This repository is the public reproducibility package for **“Exploratory cross-dataset evaluation of a locked FAP/adenoma-associated epithelial stem/progenitor transcriptional module.”** It was assembled from the Stage 11C release candidate on 2026-08-01. It contains no raw expression matrices, FASTQ, CEL, H5/H5AD, RDS, caches, credentials, SSH keys, or directly identifying participant information.

## Scientific scope

The project is an exploratory cross-dataset transcriptomic evaluation of epithelial-state changes across normal mucosa, adenoma/polyp and colorectal cancer. The prespecified Stage 6A patient-level discovery screen was negative. M02 is therefore retained only as the frozen **FAP/adenoma-associated exploratory epithelial stem/progenitor module** under the Level C claim ceiling. Spatial patient-level inference was not estimable; Stage 10F and 10G were skipped.

## Recommended reading order

1. `README.md` and `CODEBOOK.md`.
2. `DATA_ACCESSIONS.tsv`, `PATIENT_SAMPLE_MAP.tsv`, and `DATASET_INDEPENDENCE.tsv`.
3. `GENE_ID_MAPPING.tsv` and `COHORT_COMMON_GENESETS.tsv`.
4. `source_data/tables/main/Table2_Frozen_results.tsv` and `source_data/NEGATIVE_AND_NOT_ESTIMABLE.tsv`.
5. `SCRIPT_OUTPUT_MAP.tsv`, `environment/`, and `REPRODUCIBILITY_CHECK.md`.
6. `SECURITY_SCAN.md`, then verify `SHA256SUMS`.

## Reproduction

Run `python reproduce_minimal.py` from this directory. It uses only the Python standard library and verifies checksums, required files, the 36-gene canonical module, dataset independence, script coverage, Stage 10F/10G closure, and explicit negative/NOT_ESTIMABLE records. This is a minimum provenance reproduction, not a rerun of the bioinformatics analyses.

## Path placeholders

Release copies replace private absolute paths with `${CRC_PROJECT_ROOT}`, `${CRC_LOCAL_PROJECT_ROOT}`, `${SYSTEM_TEMP}`, `${SSH_HOST}`, or `${SSH_USER}`. Original executed-script SHA256 values remain in `SCRIPT_OUTPUT_MAP.tsv`; sanitized release-copy hashes are also reported.

## Integrity and exclusions

The frozen Stage 11A manifest contains 565 entries across 26 stage labels. Raw data are accessed through the official URLs and checksums in `DATA_ACCESSIONS.tsv` and are deliberately absent. `NA` is never converted to zero; `NOT_ESTIMABLE` is not interpreted as a negative result.

## Licensing and citation

- Software and analysis scripts are released under the MIT License (`LICENSE`).
- Original project documentation and small source-data tables are released under CC BY 4.0 (`LICENSE-DATA.md`).
- Third-party datasets, repository metadata and reproduced source material remain under their original terms; this repository does not relicense or redistribute the raw public datasets.
- Citation metadata are provided in `CITATION.cff`. Use the Zenodo DOI in the repository release metadata once the archive is registered.

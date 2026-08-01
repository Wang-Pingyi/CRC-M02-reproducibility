# Security scan

Status: **PASS**

Scope: all 225 files in the local release candidate (approximately 12.2 MB).

## Exclusions and findings

- Credential assignments (`password`, `api_key`): 0 matches.
- GitHub/PAT-like tokens: 0 matches.
- PEM/OpenSSH/RSA/EC private-key blocks: 0 matches.
- Private SSH host and username identifiers: 0 matches.
- Windows user-profile absolute paths: 0 matches.
- Private server project/storage absolute paths: 0 matches.
- `.env`, private-key, FASTQ/FQ, CEL, H5/H5AD, RDS, MTX, BAM or CRAM files: 0 files.
- Files larger than 20 MB: 0 files.

## Portability treatment

Release copies of historical scripts and small source tables use documented placeholders for private paths. `SCRIPT_OUTPUT_MAP.tsv` preserves each original frozen SHA256 and separately reports the sanitized release-copy SHA256. The original repository files were not edited by this portability treatment.

## Identity scope

Only public, de-identified dataset patient/donor/sample labels needed to reproduce grouping and pairing are included. No direct patient identity information is present. Official repository URLs and public accession/DOI identifiers are retained deliberately.

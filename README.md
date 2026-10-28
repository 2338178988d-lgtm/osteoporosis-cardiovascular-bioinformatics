> Public code release draft — NOT YET PUBLISHED

# Osteoporosis–cardiovascular comorbidity bioinformatics analysis

This directory is a reviewable draft of the code and reproducibility materials
for the planned public repository:
`2338178988d-lgtm/osteoporosis-cardiovascular-bioinformatics`.
It is not a final repository and has not been uploaded to GitHub, Zenodo, or
any other public archive. Project-authored code is intended to use the
BSD-3-Clause license; third-party materials retain their own terms.

## Included

- `analysis/` — analysis scripts, module READMEs, download/verification scripts,
  and figure-generation code.
- `figures/` — selected rendered figures included for provenance review.
- `manifests/` — endpoint and data-source metadata.

## Deliberately excluded

- Raw and processed participant-level or expression data.
- GWAS summary-statistics files, FinnGen files, NHANES files, GEO matrices,
  reference panels, and PLINK binaries.
- `.venv/`, `.cache/`, `.wisp/`, `submission/`, `uploads/`, `runs/`, archives,
  and internal project state.
- Author identities, ethics determinations, submission declarations, and
  private access credentials.

## Current reproducibility status

The frozen inventory records 96 planned inputs, 91 present and size-validated,
and 5 missing or unavailable. The complete pipeline has not been independently
rerun in the current project. Raw inputs should be obtained from their official
accession/URL sources using the project download and verification scripts.

## Before publication

1. Confirm third-party data, software, font, and figure-runtime terms.
2. Confirm the final BSD-3-Clause scope for project-authored code.
3. Add a pinned release tag/commit and archive the release with Zenodo if a DOI
   is required.
4. Run the download, manuscript, and clean-environment verification checks.
5. Obtain author confirmation for the final Data Availability and Code
   Availability statements.

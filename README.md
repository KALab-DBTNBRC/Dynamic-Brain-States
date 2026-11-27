# Dynamic Brain States

Companion code for **Familial and Genetic Contributions to Dynamic Brain-State Organization: Association with Alcohol Consumption in a Twin Cohort**.

## Repository guide

| Folder | Manuscript analysis |
|---|---|
| `01_component_selection` | ICA quality features, classification, retained component selection (steps 01–04) |
| `02_clustering` | Concatenation and consensus stability (steps 05–08; Supplementary Figure S2 / Table S1) |
| `03_spatial_characterization` | Centroid visualization, reconstruction and atlas overlap (steps 10–14; Figure 1, Supplementary Figures S1/S3) |
| `04_statistical_analysis` | ACE models and continuous/ordinal regression (steps 15 and 17) |
| `05_figures` | Main Figures 1–3 and Supplementary Figures S1/S3 |
| `06_tables` | Main Tables 1–2 and Supplementary Tables S2–S5 |
| `supplementary_data` | Supplied Supplementary Files S1–S3 |
| `legacy` | Earlier clustering, Pearson heritability and alternate ACE table scripts; not the final methods |

Run from the repository root unless a script specifies otherwise. Read each script's usage header. The step numbers are identifiers, not a fully automated execution sequence: centroid reconstruction (12) may precede resampling (11). Python 3.12.10, R 4.4.1 and MATLAB R2023b are the manuscript-reported versions. Python packages include numpy, scipy, pandas, nibabel, nilearn, matplotlib, numba, statsmodels, openpyxl, python-calamine, Pillow and python-docx. R scripts declare their packages; MATLAB requires Statistics and Machine Learning Toolbox and, for the legacy clustering script, TbCAPs.

## Reproducibility status

This is an organized source snapshot, not a verified end-to-end pipeline. Read REVIEW.md before public release. Several scripts require local paths and missing inputs; no model was rerun. The journal submission and cover letter are deliberately outside the repository. Confirm HCP data-use permissions before publishing supplementary CSVs. No license has been assigned; authors should select one.

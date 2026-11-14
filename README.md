# Dynamic-Brain-States

Analysis code for a twin study of dynamic brain states. Starting from group-ICA time courses of resting-state fMRI in monozygotic (MZ) and dizygotic (DZ) twins, the pipeline:

1. removes noise components,
2. identifies recurring brain states by consensus clustering (K chosen from consensus metrics, K = 4 in the later steps),
3. maps each state onto functional networks and anatomical regions, and
4. estimates the heritability of state dynamics (fractional occupancy, dwell time, appearance rate) and their association with alcohol use and odor identification.

The repository contains code only. No imaging or phenotype data are included.

## Repository layout

```
Dynamic-Brain-States/
├── src/
│   ├── 01_ic_noise_classification/      Steps 1-3   (Python)
│   ├── 02_preprocessing/                Steps 4-5   (Python, MATLAB)
│   ├── 03_consensus_clustering/         Steps 6-9   (MATLAB, Python)
│   ├── 04_brain_state_characterization/ Steps 10-14 (Python)
│   └── 05_twin_analysis/                Steps 15-17 (R, Python)
├── requirements.txt                     Python dependencies
├── LICENSE
└── README.md
```

No script contains a hard-coded path. Every input and output location is passed as an argument, so the code runs from any directory and on any machine.

## Pipeline

| Step | Script | Language | What it does | Main inputs | Main outputs |
|------|--------|----------|--------------|-------------|--------------|
| 1 | `01_ic_noise_classification/step01_calculate_fractions.py` | Python | Fraction of each of 300 ICs lying in CSF and at the brain edge | ICA maps, CSF mask, edge mask | `edge_fractions.csv`, `csf_fractions.csv` |
| 2 | `01_ic_noise_classification/step02_motion_highfreq_features.py` | Python | Max correlation with motion regressors and high-frequency power per IC and subject | IC time courses (`.txt`), motion regressors (`.xlsx`) | `group_component_features.csv`, `aggregated_component_features.csv` |
| 3 | `01_ic_noise_classification/step03_classify_components.py` | Python | Labels each IC Signal or Noise using fixed thresholds | Outputs of steps 1 and 2 | `classification_results.csv` |
| 4 | `02_preprocessing/step04_select_rois.py` | Python | Reduces each subject's 300-column time series to the 64 retained ROIs | Per-subject `.mat` files | Cleaned `.mat` files |
| 5 | `02_preprocessing/step05_concatenate_subjects.m` | MATLAB | Stacks all subjects into one matrix | Cleaned `.mat` files | `all_subjects_twins.mat` |
| 6 | `03_consensus_clustering/step06_consensus_metrics.m` | MATLAB | Exact consensus-clustering metrics (PAC, CDF area, ΔK) for K = 2-20 | `X_prototypes` (240 subjects × 200 prototypes × 64 ROIs) | Metrics `.mat`/`.csv`, quick-look figures |
| 7 | `03_consensus_clustering/consensus_core.py` | Python | Library used by step 8 (correlation k-means, consensus metrics) | | |
| 8 | `03_consensus_clustering/step08_consensus_metrics_python.py` | Python | Python re-run of step 6 (main, repeat and null runs) | `P` array (240 × 200 × 64) | `consensus_metrics_<real\|repeat\|null>_python.mat` |
| 9 | `03_consensus_clustering/step09_run_clustering.m` | MATLAB | Consensus clustering of the concatenated 64-ROI data with TbCAPs | `all_subjects` matrix | `Consensus` struct (`.mat`) |
| 10 | `04_brain_state_characterization/step10_plot_centroids.py` | Python | Bar plot of each state centroid | Centroids `C` (`.mat`) | `centroid_state_<k>.png` |
| 11 | `04_brain_state_characterization/step11_resample_to_yeo.py` | Python | Resamples the 4D centroid image to the Yeo atlas grid | Step 12 output, Yeo atlas | Resampled NIfTI |
| 12 | `04_brain_state_characterization/step12_centroids_to_4d_nifti.py` | Python | Projects centroids onto the 64 IC maps to get one brain map per state | Centroids (v7.3 `.mat`), IC maps | `centroids_4D.nii.gz` |
| 13 | `04_brain_state_characterization/step13_yeo_network_overlap.py` | Python | Spatial correlation of each state with each Yeo 17 network | Yeo atlas, step 11 output | `correlation_matrix_with_p_values.csv` |
| 14 | `04_brain_state_characterization/step14_atlas_region_overlap.py` | Python | Contribution of Harvard-Oxford regions to each state | Step 12 output, Harvard-Oxford atlases | Two `.xlsx` workbooks |
| 15 | `05_twin_analysis/step15_ace_models.R` | R | ACE / AE / CE / E twin models (OpenMx) for 23 phenotypes | Twin workbook (`.xlsx`) | `ACE_Master_Results.{rds,csv}`, `ACE_Model_Diagnostics.csv` |
| 16 | `05_twin_analysis/step16_heritability_statewise.py` | Python | State-wise MZ vs DZ similarity and Falconer heritability | `fo.csv`, `dt.csv`, `ar.csv` | Results CSV |
| 17 | `05_twin_analysis/step17_glm_models.R` | R | Brain dynamics vs alcohol use and odor identification (mixed models, FDR) | MZ and DZ workbooks | Results `.xlsx` |

Steps run in numeric order. Step 12 must run before steps 11, 13 and 14, because it creates the 4D centroid image they use.

Some inputs are produced outside this repository: the group-ICA maps and time courses (FSL MELODIC), the motion regressors, the subject-level prototypes (`X_prototypes`/`P`) and the final centroid matrix `C` (`BrainStates_MedoidPartition.mat`). The script that creates the publication figures for step 6 is also not included.

## Setup

### Python (3.10 or newer)

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### MATLAB (R2019b or newer)

Steps 5, 6 and 9 are functions that use `arguments` blocks. Step 6 needs the Statistics and Machine Learning Toolbox. Step 9 needs the TbCAPs toolbox (`CAP_ConsensusClustering_efficient`) on the MATLAB path.

### R (4.x)

Steps 15 and 17 install their packages on first run: `OpenMx`, `readxl`, `knitr` (step 15) and `readxl`, `dplyr`, `tibble`, `lme4`, `lmerTest`, `sandwich`, `lmtest`, `openxlsx` (step 17).

## Usage

Paths below are examples. Use any locations you like. Run the commands from the repository root.

```bash
# Steps 1-3: noise classification
python src/01_ic_noise_classification/step01_calculate_fractions.py \
    --ica-maps data/melodic_IC300_sum.nii.gz \
    --csf-mask data/brain_extracted_csf_mask.nii.gz \
    --edge-mask data/outer_edge_mask.nii.gz \
    --output-dir results/01_noise

python src/01_ic_noise_classification/step02_motion_highfreq_features.py \
    --time-series-dir data/ic_timeseries --regressors-dir data/regressors \
    --output-dir results/01_noise          # add --fs <1/TR> if your TR differs; see Known issues

python src/01_ic_noise_classification/step03_classify_components.py \
    --input-dir results/01_noise --output-dir results/01_noise

# Step 4: keep 64 of 300 ROIs
python src/02_preprocessing/step04_select_rois.py \
    --input-dir data/twin_300 --output-dir data/twin_64
```

```matlab
% Steps 5, 6, 9 (MATLAB)
addpath('src/02_preprocessing', 'src/03_consensus_clustering');
step05_concatenate_subjects({'data/twin_64/dizygotic', 'data/twin_64/monozygotic'}, 'data/all_subjects_twins64.mat');
step06_consensus_metrics('data/X_prototypes.mat', 'results/03_consensus');          % real data
step06_consensus_metrics('data/X_prototypes.mat', 'results/03_consensus', true);    % ROI-permuted null
step09_run_clustering('data/all_subjects_twins64.mat', 'results/03_consensus/ConsensusResults.mat');
```

```bash
# Step 8: Python version of the consensus metrics (~30 min per run on 2 cores)
python src/03_consensus_clustering/step08_consensus_metrics_python.py data/P.mat results/03_consensus
python src/03_consensus_clustering/step08_consensus_metrics_python.py data/P.mat results/03_consensus --repeat
python src/03_consensus_clustering/step08_consensus_metrics_python.py data/P.mat results/03_consensus --null

# Steps 10, 12, 11, 13, 14: characterize the states
python src/04_brain_state_characterization/step10_plot_centroids.py \
    --centroids-mat data/centroids.mat --output-dir results/04_states/plots
python src/04_brain_state_characterization/step12_centroids_to_4d_nifti.py \
    --centroids-mat data/BrainStates_MedoidPartition.mat --ic-maps data/melodic_IC_64_final.nii.gz \
    --output results/04_states/centroids_4D.nii.gz
python src/04_brain_state_characterization/step11_resample_to_yeo.py \
    --centroids results/04_states/centroids_4D.nii.gz --yeo-atlas data/atlas/Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz \
    --output results/04_states/centroids_resampled_to_Yeo.nii.gz
python src/04_brain_state_characterization/step13_yeo_network_overlap.py \
    --yeo-atlas data/atlas/Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz \
    --centroids results/04_states/centroids_resampled_to_Yeo.nii.gz --output-dir results/04_states/yeo
python src/04_brain_state_characterization/step14_atlas_region_overlap.py \
    --centroids results/04_states/centroids_4D.nii.gz \
    --cortical-atlas data/atlas/HarvardOxford-cort-maxprob-thr25-2mm.nii.gz \
    --subcortical-atlas data/atlas/HarvardOxford-sub-maxprob-thr25-2mm.nii.gz \
    --output-dir results/04_states/atlas

# Steps 15-17: twin analyses
Rscript src/05_twin_analysis/step15_ace_models.R data/Twin_Output.xlsx results/05_twin/ace
python src/05_twin_analysis/step16_heritability_statewise.py \
    --data-dir data/heritability --output results/05_twin/Twin_brain_dynamics_statewise_results.csv
Rscript src/05_twin_analysis/step17_glm_models.R data/Mono_GLM_Sheet-all-states.xlsx data/DI_GLM_sheet-all-states.xlsx results/05_twin/glm
```

Run `python <script> --help` for the full argument list of any Python script.

Keep input data in `data/` and write outputs to `results/`. Both folders are listed in `.gitignore`, so imaging data and phenotypes are not committed by accident.

## Expected input formats

- **Step 2:** `<subject>.txt` holds 4800 × 300 values (time points × ICs). `<subject>_*.xlsx` holds the motion regressors, with the first 12 columns used. Subject IDs are the text before the first `_` (regressors) or before the file extension (time series).
- **Steps 4-5:** each `.mat` file has a 2-D variable `ts_data` (time points × ROIs).
- **Step 16:** each CSV has `TwinPairID`, `ZygosityGT1` (starting with `MZ` or `DZ`) and one column per state (`State 1`, `state_2`, ...).
- **Step 15:** the workbook has sheets `<phenotype>_MZ` and `<phenotype>_DZ`, each with columns `Twin1` and `Twin2`.

## Known issues

These were found while cleaning up the code. They are documented rather than changed, because fixing them could change published results.

1. **Step 3 may drop IC 0 and shift all indices by one.** The script reads the CSVs with `skiprows=N, header=0`, which skips N rows and then treats the next row as the header. Steps 1 and 2 already write their header rows, so the first data row (IC 0) is consumed as the header. Running steps 1-3 on synthetic data gave 299 classified components instead of 300. Check this against your real files.
2. **Step 2 sampling frequency.** The original comment said `fs` is `1/0.72` (HCP TR), but the value used was `0.5`. The default is kept at 0.5 so results match the earlier analysis. Pass `--fs` to change it. The 0.1 Hz high-frequency cut-off depends on this value.
3. **Step 9 uses 300 ROIs in a 64-ROI pipeline.** The original script defined `nROIs = 300` but never used it. The variable was removed. The clustering input is the 64-ROI matrix.
4. **Step 6 refers to a script that is not in the repository.** The prototype generation step (`X_prototypes`) and the publication-figure script are not included.

## License

MIT. See [LICENSE](LICENSE).

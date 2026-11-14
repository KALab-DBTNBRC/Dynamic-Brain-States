"""Step 2: motion correlation and high-frequency content of each IC time course.

For every subject present in both input directories, computes per IC:
  * the maximum (by absolute value) Pearson correlation with the motion regressors
  * the fraction of Welch spectral power above --hf-threshold (Hz)

Inputs:
    --time-series-dir   one "<subject>.txt" per subject; n_timepoints x n_components values
    --regressors-dir    one "<subject>_*.xlsx" per subject; first 12 columns = motion parameters

Outputs (written to --output-dir):
    group_component_features.csv        per subject and IC
    aggregated_component_features.csv   mean / std / median per IC across subjects

Example:
    python step02_motion_highfreq_features.py \
        --time-series-dir data/ic_timeseries \
        --regressors-dir data/regressors \
        --output-dir results/01_noise_classification
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import welch
from scipy.stats import pearsonr


def max_rp_correlation(time_courses, motion_params):
    """Per IC, the correlation with the motion parameter having the largest |r|."""
    max_corrs = []
    for tc in time_courses.T:
        corrs = [pearsonr(tc, mp)[0] for mp in motion_params.T]
        max_corrs.append(max(corrs, key=abs))
    return np.array(max_corrs)


def high_frequency_content(time_courses, fs=1.0, threshold=0.1):
    """Per IC, the fraction of Welch spectral power at frequencies above threshold."""
    ratios = []
    for tc in time_courses.T:
        f, pxx = welch(tc, fs=fs)
        ratios.append(np.sum(pxx[f > threshold]) / np.sum(pxx))
    return np.array(ratios)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--time-series-dir", required=True, type=Path, help="Directory with <subject>.txt IC time courses")
    p.add_argument("--regressors-dir", required=True, type=Path, help="Directory with <subject>_*.xlsx motion regressors")
    p.add_argument("--output-dir", required=True, type=Path, help="Directory for the output CSV files")
    p.add_argument("--n-timepoints", type=int, default=4800, help="Rows per subject time-course file (default: 4800)")
    p.add_argument("--n-components", type=int, default=300, help="Number of ICs (default: 300)")
    p.add_argument("--fs", type=float, default=0.5,
                   help="Sampling frequency in Hz, i.e. 1/TR (default: 0.5, the value used in the original analysis)")
    p.add_argument("--hf-threshold", type=float, default=0.1, help="High-frequency cut-off in Hz (default: 0.1)")
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    excel_files = sorted(f.name for f in args.regressors_dir.glob("*.xlsx"))
    time_series_files = sorted(f.name for f in args.time_series_dir.glob("*.txt"))

    # Subject IDs: text before the first underscore (regressors) / before the extension (time series)
    subjects_from_excel = {f.split("_")[0] for f in excel_files}
    subjects_from_time_series = {f.split(".")[0] for f in time_series_files}
    common_subjects = sorted(subjects_from_excel & subjects_from_time_series)
    print(f"Common subjects: {len(common_subjects)}")

    expected_size = args.n_timepoints * args.n_components
    results_list = []

    for subj in common_subjects:
        time_courses = np.loadtxt(args.time_series_dir / f"{subj}.txt")
        if time_courses.size != expected_size:
            print(f"Unexpected size for subject {subj}: expected {expected_size}, got {time_courses.size}. Skipping.")
            continue
        time_courses = time_courses.reshape((args.n_timepoints, args.n_components))

        excel_file = next((f for f in excel_files if f.startswith(subj)), None)
        if not excel_file:
            print(f"No regressor file found for subject {subj}. Skipping.")
            continue

        # Motion parameters are assumed to be the first 12 columns
        motion_params = pd.read_excel(args.regressors_dir / excel_file).iloc[:, :12].values
        if motion_params.shape[0] != time_courses.shape[0]:
            print(f"Row mismatch for subject {subj}: {time_courses.shape[0]} time points vs "
                  f"{motion_params.shape[0]} motion rows. Skipping.")
            continue

        max_corrs = max_rp_correlation(time_courses, motion_params)
        hf_ratios = high_frequency_content(time_courses, fs=args.fs, threshold=args.hf_threshold)

        for ic_index in range(time_courses.shape[1]):
            results_list.append({
                "Subject": subj,
                "IC Index": ic_index,
                "Max RP Correlation": max_corrs[ic_index],
                "High-Frequency Content": hf_ratios[ic_index],
            })
        print(f"Processed subject {subj}")

    if not results_list:
        raise SystemExit("No subjects were processed. Check the input directories and file naming.")

    results_df = pd.DataFrame(results_list)
    individual_path = args.output_dir / "group_component_features.csv"
    results_df.to_csv(individual_path, index=False)
    print(f"Individual results saved to {individual_path}")

    aggregated = results_df.groupby("IC Index").agg({
        "Max RP Correlation": ["mean", "std", "median"],
        "High-Frequency Content": ["mean", "std", "median"],
    }).reset_index()
    aggregated_path = args.output_dir / "aggregated_component_features.csv"
    aggregated.to_csv(aggregated_path, index=False)
    print(f"Aggregated results saved to {aggregated_path}")


if __name__ == "__main__":
    main()

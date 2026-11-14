"""Step 10: bar plot of every brain-state centroid (one PNG per state).

Input:  .mat file containing C (n_states x n_rois), e.g. 4 x 64.
Output: <output-dir>/centroid_state_<k>.png

Example:
    python step10_plot_centroids.py \
        --centroids-mat results/centroids.mat \
        --output-dir results/centroid_plots
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # write files only; no display needed
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--centroids-mat", required=True, type=Path, help=".mat file containing the centroid matrix C")
    p.add_argument("--output-dir", required=True, type=Path, help="Directory for the PNG files")
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    C = loadmat(args.centroids_mat)["C"]  # shape = (n_states, n_rois)
    print("Centroids shape:", C.shape)

    for i in range(C.shape[0]):
        fig, ax = plt.subplots(figsize=(12, 4))
        ax.bar(np.arange(1, C.shape[1] + 1), C[i, :])
        ax.set_title(f"Brain State {i+1} - ROI Activation Pattern")
        ax.set_xlabel("ROI Index")
        ax.set_ylabel("Mean Activation")
        fig.savefig(args.output_dir / f"centroid_state_{i+1}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    print(f"Saved {C.shape[0]} centroid plots to {args.output_dir}")


if __name__ == "__main__":
    main()

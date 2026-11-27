"""Step 1: fraction of each group-ICA component that lies in CSF and at the brain edge.

For every component map, the fraction is sum(map * mask) / sum(map). The CSF and edge
masks are resampled (linear) to the grid of the ICA maps first.

Outputs (written to --output-dir):
    edge_fractions.csv
    csf_fractions.csv

Example:
    python step01_calculate_fractions.py \
        --ica-maps melodic_IC300_sum.nii.gz \
        --csf-mask brain_extracted_csf_mask.nii.gz \
        --edge-mask outer_edge_mask.nii.gz \
        --output-dir results/01_noise_classification
"""
import argparse
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy.ndimage import zoom


def resample_to_shape(mask, target_shape):
    """Resample a 3D mask to target_shape using linear interpolation."""
    zoom_factors = [float(t) / s for t, s in zip(target_shape, mask.shape)]
    return zoom(mask, zoom_factors, order=1)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ica-maps", required=True, type=Path, help="4D NIfTI with the group ICA spatial maps")
    p.add_argument("--csf-mask", required=True, type=Path, help="NIfTI CSF mask")
    p.add_argument("--edge-mask", required=True, type=Path, help="NIfTI outer-edge mask")
    p.add_argument("--output-dir", required=True, type=Path, help="Directory for the output CSV files")
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    ica_maps = nib.load(args.ica_maps).get_fdata()
    csf_mask = nib.load(args.csf_mask).get_fdata()
    edge_mask = nib.load(args.edge_mask).get_fdata()

    # Resample masks to match the shape of the ICA maps
    target_shape = ica_maps.shape[:3]
    csf_mask_resampled = resample_to_shape(csf_mask, target_shape)
    edge_mask_resampled = resample_to_shape(edge_mask, target_shape)

    edge_fractions = []
    csf_fractions = []
    for i in range(ica_maps.shape[-1]):
        ic_map = ica_maps[..., i]
        edge_fractions.append(np.sum(ic_map * edge_mask_resampled) / np.sum(ic_map))
        csf_fractions.append(np.sum(ic_map * csf_mask_resampled) / np.sum(ic_map))

    pd.Series(edge_fractions, name="Edge Fraction").to_csv(args.output_dir / "edge_fractions.csv", index=False)
    pd.Series(csf_fractions, name="CSF Fraction").to_csv(args.output_dir / "csf_fractions.csv", index=False)
    print(f"Saved edge_fractions.csv and csf_fractions.csv to {args.output_dir}")


if __name__ == "__main__":
    main()

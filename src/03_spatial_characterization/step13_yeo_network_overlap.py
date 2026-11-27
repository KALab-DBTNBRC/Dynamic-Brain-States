"""Step 13: spatial correlation between each brain-state centroid and each Yeo network.

Builds a 4D image with one binary volume per Yeo label, then correlates (Pearson) every
Yeo network with every centroid across voxels.

Inputs:
    --yeo-atlas   Yeo 17-network atlas NIfTI
    --centroids   4D centroid NIfTI already resampled to the atlas grid (Step 11)

Outputs (written to --output-dir):
    ref_4d_img.nii.gz                      one binary volume per Yeo label
    correlation_matrix_with_p_values.csv   Reference_ID, Target_ID, Correlation, P-value

Example:
    python step13_yeo_network_overlap.py \
        --yeo-atlas Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz \
        --centroids centroids_resampled_to_Yeo.nii.gz \
        --output-dir results/yeo_overlap
"""
import argparse
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn.maskers import NiftiMasker
from scipy.stats import pearsonr


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--yeo-atlas", required=True, type=Path, help="Yeo 17-network atlas NIfTI")
    p.add_argument("--centroids", required=True, type=Path, help="Centroid NIfTI resampled to the atlas space")
    p.add_argument("--output-dir", required=True, type=Path, help="Directory for the outputs")
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    ref_img = nib.load(args.yeo_atlas)
    target_img = nib.load(args.centroids)
    assert ref_img.shape[:3] == target_img.shape[:3], (
        "Spatial dimensions of the images do not match. Resample the images to the same space and dimensions.")

    # One binary mask per label, stacked into a 4D image
    ref_data = ref_img.get_fdata()
    labels = np.unique(ref_data)
    labels = labels[labels != 0]  # exclude background
    ref_4d_data = np.squeeze(np.stack([(ref_data == label).astype(int) for label in labels], axis=-1))
    ref_4d_img = nib.Nifti1Image(ref_4d_data, ref_img.affine, ref_img.header)
    nib.save(ref_4d_img, args.output_dir / "ref_4d_img.nii.gz")
    print(f"Shape of ref_4d_img data: {ref_4d_img.shape}")  # (x, y, z, 17)

    masker = NiftiMasker()
    ref_voxels = masker.fit_transform(ref_4d_img).T     # (n_voxels, n_labels)
    target_voxels = masker.transform(target_img).T      # (n_voxels, n_states)
    print(f"Shape of ref_data: {ref_voxels.shape}")
    print(f"Shape of target_data: {target_voxels.shape}")

    rows = []
    for i in range(ref_voxels.shape[1]):
        for j in range(target_voxels.shape[1]):
            corr, p_val = pearsonr(ref_voxels[:, i], target_voxels[:, j])
            rows.append([f"ref_{i+1}", f"target_{j+1}", corr, round(p_val, 4)])

    df = pd.DataFrame(rows, columns=["Reference_ID", "Target_ID", "Correlation", "P-value"])
    out_csv = args.output_dir / "correlation_matrix_with_p_values.csv"
    df.to_csv(out_csv, index=False)
    print(f"Correlation matrix with p-values saved to {out_csv}")


if __name__ == "__main__":
    main()

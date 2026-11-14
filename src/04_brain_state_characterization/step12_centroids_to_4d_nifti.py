"""Step 12: project the brain-state centroids onto the IC spatial maps.

Each centroid is a vector of 64 IC weights. The weighted sum of the 64 IC maps gives one
brain map per state, saved as a 4D NIfTI (X x Y x Z x n_states).

Input:
    --centroids-mat   MATLAB v7.3 (HDF5) file containing C (n_states x 64)
    --ic-maps         4D NIfTI with the 64 final IC maps

Example:
    python step12_centroids_to_4d_nifti.py \
        --centroids-mat BrainStates_MedoidPartition.mat \
        --ic-maps melodic_IC_64_final.nii.gz \
        --output centroids_4D.nii.gz
"""
import argparse
import sys
from pathlib import Path

import h5py
import nibabel as nib
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--centroids-mat", required=True, type=Path, help="MATLAB v7.3 .mat file containing C")
    p.add_argument("--ic-maps", required=True, type=Path, help="4D NIfTI with the IC spatial maps")
    p.add_argument("--output", required=True, type=Path, help="Output 4D NIfTI path")
    p.add_argument("--expected-nodes", type=int, default=64, help="Expected number of ICs/nodes (default: 64)")
    return p.parse_args()


def load_centroids(mat_path):
    """Read C from a v7.3 MATLAB file, transposed to (n_states, n_nodes)."""
    print(f"Loading centroids from {mat_path}...")
    with h5py.File(mat_path, "r") as f:
        print("Keys/datasets in the file:", list(f.keys()))
        if "C" not in f:
            raise KeyError(f"Could not find 'C' at the top level of {mat_path}; found: {list(f.keys())}")
        # MATLAB is column-major, h5py row-major, hence the transpose
        return np.array(f["C"]).T


def main():
    args = parse_args()

    try:
        C = load_centroids(args.centroids_mat)
        K, n_nodes = C.shape

        if not args.ic_maps.exists():
            raise FileNotFoundError(f"Component map file not found: {args.ic_maps}")
        ic_maps_img = nib.load(args.ic_maps)
        ic_maps_data = ic_maps_img.get_fdata()

        n_ic_volumes = ic_maps_data.shape[-1]
        if n_nodes != args.expected_nodes or n_ic_volumes != args.expected_nodes:
            raise ValueError(f"Dimension mismatch: expected {args.expected_nodes} nodes. "
                             f"C has {n_nodes} columns, and the NIfTI has {n_ic_volumes} volumes.")
        print(f"Dimensions confirmed: {K} states x {n_nodes} nodes.")

        # (K x 64) @ (64 x voxels) = (K x voxels)
        n_voxels = np.prod(ic_maps_data.shape[:3])
        ic_maps_reshaped = ic_maps_data.reshape(n_voxels, n_ic_volumes)
        state_maps_2d = C @ ic_maps_reshaped.T

        new_shape = list(ic_maps_img.shape[:3]) + [K]
        state_maps_4d = state_maps_2d.T.reshape(new_shape)

        args.output.parent.mkdir(parents=True, exist_ok=True)
        nib.save(nib.Nifti1Image(state_maps_4d, ic_maps_img.affine), args.output)
        print(f"4D centroid image saved as: {args.output}")
    except Exception as e:
        sys.exit(f"An error occurred during conversion: {e}")


if __name__ == "__main__":
    main()

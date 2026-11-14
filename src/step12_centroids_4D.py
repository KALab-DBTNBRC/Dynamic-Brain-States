import nibabel as nib
import numpy as np
import os
import h5py

# --- 1. USER INPUTS (ADJUST PATHS IF NECESSARY) ---
MATLAB_FILE = 'IC_Labeling/IC_labelling/1st_Brain_State_Pipeline_Twins/BrainStates_MedoidPartition.mat'
IC_MAPS_FILE = 'IC_Labeling/IC_labelling/melodic_IC_64_final.nii.gz'
OUTPUT_FILE = 'centroids_4D.nii.gz'
EXPECTED_NODES = 64

# --- 2. LOAD DATA ---
try:
    print(f"Loading Centroids from {MATLAB_FILE}...")

    # Step 2a: Open v7.3 MATLAB file with h5py and locate the 'C' matrix
    with h5py.File(MATLAB_FILE, 'r') as f:
        print("Keys/datasets in the file:", list(f.keys()))
        # Try accessing 'C' directly; adjust the path if needed.
        try:
            C = np.array(f['C'])
        except KeyError:
            # If 'C' is nested, print out group structure and look for datasets
            for key in f.keys():
                print(f"Available sub-items in '{key}':", list(f[key].keys()))
            raise KeyError("Could not find 'C' at the expected top level. Check printed groups above.")
        # Transpose, as MATLAB (column-major) to Python (row-major) difference
        C = C.T

    # Ensure C is an array before getting shape
    if not isinstance(C, np.ndarray):
        C = np.array(C)

    K, n_nodes = C.shape

    # --- 2b: Load Component Maps (Reference image for spatial info) ---
    if not os.path.exists(IC_MAPS_FILE):
        raise FileNotFoundError(f"Error: Component map file not found at {IC_MAPS_FILE}. Please check the path.")

    ic_maps_img = nib.load(IC_MAPS_FILE)
    ic_maps_data = ic_maps_img.get_fdata()

    # --- 3. CRITICAL DIMENSION CHECK ---
    n_ic_volumes = ic_maps_data.shape[-1]

    if n_nodes != EXPECTED_NODES or n_ic_volumes != EXPECTED_NODES:
        raise ValueError(f"DIMENSION MISMATCH: Expected {EXPECTED_NODES} nodes. "
                         f"C has {n_nodes} columns, and the NIfTI has {n_ic_volumes} volumes.")

    print(f"Dimensions confirmed: {K} states x {n_nodes} nodes.")

    # --- 4. MATRIX MULTIPLICATION (Projection) ---
    n_voxels = np.prod(ic_maps_data.shape[:3])
    ic_maps_reshaped = ic_maps_data.reshape(n_voxels, n_ic_volumes)

    # Project the C matrix onto the Component Maps. (K x 64) @ (64 x Voxels) = (K x Voxels)
    state_maps_2d = C @ ic_maps_reshaped.T

    # --- 5. CONVERT TO 4D NIFTI ---
    # Reshape the K x Voxels matrix back into a 4D image (X x Y x Z x K)
    new_shape = list(ic_maps_img.shape[:3]) + [K]
    state_maps_4d = state_maps_2d.T.reshape(new_shape)

    # Create the new NIfTI image, preserving the spatial information (affine)
    centroids_4D_img = nib.Nifti1Image(state_maps_4d, ic_maps_img.affine)

    # --- 6. SAVE FILE ---
    nib.save(centroids_4D_img, OUTPUT_FILE)

    print(f"\n✅ Step 1 Complete. 4D Centroids image saved as: {OUTPUT_FILE}")

except Exception as e:
    print(f"An error occurred during conversion: {e}")

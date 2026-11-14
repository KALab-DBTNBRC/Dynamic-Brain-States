import nibabel as nib
from nilearn.image import resample_to_img
import os

# --- INPUT FILES ---
# Your 4D Centroids file (Output from Step 1)
ICA_CENTROIDS = 'IC_Labeling/IC_labelling/Step_1_C_to_4D_nifti/centroids_4D.nii.gz'

# The 1mm Yeo Atlas file path from your project structure
YEO_ATLAS = 'IC_Labeling/IC_labelling/atlas_RSNs/Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz'

# --- OUTPUT FILE ---
OUTPUT_FILE = 'centroids_resampled_to_Yeo.nii.gz'

print("Loading images and resampling centroids...")

if not os.path.exists(ICA_CENTROIDS):
    raise FileNotFoundError(f"Error: Centroids file not found at {ICA_CENTROIDS}. Did Step 1 complete?")
    
ica_img = nib.load(ICA_CENTROIDS)
mask_img = nib.load(YEO_ATLAS)

# Resample centroids to match the spatial characteristics of the Yeo Atlas
# 'linear' interpolation is appropriate for correlation maps.
ica_resampled = resample_to_img(ica_img, mask_img, interpolation='linear')

# Save the resampled image
nib.save(ica_resampled, OUTPUT_FILE)

print(f"✅ Resampling complete. File saved as: {OUTPUT_FILE}")
import os
import nibabel as nib
import numpy as np
from scipy.ndimage import zoom
import pandas as pd

# Function to resample the mask to match the target shape
def resample_to_shape(mask, target_shape):
    zoom_factors = [float(t) / s for t, s in zip(target_shape, mask.shape)]
    return zoom(mask, zoom_factors, order=1)

# Change to the desired directory
os.chdir(r"./src/Step_1_Noise and signal codes n results")

# Load ICA spatial maps
ica_maps = nib.load('melodic_IC300_sum.nii.gz').get_fdata()

# Load masks
#csf_mask = nib.load('csf_mask_pve_2.nii.gz').get_fdata()
csf_mask = nib.load('brain_extracted_csf_mask.nii.gz').get_fdata()
#edge_mask = nib.load('edge_mask.nii.gz').get_fdata()
edge_mask = nib.load('outer_edge_mask.nii.gz').get_fdata()

# Resample masks to match the shape of the ICA maps
target_shape = ica_maps.shape[:3]
csf_mask_resampled = resample_to_shape(csf_mask, target_shape)
edge_mask_resampled = resample_to_shape(edge_mask, target_shape)

# Initialize lists to store fractions
edge_fractions = []
csf_fractions = []

# Calculate fractions for each IC
for i in range(ica_maps.shape[-1]):
    ic_map = ica_maps[..., i]
    
    edge_fraction = np.sum(ic_map * edge_mask_resampled) / np.sum(ic_map)
    csf_fraction = np.sum(ic_map * csf_mask_resampled) / np.sum(ic_map)
    
    edge_fractions.append(edge_fraction)
    csf_fractions.append(csf_fraction)

# Convert lists to pandas Series for easy handling
edge_fractions = pd.Series(edge_fractions, name='Edge Fraction')
csf_fractions = pd.Series(csf_fractions, name='CSF Fraction')

# Save the results
edge_fractions.to_csv('edge_fractions.csv', index=False)
csf_fractions.to_csv('csf_fractions.csv', index=False)

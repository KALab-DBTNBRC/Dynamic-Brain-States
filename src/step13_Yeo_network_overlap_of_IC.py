import nibabel as nib
import numpy as np
import pandas as pd
from nilearn.input_data import NiftiMasker
from scipy.stats import pearsonr

# Load the reference label image and target image
ref_img_path = 'IC_Labeling/IC_labelling/atlas_RSNs/Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz'

# Please note that this image is resampled image to the atlas sapce. 
target_img_path = 'IC_Labeling/IC_labelling/Step2_Functional_Decoding/centroids_resampled_to_Yeo.nii.gz'

ref_img = nib.load(ref_img_path)
target_img = nib.load(target_img_path)

# Ensure the spatial dimensions are the same
assert ref_img.shape[:3] == target_img.shape[:3], "Spatial dimensions of the images do not match. Resample the images to the same space and dimensions."

# Extract each label as a separate image and concatenate them to form a 4D image
labels = np.unique(ref_img.get_fdata())
labels = labels[labels != 0]  # Exclude the background label (0)

# Initialize an empty list to hold the binary masks
label_images = []

# Create binary masks for each label
ref_data = ref_img.get_fdata()
for label in labels:
    label_mask = (ref_data == label).astype(int)
    label_images.append(label_mask)

# Stack the label images to form a 4D image
ref_4d_data = np.stack(label_images, axis=-1)

# Ensure the data is 4D (x, y, z, labels) by squeezing any extra dimensions
ref_4d_data = np.squeeze(ref_4d_data)
ref_4d_img = nib.Nifti1Image(ref_4d_data, ref_img.affine, ref_img.header)

# Save the 4D reference image for verification (optional)
nib.save(ref_4d_img, 'ref_4d_img.nii.gz')

# Verify the dimensions of the 4D reference image
print(f"Shape of ref_4d_img data: {ref_4d_img.shape}")  # Expected: (x, y, z, 17)

# Use NiftiMasker to extract the time series for each voxel
masker = NiftiMasker()
ref_data = masker.fit_transform(ref_4d_img).T  # Transpose to get shape (n_voxels, 17)
target_data = masker.transform(target_img).T  # Transpose to get shape (n_voxels, 100)

# Verify the shapes
print(f"Shape of ref_data: {ref_data.shape}")  # Expected: (n_voxels, 17)
print(f"Shape of target_data: {target_data.shape}")  # Expected: (n_voxels, 100)

# Initialize the correlation and p-value matrices with correct dimensions (17, 100)
correlations = np.zeros((ref_data.shape[1], target_data.shape[1]))
p_values = np.zeros((ref_data.shape[1], target_data.shape[1]))

# Compute the correlation and p-values
for i in range(ref_data.shape[1]):  # Iterate over each volume of the reference
    for j in range(target_data.shape[1]):  # Iterate over each volume of the target
        corr, p_val = pearsonr(ref_data[:, i], target_data[:, j])
        correlations[i, j] = corr
        p_values[i, j] = p_val

# Prepare data for saving to CSV
data = []
for i in range(ref_data.shape[1]):
    for j in range(target_data.shape[1]):
        data.append([f'ref_{i+1}', f'target_{j+1}', correlations[i, j], round(p_values[i, j], 4)])

# Create a DataFrame
df = pd.DataFrame(data, columns=['Reference_ID', 'Target_ID', 'Correlation', 'P-value'])

# Save to CSV
df.to_csv('correlation_matrix_with_p_values.csv', index=False)

print("Correlation matrix with p-values saved to correlation_matrix_with_p_values.csv")

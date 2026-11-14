import os
import numpy as np
import pandas as pd
from scipy.signal import welch
from scipy.stats import pearsonr

# Function to calculate the maximum RP correlation
def max_rp_correlation(time_courses, motion_params):
    max_corrs = []
    for tc in time_courses.T:
        corrs = [pearsonr(tc, mp)[0] for mp in motion_params.T]
        max_corrs.append(max(corrs, key=abs))
    return np.array(max_corrs)

# Function to calculate high-frequency content
def high_frequency_content(time_courses, fs=1.0, threshold=0.1):
    high_freq_ratios = []
    for tc in time_courses.T:
        f, Pxx = welch(tc, fs=fs)
        high_freq_power = np.sum(Pxx[f > threshold])
        total_power = np.sum(Pxx)
        high_freq_ratios.append(high_freq_power / total_power)
    return np.array(high_freq_ratios)

# Define the directories
time_series_dir = r'path_to_your_time_series_directory'  # Replace with the actual path to your time series directory
excel_dir = r'path_to_your_excel_directory'  # Replace with the actual path to your Excel directory

# List all Excel and .txt files in their respective directories
excel_files = [f for f in os.listdir(excel_dir) if f.endswith('.xlsx')]
time_series_files = [f for f in os.listdir(time_series_dir) if f.endswith('.txt')]

# Extract subject IDs from filenames (before the first underscore)
subjects_from_excel = {f.split('_')[0] for f in excel_files}
subjects_from_time_series = {f.split('.')[0] for f in time_series_files}

# Find common subjects
common_subjects = subjects_from_excel.intersection(subjects_from_time_series)
print(f"Common subjects: {common_subjects}")

# Placeholder for aggregated results
results_list = []

# Sampling frequency (fs) based on your TR (Repetition Time)
fs = 0.5  # Replace with the actual sampling frequency of your data (1/0.72)

# Iterate over common subjects
for subj in common_subjects:
    # Load subject-specific IC time series
    subject_file = os.path.join(time_series_dir, f'{subj}.txt')
    time_courses = np.loadtxt(subject_file)
    
    # Reshape the time courses to 4800 rows and 300 columns
    if time_courses.size != 4800 * 300:
        print(f"Unexpected size for time courses of subject {subj}. Expected {4800 * 300}, but got {time_courses.size}.")
        continue
    
    time_courses = time_courses.reshape((4800, 300))
    print(f"Time courses shape for subject {subj}: {time_courses.shape}")
    
    # Find the corresponding Excel file
    excel_file = next((f for f in excel_files if f.startswith(subj)), None)
    if not excel_file:
        print(f"No corresponding Excel file found for subject {subj}.")
        continue
    
    # Load subject-specific motion parameters from the Excel file
    excel_path = os.path.join(excel_dir, excel_file)
    subj_regressors = pd.read_excel(excel_path)
    
    # Extract motion parameters (assumes they are the first 12 columns)
    motion_params = subj_regressors.iloc[:, :12].values
    print(f"Motion parameters shape for subject {subj}: {motion_params.shape}")
    
    # Check if the number of rows match
    print(f"Number of rows in time courses for subject {subj}: {time_courses.shape[0]}")
    print(f"Number of rows in motion parameters for subject {subj}: {motion_params.shape[0]}")
    
    if motion_params.shape[0] != time_courses.shape[0]:
        print(f"Mismatch in rows between time courses and motion parameters for subject {subj}.")
        continue
    
    # Calculate maximum RP correlation
    max_corrs = max_rp_correlation(time_courses, motion_params)
    
    # Calculate high-frequency content
    high_freq_ratios = high_frequency_content(time_courses, fs=fs)
    
    for ic_index in range(time_courses.shape[1]):
        results_list.append({
            'Subject': subj,
            'IC Index': ic_index,
            'Max RP Correlation': max_corrs[ic_index],
            'High-Frequency Content': high_freq_ratios[ic_index]
        })

# Convert results to DataFrame and save
results_df = pd.DataFrame(results_list)


# Save the individual results
individual_results_path = r'path_to_your_output_directory\group_component_features.csv'  # Replace with the actual path
try:
    results_df.to_csv(individual_results_path, index=False)
    print(f"Individual results saved to {individual_results_path}")
except PermissionError as e:
    print(f"Permission error while saving individual results: {e}")
    print(f"Please ensure the file is not open and you have write permissions to the directory.")

# Aggregate the results across all subjects
aggregated_results = results_df.groupby('IC Index').agg({
    'Max RP Correlation': ['mean', 'std', 'median'],
    'High-Frequency Content': ['mean', 'std', 'median']
}).reset_index()

# Save the aggregated results
aggregated_results_path = r'path_to_your_output_directory\aggregated_component_features.csv'  # Replace with the actual path
try:
    aggregated_results.to_csv(aggregated_results_path, index=False)
    print(f"Aggregated results saved to {aggregated_results_path}")
except PermissionError as e:
    print(f"Permission error while saving aggregated results: {e}")
    print(f"Please ensure the file is not open and you have write permissions to the directory.")

# Determine if ICs are noise based on aggregated values
# This is a placeholder for your specific noise determination logic
#noise_thresholds = {
#    'Max RP Correlation Mean': 0.5,  # example threshold
#    'High-Frequency Content Mean': 0.5  # example threshold
#}

#aggregated_results['Is Noise'] = (
#    (aggregated_results[('Max RP Correlation', 'mean')] > noise_thresholds['Max RP Correlation Mean']) |
#    (aggregated_results[('High-Frequency Content', 'mean')] > noise_thresholds['High-Frequency Content Mean'])
#)

#Save the results with noise determination
#noise_determined_path = r'C:\Users\Tanya\Desktop\aggregated_component_features_with_noise.csv'
#try:
#    aggregated_results.to_csv(noise_determined_path, index=False)
#    print(f"Results with noise determination saved to {noise_determined_path}")
#except PermissionError as e:
#    print(f"Permission error while saving noise determination results: {e}")
#    print(f"Please ensure the file is not open and you have write permissions to the directory.")

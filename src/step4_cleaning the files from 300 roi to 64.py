import os
import scipy.io as sio
import numpy as np

# === CONFIGURATION ===
input_folder = "add_twin_matlab_files_here"  # <-- change this to the folder containing your .mat files
output_folder = "./src/Cleaned_twins"

columns_to_remove_within_71 = [9, 12, 27, 43, 46, 58, 69]  # Example: columns to remove (1-based indexing)
start_remove_from = 72
end_remove_to = 300
data_key = "ts_data"  # <-- change this to the key inside your .mat file (e.g., 'data', 'X', etc.)

# === SCRIPT ===
os.makedirs(output_folder, exist_ok=True)

for filename in os.listdir(input_folder):
    if filename.endswith(".mat"):
        file_path = os.path.join(input_folder, filename)
        mat_data = sio.loadmat(file_path)

        # Extract the main array
        if data_key not in mat_data:
            print(f"⚠️ Skipping {filename} — key '{data_key}' not found.")
            continue

        data = mat_data[data_key]

        # Ensure it's at least 2D
        if data.ndim < 2:
            print(f"⚠️ Skipping {filename} — not 2D data.")
            continue

        # Get all column indices
        total_cols = data.shape[1]

        # Columns to remove (convert 1-based to 0-based)
        remove_cols = [i - 1 for i in columns_to_remove_within_71 if i <= total_cols]
        remove_cols += list(range(start_remove_from - 1, min(end_remove_to, total_cols)))

        # Keep only the columns we want
        keep_cols = [i for i in range(total_cols) if i not in remove_cols]
        cleaned_data = data[:, keep_cols]

        # Replace old data with cleaned data
        mat_data[data_key] = cleaned_data

        # Save cleaned file
        output_path = os.path.join(output_folder, filename)
        sio.savemat(output_path, mat_data)

        print(f"✅ Cleaned {filename}: {data.shape[1]} → {cleaned_data.shape[1]} columns")

print(" All .mat files cleaned and saved successfully!")
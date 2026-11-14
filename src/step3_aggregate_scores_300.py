import pandas as pd
import os

# Define the directories
data_dir = r"./src/Step_1_Noise and signal codes n results"
output_dir = r"./src/Step_1_Noise and signal codes n results"

# Function to load CSV file and handle errors
def load_csv(file_path, skip_rows=0, header=None):
    try:
        data = pd.read_csv(file_path, skiprows=skip_rows, header=header)
        print(f"Successfully loaded {file_path}")
        return data
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None
    except pd.errors.EmptyDataError:
        print(f"Empty data: {file_path}")
        return None
    except pd.errors.ParserError:
        print(f"Parsing error: {file_path}")
        return None
    except Exception as e:
        print(f"An error occurred while loading {file_path}: {e}")
        return None

# Load the CSV files without skipping rows to inspect column names
agg_component_features_path = os.path.join(data_dir, 'aggregated_component_features.csv')
edge_fractions_path = os.path.join(data_dir, 'edge_fractions.csv')
csf_fractions_path = os.path.join(data_dir, 'csf_fractions.csv')

agg_component_features_raw = load_csv(agg_component_features_path)  # Load without skipping rows
edge_fractions_raw = load_csv(edge_fractions_path)  # Load without skipping rows
csf_fractions_raw = load_csv(csf_fractions_path)  # Load without skipping rows

# Print column names for inspection
print("Columns in aggregated_component_features_raw:")
print(agg_component_features_raw.columns)

print("Columns in edge_fractions_raw:")
print(edge_fractions_raw.columns)

print("Columns in csf_fractions_raw:")
print(csf_fractions_raw.columns)

# Load the CSV files again, skipping the appropriate rows for actual data
agg_component_features = load_csv(agg_component_features_path, skip_rows=2, header=0)  # Skip first 2 rows for data
edge_fractions = load_csv(edge_fractions_path, skip_rows=1, header=0)  # Skip first row for data
csf_fractions = load_csv(csf_fractions_path, skip_rows=1, header=0)  # Skip first row for data

# Rename columns in edge_fractions and csf_fractions for merging
edge_fractions.columns = ['Edge Fraction']
csf_fractions.columns = ['CSF Fraction']

# Add IC Index to edge_fractions and csf_fractions for merging
edge_fractions['IC Index'] = range(len(edge_fractions))
csf_fractions['IC Index'] = range(len(csf_fractions))

# Extract relevant columns from aggregated_component_features
agg_component_features.columns = ['IC Index', 'Max RP Correlation mean', 'Max RP Correlation std', 'Max RP Correlation median', 'High-Frequency Content mean', 'High-Frequency Content std', 'High-Frequency Content median']
features_df = agg_component_features[['IC Index', 'Max RP Correlation mean', 'High-Frequency Content mean']]
features_df = features_df.rename(columns={'Max RP Correlation mean': 'Max RP Correlation', 'High-Frequency Content mean': 'High-Frequency Content'})

# Verify the structure of the edge_fractions and csf_fractions dataframes
if 'Edge Fraction' not in edge_fractions.columns or 'CSF Fraction' not in csf_fractions.columns:
    print("Expected columns are not as expected in one or more data files. Exiting.")
    exit()


# Check if the first 100 rows are present
if len(edge_fractions) <= 100 or len(csf_fractions) <= 100:
    print("One or more data files do not contain at least 100 rows. Exiting.")
    exit()

# Merge Edge Fraction and CSF Fraction to features_df
features_df = features_df.merge(edge_fractions[['IC Index', 'Edge Fraction']], on='IC Index')
features_df = features_df.merge(csf_fractions[['IC Index', 'CSF Fraction']], on='IC Index')

# Define thresholds for each feature
rp_threshold = 0.5
hf_threshold = 0.3
edge_threshold = 0.4
csf_threshold = 0.3

# Classify each IC as signal or noise
def classify_ic(row):
    if (row['Max RP Correlation'] > rp_threshold or
        row['High-Frequency Content'] > hf_threshold or
        row['Edge Fraction'] > edge_threshold or
        row['CSF Fraction'] > csf_threshold):
        return 'Noise'
    else:
        return 'Signal'

features_df['Classification'] = features_df.apply(classify_ic, axis=1)

# Save classification results
output_path = os.path.join(output_dir, 'classification_results.csv')
features_df.to_csv(output_path, index=False)

print("Classification results saved to:", output_path)

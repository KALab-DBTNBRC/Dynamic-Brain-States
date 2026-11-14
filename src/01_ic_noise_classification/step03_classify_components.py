"""Step 3: classify each IC as Signal or Noise from the Step 1 and Step 2 outputs.

An IC is labelled Noise if ANY of these holds:
    Max RP Correlation     > 0.5
    High-Frequency Content > 0.3
    Edge Fraction          > 0.4
    CSF Fraction           > 0.3

Inputs (all in --input-dir):
    aggregated_component_features.csv   (Step 2)
    edge_fractions.csv, csf_fractions.csv   (Step 1)

Output (written to --output-dir):
    classification_results.csv

Example:
    python step03_classify_components.py \
        --input-dir results/01_noise_classification \
        --output-dir results/01_noise_classification
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

RP_THRESHOLD = 0.5
HF_THRESHOLD = 0.3
EDGE_THRESHOLD = 0.4
CSF_THRESHOLD = 0.3


def load_csv(file_path, skip_rows=0, header=None):
    """Read a CSV, returning None (with a message) on failure."""
    try:
        data = pd.read_csv(file_path, skiprows=skip_rows, header=header)
        print(f"Successfully loaded {file_path}")
        return data
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except pd.errors.EmptyDataError:
        print(f"Empty data: {file_path}")
    except pd.errors.ParserError:
        print(f"Parsing error: {file_path}")
    except Exception as e:
        print(f"An error occurred while loading {file_path}: {e}")
    return None


def classify_ic(row):
    if (row["Max RP Correlation"] > RP_THRESHOLD or
            row["High-Frequency Content"] > HF_THRESHOLD or
            row["Edge Fraction"] > EDGE_THRESHOLD or
            row["CSF Fraction"] > CSF_THRESHOLD):
        return "Noise"
    return "Signal"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True, type=Path, help="Directory with the Step 1 and Step 2 CSV files")
    p.add_argument("--output-dir", required=True, type=Path, help="Directory for classification_results.csv")
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    agg_path = args.input_dir / "aggregated_component_features.csv"
    edge_path = args.input_dir / "edge_fractions.csv"
    csf_path = args.input_dir / "csf_fractions.csv"

    # NOTE (unverified, kept from the original analysis): skiprows=N together with header=0 skips N
    # rows and then treats the next row as the header. Steps 1 and 2 write 2 (aggregate) and 1
    # (edge/CSF) header rows, so this consumes the first data row (IC 0) as the header and shifts
    # the IC index by one. Check against your real CSVs before relying on the classification.
    agg = load_csv(agg_path, skip_rows=2, header=0)
    edge = load_csv(edge_path, skip_rows=1, header=0)
    csf = load_csv(csf_path, skip_rows=1, header=0)
    if agg is None or edge is None or csf is None:
        sys.exit("Could not load all input files.")

    edge.columns = ["Edge Fraction"]
    csf.columns = ["CSF Fraction"]
    edge["IC Index"] = range(len(edge))
    csf["IC Index"] = range(len(csf))

    agg.columns = ["IC Index",
                   "Max RP Correlation mean", "Max RP Correlation std", "Max RP Correlation median",
                   "High-Frequency Content mean", "High-Frequency Content std", "High-Frequency Content median"]
    features = agg[["IC Index", "Max RP Correlation mean", "High-Frequency Content mean"]].rename(columns={
        "Max RP Correlation mean": "Max RP Correlation",
        "High-Frequency Content mean": "High-Frequency Content",
    })

    if len(edge) <= 100 or len(csf) <= 100:
        sys.exit("One or more data files contain 100 rows or fewer.")

    features = features.merge(edge[["IC Index", "Edge Fraction"]], on="IC Index")
    features = features.merge(csf[["IC Index", "CSF Fraction"]], on="IC Index")
    features["Classification"] = features.apply(classify_ic, axis=1)

    output_path = args.output_dir / "classification_results.csv"
    features.to_csv(output_path, index=False)
    print("Classification results saved to:", output_path)


if __name__ == "__main__":
    main()

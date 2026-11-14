"""Step 4: reduce each subject's 300-ROI time series to the 64 retained ROIs.

Columns are removed from the matrix stored under --data-key in each .mat file:
  * a few columns inside the first 71 (1-based): 9, 12, 27, 43, 46, 58, 69
  * all columns from 72 to 300

Example:
    python step04_select_rois.py \
        --input-dir data/twin_timeseries_300 \
        --output-dir data/twin_timeseries_64
"""
import argparse
from pathlib import Path

import scipy.io as sio

COLUMNS_TO_REMOVE_WITHIN_71 = [9, 12, 27, 43, 46, 58, 69]  # 1-based
START_REMOVE_FROM = 72  # 1-based, inclusive
END_REMOVE_TO = 300     # 1-based, inclusive


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True, type=Path, help="Directory with the per-subject .mat files")
    p.add_argument("--output-dir", required=True, type=Path, help="Directory for the cleaned .mat files")
    p.add_argument("--data-key", default="ts_data", help="Variable name inside each .mat file (default: ts_data)")
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for file_path in sorted(args.input_dir.glob("*.mat")):
        # Drop MATLAB metadata keys (__header__, ...) so they are not re-written
        mat_data = {k: v for k, v in sio.loadmat(file_path).items() if not k.startswith("__")}

        if args.data_key not in mat_data:
            print(f"Skipping {file_path.name}: key '{args.data_key}' not found.")
            continue

        data = mat_data[args.data_key]
        if data.ndim < 2:
            print(f"Skipping {file_path.name}: not 2D data.")
            continue

        total_cols = data.shape[1]
        remove_cols = [i - 1 for i in COLUMNS_TO_REMOVE_WITHIN_71 if i <= total_cols]
        remove_cols += list(range(START_REMOVE_FROM - 1, min(END_REMOVE_TO, total_cols)))

        keep_cols = [i for i in range(total_cols) if i not in remove_cols]
        mat_data[args.data_key] = data[:, keep_cols]

        sio.savemat(args.output_dir / file_path.name, mat_data)
        print(f"Cleaned {file_path.name}: {total_cols} -> {len(keep_cols)} columns")

    print("Done.")


if __name__ == "__main__":
    main()

"""Step 11: resample the 4D centroid image to the grid of the Yeo 17-network atlas.

Linear interpolation is used because the centroid maps are continuous.

Example:
    python step11_resample_to_yeo.py \
        --centroids centroids_4D.nii.gz \
        --yeo-atlas Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz \
        --output centroids_resampled_to_Yeo.nii.gz
"""
import argparse
from pathlib import Path

import nibabel as nib
from nilearn.image import resample_to_img


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--centroids", required=True, type=Path, help="4D centroid NIfTI (output of Step 12)")
    p.add_argument("--yeo-atlas", required=True, type=Path, help="1 mm Yeo 17-network atlas NIfTI")
    p.add_argument("--output", required=True, type=Path, help="Output NIfTI path")
    return p.parse_args()


def main():
    args = parse_args()
    for path in (args.centroids, args.yeo_atlas):
        if not path.exists():
            raise FileNotFoundError(f"Input file not found: {path}")

    print("Loading images and resampling centroids...")
    ica_resampled = resample_to_img(nib.load(args.centroids), nib.load(args.yeo_atlas), interpolation="linear")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    nib.save(ica_resampled, args.output)
    print(f"Resampling complete. File saved as: {args.output}")


if __name__ == "__main__":
    main()

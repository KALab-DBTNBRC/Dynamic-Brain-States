"""Step 14: contribution of Harvard-Oxford atlas regions to each brain-state map.

For every component (volume) of the 4D centroid image, and for every atlas region, reports
the mean loading, voxel count, total contribution (mean x voxels) and the region's share of
the component's summed absolute contribution. One Excel sheet per component.

Inputs:
    --centroids           4D centroid NIfTI (output of Step 12)
    --cortical-atlas      HarvardOxford-cort-maxprob-thr25-2mm.nii.gz
    --subcortical-atlas   HarvardOxford-sub-maxprob-thr25-2mm.nii.gz

Outputs (written to --output-dir):
    ica_HO_cortical_LRsplit_all_components_new.xlsx   cortical regions, split left/right
    ica_HO_subcortical_all_components_new.xlsx        subcortical regions

Example:
    python step14_atlas_region_overlap.py \
        --centroids centroids_4D.nii.gz \
        --cortical-atlas HarvardOxford-cort-maxprob-thr25-2mm.nii.gz \
        --subcortical-atlas HarvardOxford-sub-maxprob-thr25-2mm.nii.gz \
        --output-dir results/atlas_overlap
"""
import argparse
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn.image import resample_to_img

cortical_labels = {
    1: "frontal pole", 2: "insular cortex", 3: "superior frontal gyrus", 4: "middle frontal gyrus",
    5: "inferior frontal gyrus, pars triangularis", 6: "inferior frontal gyrus, pars opercularis",
    7: "precentral gyrus", 8: "temporal pole", 9: "superior temporal gyrus, anterior division",
    10: "superior temporal gyrus, posterior division", 11: "middle temporal gyrus, anterior division",
    12: "middle temporal gyrus, posterior division", 13: "middle temporal gyrus, temporooccipital part",
    14: "inferior temporal gyrus, anterior division", 15: "inferior temporal gyrus, posterior division",
    16: "inferior temporal gyrus, temporooccipital part", 17: "postcentral gyrus",
    18: "superior parietal lobule", 19: "supramarginal gyrus, anterior division",
    20: "supramarginal gyrus, posterior division", 21: "angular gyrus",
    22: "lateral occipital cortex, superior division", 23: "lateral occipital cortex, inferior division",
    24: "intracalcarine cortex", 25: "frontal medial cortex", 26: "juxtapositional lobule cortex",
    27: "subcallosal cortex", 28: "paracingulate gyrus", 29: "cingulate gyrus, anterior division",
    30: "cingulate gyrus, posterior division", 31: "precuneous cortex", 32: "cuneal cortex",
    33: "frontal orbital cortex", 34: "parahippocampal gyrus, anterior division",
    35: "parahippocampal gyrus, posterior division", 36: "lingual gyrus",
    37: "temporal fusiform cortex, anterior division", 38: "temporal fusiform cortex, posterior division",
    39: "temporal occipital fusiform cortex", 40: "occipital fusiform gyrus",
    41: "frontal operculum cortex", 42: "central opercular cortex", 43: "parietal operculum cortex",
    44: "planum polare", 45: "heschls gyrus", 46: "planum temporale",
    47: "supracalcarine cortex", 48: "occipital pole"
}

subcortical_labels = {
    1: "left cerebral white matter", 2: "left cerebral cortex", 3: "left lateral ventricle",
    4: "left thalamus", 5: "left caudate", 6: "left putamen", 7: "left pallidum", 8: "brain-stem",
    9: "left hippocampus", 10: "left amygdala", 11: "left accumbens", 12: "right cerebral white matter",
    13: "right cerebral cortex", 14: "right lateral ventricle", 15: "right thalamus",
    16: "right caudate", 17: "right putamen", 18: "right pallidum", 19: "right hippocampus",
    20: "right amygdala", 21: "right accumbens"
}

def extract_to_excel_with_contribution_split_hemisphere(ica_map_path, atlas_path, label_dict, output_excel, split_hemisphere=False):
    ica_img = nib.load(ica_map_path)
    ica_data = ica_img.get_fdata()
    n_components = ica_data.shape[3]

    atlas_img = nib.load(atlas_path)
    if ica_img.shape[:3] != atlas_img.shape:
        print(f"Resampling {atlas_path} to match ICA map dimensions...")
        atlas_img = resample_to_img(atlas_img, ica_img.slicer[..., 0], interpolation='nearest')
    atlas_data = atlas_img.get_fdata().astype(int)

    # Prepare coordinate map if splitting by hemisphere
    if split_hemisphere:
        affine = atlas_img.affine
        shape = atlas_data.shape
        coords = np.indices(shape).reshape(3, -1)
        xyz = affine[:3,:3] @ coords + affine[:3,3:4]
        x_coords = xyz[0].reshape(shape)

    with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
        for comp_idx in range(n_components):
            print(f"Processing component {comp_idx + 1}/{n_components}...")
            comp_data = ica_data[..., comp_idx]
            rows = []

            for region_id, region_name in label_dict.items():
                region_mask = atlas_data == region_id

                if split_hemisphere:
                    left_mask = np.logical_and(region_mask, x_coords < 0)
                    right_mask = np.logical_and(region_mask, x_coords > 0)
                    
                    for hemi, hemi_mask in zip(['left', 'right'], [left_mask, right_mask]):
                        num_voxels = np.sum(hemi_mask)

                        if num_voxels > 0:
                            region_values = comp_data[hemi_mask]
                            mean_loading = np.mean(region_values)
                            total_contribution = mean_loading * num_voxels

                            rows.append({
                                'Component': comp_idx + 1,
                                'Region_ID': region_id,
                                'Region_Name': f"{hemi} {region_name}",
                                'Mean_ICA_Loading': mean_loading,
                                'Num_Voxels': num_voxels,
                                'Total_Contribution': total_contribution
                            })
                else:
                    num_voxels = np.sum(region_mask)

                    if num_voxels > 0:
                        region_values = comp_data[region_mask]
                        mean_loading = np.mean(region_values)
                        total_contribution = mean_loading * num_voxels

                        rows.append({
                            'Component': comp_idx + 1,
                            'Region_ID': region_id,
                            'Region_Name': region_name,
                            'Mean_ICA_Loading': mean_loading,
                            'Num_Voxels': num_voxels,
                            'Total_Contribution': total_contribution
                        })

            df = pd.DataFrame(rows)
            df['%_of_Component'] = 100 * df['Total_Contribution'].abs() / df['Total_Contribution'].abs().sum()

            sheet_name = f"Comp_{comp_idx + 1:03d}"
            df.to_excel(writer, sheet_name=sheet_name, index=False)
            print(f"Sheet added: {sheet_name}")

    print(f"\nExcel saved: {output_excel}")

def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--centroids", required=True, type=Path, help="4D centroid NIfTI")
    p.add_argument("--cortical-atlas", required=True, type=Path, help="Harvard-Oxford cortical atlas NIfTI")
    p.add_argument("--subcortical-atlas", required=True, type=Path, help="Harvard-Oxford subcortical atlas NIfTI")
    p.add_argument("--output-dir", required=True, type=Path, help="Directory for the Excel files")
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Cortical (with hemisphere split)
    extract_to_excel_with_contribution_split_hemisphere(
        args.centroids,
        args.cortical_atlas,
        cortical_labels,
        output_excel=args.output_dir / "ica_HO_cortical_LRsplit_all_components_new.xlsx",
        split_hemisphere=True,
    )

    # Subcortical (no need to split)
    extract_to_excel_with_contribution_split_hemisphere(
        args.centroids,
        args.subcortical_atlas,
        subcortical_labels,
        output_excel=args.output_dir / "ica_HO_subcortical_all_components_new.xlsx",
        split_hemisphere=False,
    )


if __name__ == "__main__":
    main()

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, pearsonr
import re

# ============================================================
# CONFIG
# ============================================================

DATA_DIR = r"C:\Users\Tanya\Desktop\Project_Files\Step_9_model\Revision-heritability"
METRICS = ["fo", "dt", "ar"]

# ============================================================
# FUNCTIONS
# ============================================================

def load_data(metric):
    return pd.read_csv(f"{DATA_DIR}\\{metric}.csv")


def fix_and_align_states(df):
    """
    Normalize state names:
    State 1 -> state1
    State_1 -> state1
    STATE1 -> state1
    """

    rename_map = {}

    for col in df.columns:

        col_l = str(col).lower().strip()

        if col_l.startswith("state"):

            nums = re.findall(r"\d+", col_l)

            if nums:
                rename_map[col] = f"state{int(nums[0])}"

    df = df.rename(columns=rename_map)

    state_cols = sorted(
        [c for c in df.columns if c.startswith("state")],
        key=lambda x: int(x.replace("state", ""))
    )

    return df, state_cols


def build_pair_dataframe(df, state, zygosity):

    pair_rows = []

    subset = df[
        df["ZygosityGT1"]
        .astype(str)
        .str.upper()
        .str.startswith(zygosity)
    ]

    for pair_id, grp in subset.groupby("TwinPairID"):

        if len(grp) != 2:
            continue

        grp = grp.reset_index(drop=True)

        pair_rows.append({
            "TwinPairID": pair_id,
            "Twin1": grp.loc[0, state],
            "Twin2": grp.loc[1, state]
        })

    return pd.DataFrame(pair_rows)


def statewise_within_pair_distance(pair_df):

    return np.abs(
        pair_df["Twin1"] -
        pair_df["Twin2"]
    ).values


def statewise_twin_correlation(pair_df):

    r, p = pearsonr(
        pair_df["Twin1"],
        pair_df["Twin2"]
    )

    return r, p


# ============================================================
# MAIN LOOP
# ============================================================

results = []

for metric in METRICS:

    print("\n" + "=" * 60)
    print(f"BRAIN DYNAMIC METRIC: {metric}")
    print("=" * 60)

    df_raw = load_data(metric)

    df, state_cols = fix_and_align_states(df_raw)

    for state in state_cols:

        print(f"\nSTATE: {state.upper()}")
        print("-" * 40)

        # -----------------------------------
        # Create MZ and DZ twin-pair tables
        # -----------------------------------

        mz_pairs = build_pair_dataframe(
            df,
            state,
            "MZ"
        )

        dz_pairs = build_pair_dataframe(
            df,
            state,
            "DZ"
        )

        if len(mz_pairs) < 2 or len(dz_pairs) < 2:

            print("Insufficient twin pairs")
            continue

        # -----------------------------------
        # Distances
        # -----------------------------------

        mz_dist = statewise_within_pair_distance(
            mz_pairs
        )

        dz_dist = statewise_within_pair_distance(
            dz_pairs
        )

        stat, p = mannwhitneyu(
            mz_dist,
            dz_dist,
            alternative="less"
        )

        print(
            f"MZ mean distance : "
            f"{mz_dist.mean():.4f}"
        )

        print(
            f"DZ mean distance : "
            f"{dz_dist.mean():.4f}"
        )

        print(
            f"P-value (MZ < DZ): "
            f"{p:.6f}"
        )

        # -----------------------------------
        # Heritability
        # -----------------------------------

        r_mz, p_mz = statewise_twin_correlation(
            mz_pairs
        )

        r_dz, p_dz = statewise_twin_correlation(
            dz_pairs
        )

        h2 = 2 * (r_mz - r_dz)

        print(f"r_MZ = {r_mz:.3f} (p = {p_mz:.6f})")
        print(f"r_DZ = {r_dz:.3f} (p = {p_dz:.6f})")
        print(f"h²   = {h2:.3f}")

        results.append({

            "Metric": metric,
            "State": state,

            "MZ_mean_distance":
                mz_dist.mean(),

            "DZ_mean_distance":
                dz_dist.mean(),

            "Distance_pvalue":
                p,

            "r_MZ":
                r_mz,

            "r_DZ":
                r_dz,

            "r_MZ_pvalue":
                p_mz,

            "r_DZ_pvalue":
                p_dz,

            "Heritability_h2":
                h2,

            "N_MZ_pairs":
                len(mz_pairs),

            "N_DZ_pairs":
                len(dz_pairs)
        })

# ============================================================
# SAVE RESULTS
# ============================================================

results_df = pd.DataFrame(results)

results_df.to_csv(
    r"C:\Users\Tanya\Desktop\Project_Files\Step_9_model\Revision-heritability\Twin_brain_dynamics_statewise_results.csv",
    index=False
)

print("\nSTATE-WISE RESULTS SAVED SUCCESSFULLY")
print(results_df)
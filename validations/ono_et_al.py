"""
Pipeline for analyzing Relative Performance Score (RPS) data from microbial coculture experiments (Ono et al.).

This script performs:
- Data loading and preprocessing from multiple supplementary tables.
- Correlation analysis between RPS and interaction mode indicators.
- Random forest regression to assess feature importance.
- Grouped analysis by carbon source count and interaction class.
- Permutation-based ANOVA and pairwise comparisons on RPS values.
- Generation of heatmaps and summary statistics.
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from statsmodels.stats.multitest import multipletests

# Suppress specific warnings if needed; here we ignore all for simplicity
warnings.filterwarnings("ignore")

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------
SUPPLEMENT_TABLE_S11 = "supplementary-tables_wraf224.xlsx"
SUPPLEMENT_TABLE_S6 = "supplementary-tables_wraf224.xlsx"
SUPPLEMENT_TABLE_S8 = "supplementary-tables_wraf224.xlsx"
CARBON_COMP_TSV = "carboncomp_output.tsv.xz"

# ------------------------------------------------------------------------------
# Utility Functions
# ------------------------------------------------------------------------------


def prepare_pairwise_data(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create summary statistics for each unique strain pair.

    For each unique pair (order-independent), computes:
    - median RPS
    - modal interaction class
    - number of observations
    - proportion of competition and mutualism

    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns 'genome1', 'genome2', 'RPS', 'interaction_class'.

    Returns
    -------
    pd.DataFrame
        Summary per unique pair.
    """
    # Create a unique pair ID (order does not matter)
    df = df.copy()
    df["pair_id"] = df.apply(
        lambda row: "_".join(sorted([str(row["genome1"]), str(row["genome2"])])), axis=1
    )

    pair_summary = []
    for pair_id, group in df.groupby("pair_id"):
        median_rps = group["RPS"].median()

        # Most frequent interaction class (mode, take first if tie)
        interaction_counts = group["interaction_class"].value_counts()
        if not interaction_counts.empty:
            max_count = interaction_counts.max()
            modes = interaction_counts[interaction_counts == max_count].index.tolist()
            modal_interaction = modes[0]
        else:
            modal_interaction = None

        pair_summary.append(
            {
                "pair_id": pair_id,
                "RPS": median_rps,
                "modal_interaction": modal_interaction,
                "n_observations": len(group),
                "competition_proportion": (group["interaction_class"] == "Competition").mean(),
                "mutualism_proportion": (group["interaction_class"] == "Mutualism").mean(),
            }
        )

    return pd.DataFrame(pair_summary)


def anova_f_statistic(groups: list[np.ndarray]) -> float:
    """
    Calculate the one-way ANOVA F-statistic for given groups.

    Parameters
    ----------
    groups : list of np.ndarray
        Each array contains data for one group.

    Returns
    -------
    float
        F-statistic, or 0.0 if insufficient groups or data.
    """
    # Remove empty groups
    groups = [g for g in groups if len(g) > 0]

    if len(groups) < 2:
        return 0.0

    all_data = np.concatenate(groups)
    grand_mean = np.mean(all_data)

    # Between-group sum of squares
    ssb = 0
    for group in groups:
        n = len(group)
        if n > 0:
            group_mean = np.mean(group)
            ssb += n * (group_mean - grand_mean) ** 2

    # Within-group sum of squares
    ssw = 0
    for group in groups:
        if len(group) > 0:
            group_mean = np.mean(group)
            ssw += np.sum((group - group_mean) ** 2)

    k = len(groups)
    n_total = len(all_data)
    df_between = k - 1
    df_within = n_total - k

    if df_within == 0 or ssw == 0:
        return 0.0

    ms_between = ssb / df_between
    ms_within = ssw / df_within
    return ms_between / ms_within


def permutation_anova_by_interaction_class(
    df: pd.DataFrame, n_permutations: int = 9999
) -> tuple[float | None, float | None, list]:
    """
    Perform a permutation-based ANOVA on RPS grouped by modal interaction class.

    Parameters
    ----------
    df : pd.DataFrame
        Contains 'genome1', 'genome2', 'RPS', 'interaction_class'.
    n_permutations : int
        Number of permutations.

    Returns
    -------
    f_observed : float or None
        Observed F-statistic.
    p_value : float or None
        Permutation p-value.
    interaction_classes : list
        List of interaction classes present.
    """
    pair_summary = prepare_pairwise_data(df)
    print(f"\nNumber of unique strain pairs: {len(pair_summary)}")
    print("\nPair summary (first few rows):")
    print(pair_summary[["pair_id", "modal_interaction", "RPS", "n_observations"]].head())

    # Group by modal interaction
    interaction_classes = sorted(pair_summary["modal_interaction"].unique())
    groups = []
    for ic in interaction_classes:
        group_data = pair_summary[pair_summary["modal_interaction"] == ic]["RPS"].values
        if len(group_data) > 0:
            groups.append(group_data)
            print(f"  {ic}: {len(group_data)} pairs")

    if len(groups) < 2:
        print("Need at least 2 groups with data for ANOVA.")
        return None, None, interaction_classes

    f_observed = anova_f_statistic(groups)
    print(f"\nObserved F-statistic: {f_observed:.4f}")

    # Permutation test
    all_rps = np.concatenate(groups)
    group_sizes = [len(g) for g in groups]
    f_permuted = []

    print(f"\nRunning {n_permutations} permutations...")
    for i in range(n_permutations):
        if i % 1000 == 0 and i > 0:
            print(f"  Completed {i} permutations")

        permuted_rps = np.random.permutation(all_rps)
        permuted_groups = []
        start = 0
        for size in group_sizes:
            permuted_groups.append(permuted_rps[start : start + size])
            start += size
        f_permuted.append(anova_f_statistic(permuted_groups))

    f_permuted = np.array(f_permuted)
    p_value = (np.sum(f_permuted >= f_observed) + 1) / (n_permutations + 1)
    return f_observed, p_value, interaction_classes


def permutation_test_pairwise_comparisons(
    df: pd.DataFrame, n_permutations: int = 9999
) -> list[dict]:
    """
    Perform pairwise permutation tests between interaction classes on RPS.

    Parameters
    ----------
    df : pd.DataFrame
        Contains 'genome1', 'genome2', 'RPS', 'interaction_class'.
    n_permutations : int
        Number of permutations.

    Returns
    -------
    list of dict
        Each dict contains comparison statistics.
    """
    pair_summary = prepare_pairwise_data(df)
    classes = pair_summary["modal_interaction"].unique()

    print("\nPairwise comparisons:")
    results = []

    for i in range(len(classes)):
        for j in range(i + 1, len(classes)):
            class1, class2 = classes[i], classes[j]
            data1 = pair_summary[pair_summary["modal_interaction"] == class1]["RPS"].values
            data2 = pair_summary[pair_summary["modal_interaction"] == class2]["RPS"].values

            if len(data1) < 2 or len(data2) < 2:
                continue

            obs_diff = np.mean(data1) - np.mean(data2)
            combined = np.concatenate([data1, data2])
            n1 = len(data1)

            perm_diffs = []
            for _ in range(n_permutations):
                perm = np.random.permutation(combined)
                perm_diff = np.mean(perm[:n1]) - np.mean(perm[n1:])
                perm_diffs.append(perm_diff)

            perm_diffs = np.array(perm_diffs)
            p_value = (np.sum(np.abs(perm_diffs) >= np.abs(obs_diff)) + 1) / (n_permutations + 1)

            results.append(
                {
                    "class1": class1,
                    "class2": class2,
                    "mean1": np.mean(data1),
                    "mean2": np.mean(data2),
                    "n1": len(data1),
                    "n2": len(data2),
                    "diff": obs_diff,
                    "p_value": p_value,
                }
            )

    results.sort(key=lambda x: x["p_value"])
    for r in results:
        sig = (
            "***"
            if r["p_value"] < 0.001
            else "**" if r["p_value"] < 0.01 else "*" if r["p_value"] < 0.05 else ""
        )
        print(
            f"  {r['class1']} (n={r['n1']}, mean={r['mean1']:.4f}) vs "
            f"{r['class2']} (n={r['n2']}, mean={r['mean2']:.4f}): "
            f"diff={r['diff']:.4f}, p={r['p_value']:.4f}{sig}"
        )

    return results


def run_simplified_analysis(df: pd.DataFrame) -> dict:
    """
    Run a simplified permutation-based analysis on RPS values.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain 'genome1', 'genome2', 'RPS', 'interaction_class'.

    Returns
    -------
    dict
        Summary of results including pair summary, ANOVA result, pairwise results.
    """
    print("=" * 60)
    print("SIMPLIFIED PERMUTATION ANALYSIS FOR RPS VALUES")
    print("=" * 60)

    print("\n1. BASIC STATISTICS")
    print("-" * 40)

    unique_pairs = df.apply(
        lambda row: tuple(sorted([str(row["genome1"]), str(row["genome2"])])), axis=1
    ).unique()
    print(f"Total unique strain pairs: {len(unique_pairs)}")
    print(f"Total observations: {len(df)}")

    overall_means = df.groupby("interaction_class")["RPS"].agg(["mean", "std", "count", "median"])
    print("\nOverall RPS by interaction class:")
    for idx, row in overall_means.iterrows():
        print(
            f"  {idx}: mean={row['mean']:.4f}, std={row['std']:.4f}, "
            f"median={row['median']:.4f}, n={int(row['count'])}"
        )

    print("\n2. PERMUTATION ANOVA (USING UNIQUE PAIRS)")
    print("-" * 40)
    f_obs, p_value, classes = permutation_anova_by_interaction_class(df, n_permutations=5000)

    if f_obs is not None:
        print(f"\nPermutation ANOVA Result:")
        print(f"  F = {f_obs:.4f}, p = {p_value:.4f}")
        if p_value < 0.05:
            print("  -> Significant difference between interaction classes (p < 0.05)")
        else:
            print("  -> No significant difference between interaction classes")

    print("\n3. PAIRWISE COMPARISONS")
    print("-" * 40)
    pairwise_results = permutation_test_pairwise_comparisons(df, n_permutations=5000)

    print("\n4. TREND ANALYSIS")
    print("-" * 40)

    pair_summary = prepare_pairwise_data(df)
    if len(pair_summary) > 2:
        rho, p = spearmanr(pair_summary["RPS"], pair_summary["competition_proportion"])
        print(f"\nCorrelation (pair level):")
        print(f"  RPS vs. Competition proportion: rho = {rho:.4f}, p = {p:.4f}")

    print("\n" + "=" * 60)
    print("KEY INTERPRETATION POINTS:")
    print("-" * 60)
    print("1. The permutation ANOVA tests if RPS differs between interaction classes")
    print("2. Each unique strain pair (n≈28) contributes one data point")
    print("3. Pairwise comparisons show which specific classes differ")
    print("4. The heatmap shows median trends despite statistical limitations")
    print("=" * 60)

    return {
        "pair_summary": pair_summary,
        "anova_result": {"f": f_obs, "p": p_value, "classes": classes},
        "pairwise_results": pairwise_results,
    }


# ------------------------------------------------------------------------------
# Data Loading Functions
# ------------------------------------------------------------------------------


def load_strain_name_map() -> dict:
    """
    Load strain name mapping from Table S11 (column 2 -> column 1).

    Returns
    -------
    dict
        Mapping from original strain ID to descriptive name.
    """
    df_map = pd.read_excel(SUPPLEMENT_TABLE_S11, sheet_name="Table S11")
    # Assuming column index 2 is the key, column index 1 is the value
    mapping = df_map.set_index(df_map.columns[2])[df_map.columns[1]].to_dict()
    return mapping


def load_interaction_modes() -> pd.DataFrame:
    """
    Load interaction mode data from Table S6 and make it symmetric (both directions).

    Returns
    -------
    pd.DataFrame
        Columns: 'genome1', 'genome2', plus mode indicator columns (e.g., 'Competition', etc.).
    """
    df_modes = pd.read_excel(SUPPLEMENT_TABLE_S6, sheet_name="Table S6")
    # Original table has a 'pair' column like "genome1:genome2"
    df_modes["genome1"] = df_modes["pair"].apply(lambda x: x.split(":")[0])
    df_modes["genome2"] = df_modes["pair"].apply(lambda x: x.split(":")[1])
    df_modes.drop("pair", axis=1, inplace=True)

    # Duplicate with swapped genomes to have both directions for merging
    df_modes_swapped = df_modes.copy()
    df_modes_swapped.rename(
        columns={"genome1": "genome2", "genome2": "genome1"}, inplace=True
    )
    df_modes = pd.concat([df_modes, df_modes_swapped], ignore_index=True)
    return df_modes


def load_carbon_interaction_data() -> pd.DataFrame:
    """
    Load carbon source count and interaction class from Table S8 and make it symmetric.

    Returns
    -------
    pd.DataFrame
        Columns: 'genome1', 'genome2', 'carbonsource_count', 'interaction_class'.
    """
    df_carbon = pd.read_excel(SUPPLEMENT_TABLE_S8, sheet_name="Table S8")
    df_carbon = df_carbon[
        ["coculture_species_GFP", "coculture_species_mScarlet", "carbonsource_count", "interaction_class"]
    ].copy()
    df_carbon.rename(
        columns={"coculture_species_GFP": "genome1", "coculture_species_mScarlet": "genome2"},
        inplace=True,
    )

    # Symmetric version
    df_carbon_swapped = df_carbon.copy()
    df_carbon_swapped.rename(
        columns={"genome1": "genome2", "genome2": "genome1"}, inplace=True
    )
    df_carbon = pd.concat([df_carbon, df_carbon_swapped], ignore_index=True)
    return df_carbon


def load_rps_data() -> pd.DataFrame:
    """
    Load the main RPS data from carboncomp_output.tsv.xz and rename the score column to 'RPS'.

    Returns
    -------
    pd.DataFrame
        Raw RPS data with columns: 'genome1', 'genome2', 'RPS'.
    """
    df_rps = pd.read_table(CARBON_COMP_TSV)
    # Rename the original 'EIT' column to 'RPS'
    df_rps.rename(columns={"EIT": "RPS"}, inplace=True)
    # Extract first two parts of the genome identifier (e.g., "GCF_000123" from "GCF_000123.1")
    df_rps["genome1"] = df_rps["genome1"].apply(lambda x: "_".join(x.split("_")[:2]))
    df_rps["genome2"] = df_rps["genome2"].apply(lambda x: "_".join(x.split("_")[:2]))
    return df_rps


# ------------------------------------------------------------------------------
# Analysis Functions
# ------------------------------------------------------------------------------


def correlate_rps_with_modes(df: pd.DataFrame, mode_columns: list) -> pd.DataFrame:
    """
    Compute Spearman correlations between RPS and each mode column.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain 'RPS' and the mode columns.
    mode_columns : list
        Names of mode indicator columns.

    Returns
    -------
    pd.DataFrame
        Results with columns: v1, v2, rho, p, q.
    """
    results = []
    for col in mode_columns:
        r, p = spearmanr(df["RPS"], df[col])
        results.append(("RPS", col, r, p))

    res_df = pd.DataFrame(results, columns=["v1", "v2", "rho", "p"])
    _, res_df["q"], _, _ = multipletests(res_df["p"], method="fdr_bh")
    return res_df


def random_forest_importance(df: pd.DataFrame, feature_cols: list, target: str = "RPS") -> pd.Series:
    """
    Fit a random forest regressor and return feature importances.

    Parameters
    ----------
    df : pd.DataFrame
        Data.
    feature_cols : list
        Feature column names.
    target : str
        Target column name.

    Returns
    -------
    pd.Series
        Feature importances indexed by feature name.
    """
    regr = RandomForestRegressor(max_depth=2, random_state=0)
    regr.fit(df[feature_cols], df[target])
    return pd.Series(regr.feature_importances_, index=feature_cols, name="variable importance")


def create_carbon_summary_and_heatmaps(
    df: pd.DataFrame, output_prefix: str = ""
) -> pd.DataFrame:
    """
    Group by carbon source count and interaction class, compute statistics,
    save to TSV, and generate heatmaps for mean and median RPS.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain 'carbonsource_count', 'interaction_class', 'RPS'.
    output_prefix : str
        Prefix for output files.

    Returns
    -------
    pd.DataFrame
        Summary statistics.
    """
    grouped = df.groupby(["carbonsource_count", "interaction_class"])["RPS"]
    agg = grouped.agg(["mean", "std", "median"]).reset_index()
    agg.columns = ["carbonsource_count", "interaction_class", "mean", "std", "median"]
    agg.to_csv(
        f"{output_prefix}csource_RPS_interaction_meanstdmedian.tsv",
        sep="\t",
        index=False,
    )

    # Heatmap for mean
    pivot_mean = df.pivot_table(
        index="carbonsource_count",
        columns="interaction_class",
        values="RPS",
        aggfunc="mean",
    )
    plt.figure()
    sns.heatmap(pivot_mean)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}heatmap_RPS_csources_interaction.svg")
    plt.close()

    # Heatmap for median
    pivot_median = df.pivot_table(
        index="carbonsource_count",
        columns="interaction_class",
        values="RPS",
        aggfunc="median",
    )
    plt.figure()
    sns.heatmap(pivot_median)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}heatmap_RPS_csources_interaction_median.svg")
    plt.close()

    return agg


# ------------------------------------------------------------------------------
# Main Pipeline
# ------------------------------------------------------------------------------


def main():
    """Execute the full analysis pipeline."""
    print("Starting RPS analysis pipeline...")

    # 1. Load strain name mapping and apply to RPS data
    strain_map = load_strain_name_map()
    df_rps = load_rps_data()
    df_rps["genome1"] = df_rps["genome1"].map(strain_map)
    df_rps["genome2"] = df_rps["genome2"].map(strain_map)
    # Drop rows where mapping failed (if any)
    df_rps.dropna(subset=["genome1", "genome2"], inplace=True)

    # 2. Load interaction modes and merge with RPS data
    df_modes = load_interaction_modes()
    df_merged = df_rps.merge(df_modes, on=["genome1", "genome2"], how="inner")
    mode_columns = df_merged.columns[-6:].tolist()  # last 6 columns are mode indicators

    # 3. Correlations
    corr_results = correlate_rps_with_modes(df_merged, mode_columns)
    corr_results.to_csv("correlations_mode_RPS.tsv", sep="\t", index=False)
    print("\nCorrelation results saved to correlations_mode_RPS.tsv")

    # 4. Random forest variable importance
    varimp = random_forest_importance(df_merged, mode_columns, target="RPS")
    print("\nRandom forest variable importance:")
    print(varimp)

    # 5. Create 'Any' column (sum of first 5 modes, excluding the last, presumably 'Neutral')
    df_merged["Any"] = df_merged[mode_columns[:-1]].sum(axis=1)
    r, p = spearmanr(df_merged["RPS"], df_merged["Any"])
    print(f"\nRPS vs. Any relation (except Neutral): rho={r:.4f}, p={p:.4f}")

    # 6. Load carbon source data, merge with RPS, and produce summaries/heatmaps
    df_carbon = load_carbon_interaction_data()
    df_carbon_rps = df_carbon.merge(
        df_rps[["genome1", "genome2", "RPS"]], on=["genome1", "genome2"], how="inner"
    )
    carbon_summary = create_carbon_summary_and_heatmaps(df_carbon_rps)
    print("\nCarbon source summary saved and heatmaps generated.")

    # 7. Run simplified permutation analysis on the carbon-merged data
    #    (which already contains RPS and interaction_class)
    analysis_results = run_simplified_analysis(df_carbon_rps)

    print("\nPipeline completed successfully.")


if __name__ == "__main__":
    main()

"""
Script for analyzing CaCo-derived metrics with perturbation experiment data.
Includes standard analyses (RPS distributions, statistics, correlations, RF)
and cross-correlation function (CCF) analysis between RPS and collection date
across experimental conditions (Coello-Camba et al. dataset).
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
from sklearn.ensemble import RandomForestRegressor
from statsmodels.tsa.stattools import ccf
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# ------------------------------------------------------------------------------
# Data loading and preparation (from coello-camba_et_all.py)
# ------------------------------------------------------------------------------
def load_and_prepare_data(carboncomp_path, metadata_path):
    """
    Load and merge genomic data with perturbation experiment metadata.

    Parameters
    ----------
    carboncomp_path : str
        Path to carboncomp_output.tsv.xz file
    metadata_path : str
        Path to metadata file

    Returns
    -------
    pandas.DataFrame
        Merged dataframe with genomic and environmental data
    """
    print("Loading and preparing perturbation experiment data...")

    # Load genomic data
    print("  Loading genomic data...")
    df = pd.read_table(carboncomp_path)

    # Load metadata
    print("  Loading metadata...")
    metadata = pd.read_table(metadata_path)

    # Create condition dictionary from metadata
    env_dict = metadata.set_index('genome')['Description'].to_dict()

    # Clean and prepare metadata
    metadata['sa1'] = metadata['genome'].apply(lambda x: x.split('_')[1])
    metadata = metadata[[
        'sa1', 'collection date', 'geographic location',
        'latitude and longitude', 'sample collection device or method',
        'sample material processing', 'sample size', 'ammonium',
        'miscellaneous parameter', 'nitrate', 'nitrite', 'phosphate',
        'silicate', 'chlorophyll', 'density', 'temperature', 'Description'
    ]].drop_duplicates().dropna()

    # Add condition information to genomic data
    print("  Merging condition information...")
    df['cond1'] = df['genome1'].apply(lambda x: env_dict.get(x.replace('.fa', '')))
    df['cond2'] = df['genome2'].apply(lambda x: env_dict.get(x.replace('.fa', '')))

    # Remove rows where conditions couldn't be determined
    df = df.dropna(subset=['cond1', 'cond2'])

    # Filter for pairs within the same condition
    print("  Filtering for same-condition pairs...")
    same_cond_df = df[df['cond1'] == df['cond2']].copy()

    # Add sample identifiers
    same_cond_df['sa1'] = same_cond_df['genome1'].apply(lambda x: x.split('_')[1])
    same_cond_df['sa2'] = same_cond_df['genome2'].apply(lambda x: x.split('_')[1])

    # Keep only pairs from the same sample
    same_cond_df = same_cond_df[same_cond_df['sa1'] == same_cond_df['sa2']]

    # Merge with environmental metadata
    print("  Merging environmental data...")
    same_cond_df = same_cond_df.merge(metadata, on='sa1', how='inner')

    # Convert collection date to integer format for analysis
    same_cond_df['collection date'] = pd.to_datetime(
        same_cond_df['collection date']
    ).dt.strftime('%Y%m%d').astype(int)

    print(f"  Final dataset: {same_cond_df.shape[0]} observations")
    print(f"  Conditions: {same_cond_df['cond1'].nunique()} unique conditions")

    return same_cond_df


# ------------------------------------------------------------------------------
# Standard statistical analyses (from coello-camba_et_all.py)
# ------------------------------------------------------------------------------
def calculate_summary_statistics(df, value_column='RPS', condition_column='cond1'):
    """
    Calculate summary statistics for RPS across conditions.

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe
    value_column : str
        Column name for the value to analyze
    condition_column : str
        Column name containing condition labels

    Returns
    -------
    tuple (pandas.DataFrame, list)
        Summary statistics per condition and the condition order list
    """
    print(f"\nCalculating summary statistics for {value_column}...")

    condition_order = [
        'N,P Day 1 single addition',
        'N,P 2 week continuous addition',
        'control mesocosm',
        'N,P,Si Day 1 single addition',
        'N,P,Si 2 week continuous addition'
    ]

    stats_list = []
    for condition in condition_order:
        if condition in df[condition_column].values:
            condition_data = df[df[condition_column] == condition][value_column]

            stats_dict = {
                'condition': condition,
                'mean': condition_data.mean(),
                'median': condition_data.median(),
                'std': condition_data.std(),
                'min': condition_data.min(),
                'max': condition_data.max(),
                'q25': condition_data.quantile(0.25),
                'q75': condition_data.quantile(0.75),
                'n_observations': len(condition_data),
                'n_unique_pairs': df[df[condition_column] == condition][['genome1', 'genome2']].drop_duplicates().shape[0]
            }
            stats_list.append(stats_dict)

    stats_df = pd.DataFrame(stats_list)
    stats_df.to_csv('RPS_summary_statistics_per_condition.tsv', sep='\t', index=False)
    print(f"  Summary statistics saved to 'RPS_summary_statistics_per_condition.tsv'")

    return stats_df, condition_order


def perform_statistical_tests(df, value_column='RPS', condition_column='cond1', condition_order=None):
    """
    Perform statistical tests between conditions.

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe
    value_column : str
        Column name for the value to analyze
    condition_column : str
        Column name containing condition labels
    condition_order : list, optional
        Order of conditions for comparisons

    Returns
    -------
    tuple
        (kruskal_results, mannwhitney_df, levene_df)
    """
    print(f"\nPerforming statistical tests for {value_column}...")

    if condition_order is None:
        condition_order = df[condition_column].unique().tolist()

    # 1. Kruskal-Wallis test across all conditions
    print("  Running Kruskal-Wallis test...")
    condition_groups = []
    for condition in condition_order:
        if condition in df[condition_column].values:
            condition_data = df[df[condition_column] == condition][value_column].dropna()
            if len(condition_data) > 0:
                condition_groups.append(condition_data.tolist())

    if len(condition_groups) >= 2:
        h_stat, p_value = stats.kruskal(*condition_groups)
        kruskal_results = {
            'test': 'Kruskal-Wallis',
            'h_statistic': h_stat,
            'p_value': p_value,
            'significant': p_value < 0.05,
            'n_conditions': len(condition_groups)
        }
        print(f"    H = {h_stat:.4f}, p = {p_value:.4f}")
        if p_value < 0.05:
            print("    -> Significant differences between condition groups")
        else:
            print("    -> No significant differences between condition groups")
    else:
        kruskal_results = None
        print("    Insufficient groups for Kruskal-Wallis test")

    # 2. Pairwise Mann-Whitney U tests
    print("\n  Running pairwise Mann-Whitney U tests...")
    mannwhitney_results = []
    for i, cond1 in enumerate(condition_order):
        for j, cond2 in enumerate(condition_order):
            if i < j:
                if cond1 in df[condition_column].values and cond2 in df[condition_column].values:
                    data1 = df[df[condition_column] == cond1][value_column].dropna()
                    data2 = df[df[condition_column] == cond2][value_column].dropna()
                    if len(data1) > 0 and len(data2) > 0:
                        u_stat, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
                        mannwhitney_results.append({
                            'condition1': cond1,
                            'condition2': cond2,
                            'U_statistic': u_stat,
                            'p_value': p_value,
                            'n1': len(data1),
                            'n2': len(data2),
                            'mean_diff': data1.mean() - data2.mean(),
                            'median_diff': data1.median() - data2.median()
                        })

    if mannwhitney_results:
        mannwhitney_df = pd.DataFrame(mannwhitney_results)
        _, q_values, _, _ = multipletests(mannwhitney_df['p_value'], method='fdr_bh')
        mannwhitney_df['q_value'] = q_values
        mannwhitney_df.to_csv('pairwise_mannwhitney_results.tsv', sep='\t', index=False)
        print(f"    Saved {len(mannwhitney_df)} pairwise comparisons to 'pairwise_mannwhitney_results.tsv'")
        sig = mannwhitney_df[mannwhitney_df['q_value'] < 0.05]
        if not sig.empty:
            print("\n    Significant pairwise comparisons (q < 0.05):")
            for _, row in sig.iterrows():
                print(f"      {row['condition1']} vs {row['condition2']}: "
                      f"U={row['U_statistic']:.0f}, q={row['q_value']:.4f}, "
                      f"Δmean={row['mean_diff']:.3f}")
    else:
        mannwhitney_df = pd.DataFrame()
        print("    No valid pairwise comparisons")

    # 3. Levene's test for homogeneity of variances
    print("\n  Running Levene's tests for homogeneity of variances...")
    levene_results = []
    for i, cond1 in enumerate(condition_order):
        for j, cond2 in enumerate(condition_order):
            if i < j:
                if cond1 in df[condition_column].values and cond2 in df[condition_column].values:
                    data1 = df[df[condition_column] == cond1][value_column].dropna()
                    data2 = df[df[condition_column] == cond2][value_column].dropna()
                    if len(data1) > 0 and len(data2) > 0:
                        levene_stat, p_value = stats.levene(data1, data2)
                        levene_results.append({
                            'condition1': cond1,
                            'condition2': cond2,
                            'levene_statistic': levene_stat,
                            'p_value': p_value,
                            'n1': len(data1),
                            'n2': len(data2),
                            'variance_ratio': data1.var() / data2.var() if data2.var() > 0 else np.nan
                        })

    if levene_results:
        levene_df = pd.DataFrame(levene_results)
        _, q_values, _, _ = multipletests(levene_df['p_value'], method='fdr_bh')
        levene_df['q_value'] = q_values
        levene_df.to_csv('levene_homogeneity_tests.tsv', sep='\t', index=False)
        print(f"    Saved {len(levene_df)} variance tests to 'levene_homogeneity_tests.tsv'")
    else:
        levene_df = pd.DataFrame()
        print("    No valid variance comparisons")

    return kruskal_results, mannwhitney_df, levene_df


def analyze_environmental_correlations(df, value_column='RPS'):
    """
    Analyze correlations between RPS and environmental variables.

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe
    value_column : str
        Column name for the value to analyze

    Returns
    -------
    pandas.DataFrame
        Correlation results
    """
    print(f"\nAnalyzing correlations between {value_column} and environmental variables...")

    env_variables = [
        'collection date', 'ammonium', 'nitrate', 'nitrite',
        'phosphate', 'silicate', 'chlorophyll', 'density', 'temperature'
    ]
    available_vars = [var for var in env_variables if var in df.columns]
    print(f"  Available environmental variables: {', '.join(available_vars)}")

    correlation_results = []
    for env_var in available_vars:
        valid_data = df[[value_column, env_var]].dropna()
        if len(valid_data) >= 3:
            rho, p_value = stats.spearmanr(valid_data[value_column], valid_data[env_var])
            correlation_results.append({
                'genomic_variable': value_column,
                'environmental_variable': env_var,
                'spearman_rho': rho,
                'p_value': p_value,
                'n_observations': len(valid_data),
                'interpretation': interpret_correlation_strength(abs(rho))
            })

    if correlation_results:
        corr_df = pd.DataFrame(correlation_results)
        _, q_values, _, _ = multipletests(corr_df['p_value'], method='fdr_bh')
        corr_df['q_value'] = q_values
        corr_df['abs_rho'] = corr_df['spearman_rho'].abs()
        corr_df = corr_df.sort_values('abs_rho', ascending=False)
        corr_df.to_csv('environmental_correlations.tsv', sep='\t', index=False)
        print(f"  Saved {len(corr_df)} correlations to 'environmental_correlations.tsv'")
        print("\n  Top correlations (absolute value):")
        for i, row in corr_df.head(5).iterrows():
            sig = "***" if row['q_value'] < 0.001 else "**" if row['q_value'] < 0.01 else "*" if row['q_value'] < 0.05 else ""
            print(f"    {row['environmental_variable']}: ρ = {row['spearman_rho']:.3f}{sig} "
                  f"(q = {row['q_value']:.4f}, n = {row['n_observations']})")
    else:
        corr_df = pd.DataFrame()
        print("  No valid correlations found")

    return corr_df


def interpret_correlation_strength(rho):
    """Interpret the strength of a correlation coefficient."""
    if rho >= 0.9:
        return "very strong"
    elif rho >= 0.7:
        return "strong"
    elif rho >= 0.5:
        return "moderate"
    elif rho >= 0.3:
        return "weak"
    elif rho >= 0.1:
        return "very weak"
    else:
        return "negligible"


def run_random_forest_analysis(df, target_column='RPS'):
    """
    Run Random Forest analysis to identify important environmental predictors.

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe
    target_column : str
        Target variable column

    Returns
    -------
    pandas.Series
        Feature importances
    """
    print(f"\nRunning Random Forest analysis for {target_column}...")

    feature_columns = [
        'collection date', 'ammonium', 'nitrate', 'nitrite',
        'phosphate', 'silicate', 'chlorophyll', 'density', 'temperature'
    ]
    available_features = [col for col in feature_columns if col in df.columns]

    if not available_features:
        print("  No environmental features available for Random Forest")
        return pd.Series()

    print(f"  Using features: {', '.join(available_features)}")

    X = df[available_features].copy()
    y = df[target_column].copy()
    valid_mask = ~(X.isna().any(axis=1) | y.isna())
    X = X[valid_mask]
    y = y[valid_mask]

    if len(X) < 10:
        print(f"  Insufficient data for Random Forest (n={len(X)})")
        return pd.Series()

    regr = RandomForestRegressor(max_depth=2, random_state=0, n_estimators=100, n_jobs=-1)
    regr.fit(X, y)

    importances = regr.feature_importances_
    varimp = pd.Series(importances, index=available_features, name='variable_importance').sort_values(ascending=False)
    varimp.to_csv('random_forest_feature_importances.tsv', sep='\t')
    print(f"  Feature importances saved to 'random_forest_feature_importances.tsv'")
    print("\n  Random Forest feature importances:")
    for feature, importance in varimp.items():
        print(f"    {feature}: {importance:.4f}")
    print(f"  R² score on training data: {regr.score(X, y):.4f}")

    return varimp


def create_visualizations(df, value_column='RPS', condition_column='cond1', condition_order=None):
    """
    Create visualizations of RPS distributions across conditions.
    """
    print(f"\nCreating visualizations for {value_column}...")

    if condition_order is None:
        condition_order = df[condition_column].unique().tolist()

    # 1. KDE plot
    print("  Creating KDE plot...")
    plt.figure(figsize=(10, 6))
    for condition in condition_order:
        if condition in df[condition_column].values:
            data = df[df[condition_column] == condition][value_column].dropna()
            if len(data) > 0:
                sns.kdeplot(data, label=condition, fill=True, alpha=0.3)
    plt.xlabel(value_column)
    plt.ylabel('Density')
    plt.title(f'Distribution of {value_column} Across Experimental Conditions')
    plt.legend(title='Condition')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'kde_{value_column}_by_condition.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'kde_{value_column}_by_condition.svg', bbox_inches='tight')
    plt.close()

    # 2. Swarm + box plot
    print("  Creating swarm/box plot...")
    plt.figure(figsize=(12, 6))
    sampled_data = pd.DataFrame()
    for condition in condition_order:
        if condition in df[condition_column].values:
            condition_df = df[df[condition_column] == condition]
            sample_size = min(200, len(condition_df))
            sampled = condition_df.sample(sample_size, replace=(sample_size > len(condition_df)))
            sampled_data = pd.concat([sampled_data, sampled])
    ax = sns.swarmplot(data=sampled_data, x=condition_column, y=value_column,
                       order=condition_order, size=2, alpha=0.5, color='gray')
    sns.boxplot(data=df, x=condition_column, y=value_column, order=condition_order,
                width=0.4, color='white', showfliers=False, ax=ax)
    plt.xlabel('Experimental Condition')
    plt.ylabel(value_column)
    plt.title(f'{value_column} Distribution Across Experimental Conditions')
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'box_swarm_{value_column}_by_condition.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'box_swarm_{value_column}_by_condition.svg', bbox_inches='tight')
    plt.close()

    # 3. Violin plot
    print("  Creating violin plot...")
    plt.figure(figsize=(10, 6))
    sns.violinplot(data=df, x=condition_column, y=value_column,
                   order=condition_order, inner='quartile', cut=0)
    plt.xlabel('Experimental Condition')
    plt.ylabel(value_column)
    plt.title(f'Violin Plot of {value_column} Across Conditions')
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'violin_{value_column}_by_condition.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'violin_{value_column}_by_condition.svg', bbox_inches='tight')
    plt.close()

    print("  Visualizations saved as PNG and SVG files")


# ------------------------------------------------------------------------------
# New CCF analysis functions (from pipeline.py)
# ------------------------------------------------------------------------------
def safe_ccf(x, y, adjusted=True, fft=True, nlags=None, alpha=None):
    """
    Safe wrapper around ccf that handles edge cases.
    """
    try:
        x = np.asarray(x)
        y = np.asarray(y)
        x = x[~np.isnan(x) & ~np.isinf(x)]
        y = y[~np.isnan(y) & ~np.isinf(y)]

        if len(x) < 2 or len(y) < 2:
            return None, None, None

        if nlags is None:
            nlags = min(len(x), len(y)) - 1
        elif nlags >= min(len(x), len(y)):
            nlags = min(len(x), len(y)) - 1

        ccf_result = ccf(x, y, adjusted=adjusted, fft=fft, nlags=nlags, alpha=alpha)
        lags = np.arange(-nlags, nlags + 1)

        return ccf_result, lags, nlags

    except Exception as e:
        print(f"Error in safe_ccf: {e}")
        return None, None, None


def parse_ccf_result_safe(ccf_values, lags, condition, n_obs):
    """
    Parse cross-correlation results with safe bounds checking.
    """
    if ccf_values is None or lags is None or len(ccf_values) == 0:
        return None

    try:
        min_len = min(len(ccf_values), len(lags))
        ccf_values = ccf_values[:min_len]
        lags = lags[:min_len]

        abs_ccf = np.abs(ccf_values)

        max_idx = np.argmax(abs_ccf) if len(abs_ccf) > 0 else 0
        max_corr = ccf_values[max_idx] if max_idx < len(ccf_values) else 0
        max_lag = lags[max_idx] if max_idx < len(lags) else 0

        positive_idx = np.argmax(ccf_values) if len(ccf_values) > 0 else 0
        positive_corr = ccf_values[positive_idx]
        positive_lag = lags[positive_idx] if positive_idx < len(lags) else 0

        negative_idx = np.argmin(ccf_values) if len(ccf_values) > 0 else 0
        negative_corr = ccf_values[negative_idx]
        negative_lag = lags[negative_idx] if negative_idx < len(lags) else 0

        mean_corr = np.mean(ccf_values) if len(ccf_values) > 0 else 0
        std_corr = np.std(ccf_values) if len(ccf_values) > 1 else 0

        n_effective = len(ccf_values)
        ci_95 = 1.96 / np.sqrt(n_effective) if n_effective > 0 else 0

        significant_pos = np.sum(ccf_values > ci_95) if n_effective > 0 else 0
        significant_neg = np.sum(ccf_values < -ci_95) if n_effective > 0 else 0

        if len(ccf_values) > 1 and len(lags) > 1:
            positive_area = np.trapezoid(np.clip(ccf_values, 0, None), lags)
            negative_area = np.trapezoid(np.clip(ccf_values, None, 0), lags)
            total_area = np.trapezoid(ccf_values, lags)
        else:
            positive_area = negative_area = total_area = 0

        decay_lag = None
        if max_corr > 0 and len(ccf_values) > max_idx:
            half_max = max_corr / 2
            decay_indices = np.where(ccf_values[max_idx:] < half_max)[0]
            if len(decay_indices) > 0:
                decay_idx = decay_indices[0] + max_idx
                decay_lag = lags[decay_idx] if decay_idx < len(lags) else None

        zero_crossings = np.where(np.diff(np.sign(ccf_values)))[0]
        first_zero_crossing = lags[zero_crossings[0]] if len(zero_crossings) > 0 and zero_crossings[0] < len(lags) else None

        return {
            'condition': condition,
            'n_observations': n_obs,
            'n_lags': len(lags) // 2,
            'max_abs_correlation': np.max(abs_ccf) if len(abs_ccf) > 0 else 0,
            'max_correlation': max_corr,
            'max_lag': max_lag,
            'positive_peak': positive_corr,
            'positive_lag': positive_lag,
            'negative_peak': negative_corr,
            'negative_lag': negative_lag,
            'mean_correlation': mean_corr,
            'std_correlation': std_corr,
            'ci_95': ci_95,
            'significant_pos_count': significant_pos,
            'significant_neg_count': significant_neg,
            'positive_area': positive_area,
            'negative_area': negative_area,
            'total_area': total_area,
            'decay_lag': decay_lag,
            'first_zero_crossing': first_zero_crossing
        }

    except Exception as e:
        print(f"Error parsing CCF for condition {condition}: {e}")
        return None


def analyze_ccf_across_conditions_fixed(subdf, order, cond1_column='cond1',
                                       x_col='collection date', y_col='RPS',
                                       adjusted=True, fft=True, nlags=None, alpha=None):
    """
    Fixed CCF analysis across conditions with proper error handling.
    """
    full_results = {
        'summary': [],
        'curves': {},
        'statistics': {},
        'successful_conditions': [],
        'failed_conditions': []
    }

    for c in order:
        subset = subdf[subdf[cond1_column] == c]
        if len(subset) < 3:
            full_results['failed_conditions'].append((c, f"Insufficient data: {len(subset)} points"))
            continue

        try:
            x_data = subset[x_col].dropna().values
            y_data = subset[y_col].dropna().values
            min_len = min(len(x_data), len(y_data))
            if min_len < 3:
                full_results['failed_conditions'].append((c, f"Only {min_len} valid observations"))
                continue

            x_data = x_data[:min_len]
            y_data = y_data[:min_len]

            ccf_result, lags, actual_nlags = safe_ccf(x_data, y_data, adjusted=adjusted,
                                                       fft=fft, nlags=nlags, alpha=alpha)
            if ccf_result is None:
                raise ValueError("CCF returned None")

            parsed = parse_ccf_result_safe(ccf_result, lags, c, min_len)
            if parsed is None:
                raise ValueError("Failed to parse CCF results")

            full_results['summary'].append(parsed)
            full_results['curves'][c] = {
                'lags': lags,
                'ccf': ccf_result,
                'ci_95': parsed['ci_95'],
                'n_obs': min_len,
                'x_data': x_data,
                'y_data': y_data
            }
            full_results['successful_conditions'].append(c)

        except Exception as e:
            full_results['failed_conditions'].append((c, str(e)))

    if full_results['summary']:
        summary_df = pd.DataFrame(full_results['summary'])
        summary_df['dominance'] = np.where(summary_df['max_correlation'] > 0, 'positive', 'negative')
        bins = [0, 0.1, 0.3, 0.6, 0.9, 1.0]
        labels = ['very weak', 'weak', 'moderate', 'strong', 'very strong']
        summary_df['strength_category'] = pd.cut(np.abs(summary_df['max_correlation']),
                                                 bins=bins, labels=labels, include_lowest=True)
        full_results['summary_df'] = summary_df
        full_results['statistics'] = run_statistical_tests(summary_df, full_results)

    return full_results


def run_statistical_tests(summary_df, full_results):
    """
    Run statistical tests on CCF results.
    """
    stats_dict = {}

    if len(summary_df) < 2:
        stats_dict['warning'] = "Insufficient conditions for statistical tests"
        return stats_dict

    try:
        # 1. Test if mean correlation differs from zero
        if len(summary_df) > 1:
            t_stat, p_value = stats.ttest_1samp(summary_df['max_correlation'], 0)
            stats_dict['mean_correlation_test'] = {
                't_statistic': t_stat,
                'p_value': p_value,
                'significant': p_value < 0.05
            }

        # 2. Test normality of correlations
        if len(summary_df) > 3:
            shapiro_stat, shapiro_p = stats.shapiro(summary_df['max_correlation'])
            stats_dict['normality_test'] = {
                'shapiro_statistic': shapiro_stat,
                'p_value': shapiro_p,
                'normal': shapiro_p > 0.05
            }

        # 3. Compare positive vs negative dominance groups
        pos_group = summary_df[summary_df['dominance'] == 'positive']['max_correlation']
        neg_group = summary_df[summary_df['dominance'] == 'negative']['max_correlation']
        if len(pos_group) > 1 and len(neg_group) > 1:
            levene_stat, levene_p = stats.levene(pos_group, neg_group)
            equal_var = levene_p > 0.05
            t_stat, p_value = stats.ttest_ind(pos_group, neg_group, equal_var=equal_var)
            stats_dict['dominance_comparison'] = {
                't_statistic': t_stat,
                'p_value': p_value,
                'significant': p_value < 0.05,
                'equal_variances': equal_var,
                'n_positive': len(pos_group),
                'n_negative': len(neg_group),
                'mean_positive': pos_group.mean(),
                'mean_negative': neg_group.mean()
            }

        # 4. Correlation between lag and correlation strength
        if len(summary_df) > 2:
            corr_coef, corr_p = stats.pearsonr(np.abs(summary_df['max_lag']), np.abs(summary_df['max_correlation']))
            stats_dict['lag_strength_correlation'] = {
                'correlation_coefficient': corr_coef,
                'p_value': corr_p,
                'significant': corr_p < 0.05
            }

        # 5. Kruskal-Wallis test for strength categories
        if len(summary_df['strength_category'].unique()) > 1:
            groups = []
            for cat in summary_df['strength_category'].unique():
                groups.append(summary_df[summary_df['strength_category'] == cat]['max_correlation'].values)
            if len(groups) > 1:
                h_stat, p_value = stats.kruskal(*groups)
                stats_dict['strength_category_test'] = {
                    'h_statistic': h_stat,
                    'p_value': p_value,
                    'significant': p_value < 0.05
                }

        # 6. Bootstrap confidence intervals for mean correlation
        if len(summary_df) > 1:
            n_bootstraps = 1000
            bootstrap_means = []
            for _ in range(n_bootstraps):
                sample = np.random.choice(summary_df['max_correlation'], size=len(summary_df), replace=True)
                bootstrap_means.append(np.mean(sample))
            ci_lower = np.percentile(bootstrap_means, 2.5)
            ci_upper = np.percentile(bootstrap_means, 97.5)
            stats_dict['bootstrap_ci'] = {
                'mean': np.mean(summary_df['max_correlation']),
                'ci_95_lower': ci_lower,
                'ci_95_upper': ci_upper,
                'n_bootstraps': n_bootstraps
            }

        # 7. Effect sizes
        if 'dominance_comparison' in stats_dict:
            pooled_std = np.sqrt(((len(pos_group)-1)*pos_group.std()**2 +
                                 (len(neg_group)-1)*neg_group.std()**2) /
                                 (len(pos_group) + len(neg_group) - 2))
            cohens_d = (pos_group.mean() - neg_group.mean()) / pooled_std
            stats_dict['effect_sizes'] = {
                'cohens_d_dominance': cohens_d,
                'interpretation': 'small' if abs(cohens_d) < 0.5 else
                                 'medium' if abs(cohens_d) < 0.8 else 'large'
            }

    except Exception as e:
        stats_dict['error'] = f"Statistical tests failed: {str(e)}"

    return stats_dict


def print_analysis_report(full_results):
    """
    Print comprehensive CCF analysis report.
    """
    print("="*80)
    print("CCF ANALYSIS REPORT")
    print("="*80)

    print(f"\n1. DATA SUMMARY:")
    print(f"   Successful conditions: {len(full_results['successful_conditions'])}")
    if full_results['successful_conditions']:
        print(f"   Conditions analyzed: {', '.join(full_results['successful_conditions'])}")
    if full_results['failed_conditions']:
        print(f"   Failed conditions: {len(full_results['failed_conditions'])}")
        for cond, reason in full_results['failed_conditions'][:5]:
            print(f"     - {cond}: {reason}")
        if len(full_results['failed_conditions']) > 5:
            print(f"     ... and {len(full_results['failed_conditions']) - 5} more")

    if 'summary_df' in full_results and not full_results['summary_df'].empty:
        print(f"\n2. CROSS-CORRELATION SUMMARY:")
        print(f"   Total conditions analyzed: {len(full_results['summary_df'])}")
        print(f"\n   Descriptive statistics:")
        print(full_results['summary_df'][['max_correlation', 'max_lag', 'n_observations']].describe().round(3))

        if 'dominance' in full_results['summary_df'].columns:
            dominance_counts = full_results['summary_df']['dominance'].value_counts()
            print(f"\n   Dominance distribution:")
            for dom, cnt in dominance_counts.items():
                print(f"     {dom}: {cnt} conditions")

        if 'strength_category' in full_results['summary_df'].columns:
            strength_counts = full_results['summary_df']['strength_category'].value_counts()
            print(f"\n   Strength categories:")
            for strength, cnt in strength_counts.items():
                print(f"     {strength}: {cnt} conditions")

    if 'statistics' in full_results and full_results['statistics']:
        print(f"\n3. STATISTICAL TESTS:")
        for test_name, test_results in full_results['statistics'].items():
            if test_name in ['warning', 'error']:
                print(f"   {test_name.upper()}: {test_results}")
                continue
            print(f"\n   {test_name.replace('_', ' ').title()}:")
            if isinstance(test_results, dict):
                for key, value in test_results.items():
                    if isinstance(value, float):
                        print(f"     {key}: {value:.4f}")
                    elif isinstance(value, bool):
                        print(f"     {key}: {'Yes' if value else 'No'}")
                    else:
                        print(f"     {key}: {value}")
            else:
                print(f"     {test_results}")

    if 'summary_df' in full_results and not full_results['summary_df'].empty:
        print(f"\n4. TOP CORRELATIONS:")
        print(f"\n   Strongest positive correlations:")
        pos_sorted = full_results['summary_df'].sort_values('positive_peak', ascending=False)
        for _, row in pos_sorted.head(3).iterrows():
            print(f"     {row['condition']}: {row['positive_peak']:.3f} at lag {row['positive_lag']}")

        print(f"\n   Strongest negative correlations:")
        neg_sorted = full_results['summary_df'].sort_values('negative_peak', ascending=True)
        for _, row in neg_sorted.head(3).iterrows():
            print(f"     {row['condition']}: {row['negative_peak']:.3f} at lag {row['negative_lag']}")

        print(f"\n   Strongest overall correlations:")
        abs_sorted = full_results['summary_df'].sort_values('max_abs_correlation', ascending=False)
        for _, row in abs_sorted.head(3).iterrows():
            print(f"     {row['condition']}: {row['max_correlation']:.3f} at lag {row['max_lag']} "
                  f"(abs: {row['max_abs_correlation']:.3f})")

    print("\n" + "="*80)
    print("END OF REPORT")
    print("="*80)


# ------------------------------------------------------------------------------
# Main execution
# ------------------------------------------------------------------------------
def main():
    """Main execution function."""
    # Configuration (update these paths as needed)
    CARBONCOMP_PATH = 'carboncomp_output.tsv.xz'
    METADATA_PATH = 'data/refactored_metadata.tsv'

    print("=" * 60)
    print("ANALYSIS OF RPS IN PERTURBATION EXPERIMENT (Coello-Camba et al.)")
    print("=" * 60)

    # Load and prepare data
    df = load_and_prepare_data(CARBONCOMP_PATH, METADATA_PATH)

    # Ensure RPS column exists (rename EIT if present)
    if 'EIT' in df.columns and 'RPS' not in df.columns:
        print("\nRenaming 'EIT' column to 'RPS'...")
        df['RPS'] = df['EIT']

    if 'RPS' not in df.columns:
        print("ERROR: No RPS column found in data")
        return

    print(f"\nDataset overview:")
    print(f"  Total observations: {len(df)}")
    print(f"  Unique conditions: {df['cond1'].nunique()}")
    print(f"  RPS range: [{df['RPS'].min():.3f}, {df['RPS'].max():.3f}]")
    print(f"  RPS mean ± std: {df['RPS'].mean():.3f} ± {df['RPS'].std():.3f}")

    # 1. Summary statistics
    stats_df, condition_order = calculate_summary_statistics(df, value_column='RPS')

    # 2. Statistical tests
    kruskal_results, mannwhitney_df, levene_df = perform_statistical_tests(
        df, value_column='RPS', condition_order=condition_order
    )

    # 3. Environmental correlations
    corr_df = analyze_environmental_correlations(df, value_column='RPS')

    # 4. Random Forest analysis
    varimp = run_random_forest_analysis(df, target_column='RPS')

    # 5. Visualizations
    create_visualizations(df, value_column='RPS', condition_order=condition_order)

    # 6. Cross-correlation (CCF) analysis
    print("\n" + "=" * 60)
    print("CROSS-CORRELATION (CCF) ANALYSIS")
    print("=" * 60)
    # Sort by condition and collection date as required for time-series
    df_sorted = df.sort_values(by=['cond1', 'collection date']).copy()
    full_results = analyze_ccf_across_conditions_fixed(
        subdf=df_sorted,
        order=condition_order,
        cond1_column='cond1',
        x_col='collection date',
        y_col='RPS',
        adjusted=True,
        fft=True,
        nlags=None,
        alpha=None
    )
    print_analysis_report(full_results)

    # Save CCF summary if available
    if 'summary_df' in full_results:
        full_results['summary_df'].to_csv('ccf_analysis_results.csv', index=False)
        print("\nCCF summary saved to: ccf_analysis_results.csv")

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()

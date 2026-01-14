"""
Script for analyzing CaCo-derived metrics with perturbation experiment data.
Analyzes RPS (Resource Partitioning Score) across different nutrient treatments
and environmental conditions (Coello-Camba et al. dataset).
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
from sklearn.ensemble import RandomForestRegressor
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

def load_and_prepare_data(carboncomp_path, metadata_path):
    """
    Load and merge genomic data with perturbation experiment metadata.
    
    Parameters:
    -----------
    carboncomp_path : str
        Path to carboncomp_output.tsv.xz file
    metadata_path : str
        Path to metadata file
    
    Returns:
    --------
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

def calculate_summary_statistics(df, value_column='RPS', condition_column='cond1'):
    """
    Calculate summary statistics for RPS across conditions.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe
    value_column : str
        Column name for the value to analyze
    condition_column : str
        Column name containing condition labels
    
    Returns:
    --------
    pandas.DataFrame
        Summary statistics per condition
    """
    print(f"\nCalculating summary statistics for {value_column}...")
    
    # Define condition order for consistent plotting
    condition_order = [
        'N,P Day 1 single addition',
        'N,P 2 week continuous addition',
        'control mesocosm',
        'N,P,Si Day 1 single addition',
        'N,P,Si 2 week continuous addition'
    ]
    
    # Calculate statistics per condition
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
    
    # Save to file
    stats_df.to_csv('RPS_summary_statistics_per_condition.tsv', sep='\t', index=False)
    print(f"  Summary statistics saved to 'RPS_summary_statistics_per_condition.tsv'")
    
    return stats_df, condition_order

def perform_statistical_tests(df, value_column='RPS', condition_column='cond1', condition_order=None):
    """
    Perform statistical tests between conditions.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe
    value_column : str
        Column name for the value to analyze
    condition_column : str
        Column name containing condition labels
    condition_order : list
        Order of conditions for comparisons
    
    Returns:
    --------
    tuple: (kruskal_results, mannwhitney_results, levene_results)
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
            if i < j:  # Avoid duplicates and self-comparisons
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
        
        # Apply FDR correction
        _, q_values, _, _ = multipletests(mannwhitney_df['p_value'], method='fdr_bh')
        mannwhitney_df['q_value'] = q_values
        
        # Save to file
        mannwhitney_df.to_csv('pairwise_mannwhitney_results.tsv', sep='\t', index=False)
        print(f"    Saved {len(mannwhitney_df)} pairwise comparisons to 'pairwise_mannwhitney_results.tsv'")
        
        # Print significant comparisons
        sig_comparisons = mannwhitney_df[mannwhitney_df['q_value'] < 0.05]
        if not sig_comparisons.empty:
            print(f"\n    Significant pairwise comparisons (q < 0.05):")
            for _, row in sig_comparisons.iterrows():
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
        
        # Apply FDR correction
        _, q_values, _, _ = multipletests(levene_df['p_value'], method='fdr_bh')
        levene_df['q_value'] = q_values
        
        # Save to file
        levene_df.to_csv('levene_homogeneity_tests.tsv', sep='\t', index=False)
        print(f"    Saved {len(levene_df)} variance tests to 'levene_homogeneity_tests.tsv'")
    else:
        levene_df = pd.DataFrame()
        print("    No valid variance comparisons")
    
    return kruskal_results, mannwhitney_df, levene_df

def analyze_environmental_correlations(df, value_column='RPS'):
    """
    Analyze correlations between RPS and environmental variables.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe
    value_column : str
        Column name for the value to analyze
    
    Returns:
    --------
    pandas.DataFrame
        Correlation results
    """
    print(f"\nAnalyzing correlations between {value_column} and environmental variables...")
    
    # Environmental variables to analyze
    env_variables = [
        'collection date', 'ammonium', 'nitrate', 'nitrite',
        'phosphate', 'silicate', 'chlorophyll', 'density', 'temperature'
    ]
    
    # Check which variables are available
    available_vars = [var for var in env_variables if var in df.columns]
    print(f"  Available environmental variables: {', '.join(available_vars)}")
    
    correlation_results = []
    
    for env_var in available_vars:
        # Remove rows with NaN in either variable
        valid_data = df[[value_column, env_var]].dropna()
        
        if len(valid_data) >= 3:  # Minimum for correlation
            # Calculate Spearman correlation
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
        
        # Apply FDR correction
        _, q_values, _, _ = multipletests(corr_df['p_value'], method='fdr_bh')
        corr_df['q_value'] = q_values
        
        # Sort by absolute correlation strength
        corr_df['abs_rho'] = corr_df['spearman_rho'].abs()
        corr_df = corr_df.sort_values('abs_rho', ascending=False)
        
        # Save to file
        corr_df.to_csv('environmental_correlations.tsv', sep='\t', index=False)
        print(f"  Saved {len(corr_df)} correlations to 'environmental_correlations.tsv'")
        
        # Print top correlations
        print(f"\n  Top correlations (absolute value):")
        for i, row in corr_df.head(5).iterrows():
            significance = "***" if row['q_value'] < 0.001 else "**" if row['q_value'] < 0.01 else "*" if row['q_value'] < 0.05 else ""
            print(f"    {row['environmental_variable']}: ρ = {row['spearman_rho']:.3f}{significance} "
                  f"(q = {row['q_value']:.4f}, n = {row['n_observations']})")
    else:
        corr_df = pd.DataFrame()
        print("  No valid correlations found")
    
    return corr_df

def interpret_correlation_strength(rho):
    """
    Interpret the strength of a correlation coefficient.
    
    Parameters:
    -----------
    rho : float
        Correlation coefficient (absolute value)
    
    Returns:
    --------
    str
        Interpretation string
    """
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
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe
    target_column : str
        Target variable column
    
    Returns:
    --------
    pandas.Series
        Feature importances
    """
    print(f"\nRunning Random Forest analysis for {target_column}...")
    
    # Define feature columns
    feature_columns = [
        'collection date', 'ammonium', 'nitrate', 'nitrite',
        'phosphate', 'silicate', 'chlorophyll', 'density', 'temperature'
    ]
    
    # Check which features are available
    available_features = [col for col in feature_columns if col in df.columns]
    
    if not available_features:
        print("  No environmental features available for Random Forest")
        return pd.Series()
    
    print(f"  Using features: {', '.join(available_features)}")
    
    # Prepare data
    X = df[available_features].copy()
    y = df[target_column].copy()
    
    # Remove rows with NaN in features or target
    valid_mask = ~(X.isna().any(axis=1) | y.isna())
    X = X[valid_mask]
    y = y[valid_mask]
    
    if len(X) < 10:
        print(f"  Insufficient data for Random Forest (n={len(X)})")
        return pd.Series()
    
    # Train Random Forest
    regr = RandomForestRegressor(
        max_depth=2,
        random_state=42,
        n_estimators=100,
        n_jobs=-1
    )
    
    regr.fit(X, y)
    
    # Get feature importances
    importances = regr.feature_importances_
    varimp = pd.Series(
        importances,
        index=available_features,
        name='variable_importance'
    ).sort_values(ascending=False)
    
    # Save to file
    varimp.to_csv('random_forest_feature_importances.tsv', sep='\t')
    print(f"  Feature importances saved to 'random_forest_feature_importances.tsv'")
    
    # Print results
    print(f"\n  Random Forest feature importances:")
    for feature, importance in varimp.items():
        print(f"    {feature}: {importance:.4f}")
    
    # Calculate R² score
    r2_score = regr.score(X, y)
    print(f"  R² score on training data: {r2_score:.4f}")
    
    return varimp

def create_visualizations(df, value_column='RPS', condition_column='cond1', condition_order=None):
    """
    Create visualizations of RPS distributions across conditions.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe
    value_column : str
        Column name for the value to visualize
    condition_column : str
        Column name containing condition labels
    condition_order : list
        Order of conditions for plotting
    """
    print(f"\nCreating visualizations for {value_column}...")
    
    if condition_order is None:
        condition_order = df[condition_column].unique().tolist()
    
    # 1. KDE plot for each condition
    print("  Creating KDE plot...")
    plt.figure(figsize=(10, 6))
    
    for condition in condition_order:
        if condition in df[condition_column].values:
            condition_data = df[df[condition_column] == condition][value_column].dropna()
            if len(condition_data) > 0:
                sns.kdeplot(condition_data, label=condition, fill=True, alpha=0.3)
    
    plt.xlabel(value_column)
    plt.ylabel('Density')
    plt.title(f'Distribution of {value_column} Across Experimental Conditions')
    plt.legend(title='Condition')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'kde_{value_column}_by_condition.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'kde_{value_column}_by_condition.svg', bbox_inches='tight')
    plt.close()
    
    # 2. Swarm plot with box plot overlay
    print("  Creating swarm/box plot...")
    plt.figure(figsize=(12, 6))
    
    # Sample data for swarm plot (to avoid overplotting)
    sampled_data = pd.DataFrame()
    for condition in condition_order:
        if condition in df[condition_column].values:
            condition_df = df[df[condition_column] == condition]
            # Sample up to 200 points per condition
            sample_size = min(200, len(condition_df))
            sampled = condition_df.sample(sample_size, replace=(sample_size > len(condition_df)))
            sampled_data = pd.concat([sampled_data, sampled])
    
    # Create swarm plot
    ax = sns.swarmplot(
        data=sampled_data,
        x=condition_column,
        y=value_column,
        order=condition_order,
        size=2,
        alpha=0.5,
        color='gray'
    )
    
    # Overlay box plot
    sns.boxplot(
        data=df,
        x=condition_column,
        y=value_column,
        order=condition_order,
        width=0.4,
        color='white',
        showfliers=False,
        ax=ax
    )
    
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
    
    sns.violinplot(
        data=df,
        x=condition_column,
        y=value_column,
        order=condition_order,
        inner='quartile',
        cut=0
    )
    
    plt.xlabel('Experimental Condition')
    plt.ylabel(value_column)
    plt.title(f'Violin Plot of {value_column} Across Conditions')
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'violin_{value_column}_by_condition.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'violin_{value_column}_by_condition.svg', bbox_inches='tight')
    plt.close()
    
    print(f"  Visualizations saved as PNG and SVG files")

def main():
    """Main execution function."""
    # Configuration
    CARBONCOMP_PATH = 'carboncomp_output.tsv.xz'  # CaCo output from running genomes from Coello-Camba et al.
    METADATA_PATH = 'data/refactored_metadata.tsv'  # Filtered MAGs used from COEL20-1 project in mOTUs database referring to HQ MAGs from Coello-Camba
    
    print("=" * 60)
    print("ANALYSIS OF RPS IN PERTURBATION EXPERIMENT (Coello-Camba et al.)")
    print("=" * 60)
    
    # Load and prepare data
    df = load_and_prepare_data(CARBONCOMP_PATH, METADATA_PATH)
    
    # Ensure RPS column exists (rename if EIT exists)
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
    
    # 1. Calculate summary statistics
    stats_df, condition_order = calculate_summary_statistics(df, value_column='RPS')
    
    # 2. Perform statistical tests
    kruskal_results, mannwhitney_df, levene_df = perform_statistical_tests(
        df, value_column='RPS', condition_order=condition_order
    )
    
    # 3. Analyze environmental correlations
    corr_df = analyze_environmental_correlations(df, value_column='RPS')
    
    # 4. Run Random Forest analysis
    varimp = run_random_forest_analysis(df, target_column='RPS')
    
    # 5. Create visualizations
    create_visualizations(df, value_column='RPS', condition_order=condition_order)
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    
    # Print summary of key findings
    print("\nKEY FINDINGS:")
    print("-" * 40)
    
    if kruskal_results and kruskal_results['significant']:
        print(f"• Significant differences in RPS across conditions (H = {kruskal_results['h_statistic']:.1f}, p = {kruskal_results['p_value']:.4f})")
    
    if not mannwhitney_df.empty:
        sig_pairs = mannwhitney_df[mannwhitney_df['q_value'] < 0.05]
        if len(sig_pairs) > 0:
            print(f"• {len(sig_pairs)} significant pairwise differences between conditions")
    
    if not corr_df.empty:
        top_corr = corr_df.iloc[0]
        if top_corr['q_value'] < 0.05:
            print(f"• Strongest correlation: RPS with {top_corr['environmental_variable']} "
                  f"(ρ = {top_corr['spearman_rho']:.3f}, q = {top_corr['q_value']:.4f})")
    
    if not varimp.empty:
        top_feature = varimp.index[0]
        print(f"• Most important environmental predictor: {top_feature} "
              f"(importance = {varimp.iloc[0]:.4f})")
    
    print(f"\nAll output files have been saved to current directory.")

if __name__ == "__main__":
    main()
"""
Script for analyzing CaCo-derived genomic metrics with experimental interaction data.
Compares RPS (Resource Partitioning Score) with experimentally measured interaction types.
Based on data from Ono et al. study.
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from sklearn.ensemble import RandomForestRegressor
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
from itertools import combinations

warnings.filterwarnings('ignore')

def load_and_prepare_data(carboncomp_path, excel_path):
    """
    Load and merge genomic and experimental data from Ono et al.
    
    Parameters:
    -----------
    carboncomp_path : str
        Path to carboncomp_output.tsv.xz file
    excel_path : str
        Path to Excel file with supplementary data
    
    Returns:
    --------
    pandas.DataFrame
        Merged dataframe with genomic and experimental metrics
    """
    print("Loading and preparing data...")
    
    # Load metadata for genome name mapping (Table S11)
    print("  Loading Table S11 for genome mapping...")
    metadata_s11 = pd.read_excel(excel_path, sheet_name='Table S11')
    metadata_s11 = metadata_s11.set_index(metadata_s11.columns[2])[metadata_s11.columns[1]].to_dict()
    
    # Load genomic data
    print("  Loading genomic data...")
    df = pd.read_table(carboncomp_path)
    
    # Clean genome identifiers and map to strain names
    df['genome1_original'] = df['genome1']
    df['genome2_original'] = df['genome2']
    
    df['genome1'] = df['genome1'].apply(lambda x: '_'.join(x.split('_')[:2]))
    df['genome2'] = df['genome2'].apply(lambda x: '_'.join(x.split('_')[:2]))
    
    df['genome1'] = df['genome1'].apply(lambda x: metadata_s11.get(x, x))
    df['genome2'] = df['genome2'].apply(lambda x: metadata_s11.get(x, x))
    
    # Load interaction metadata (Table S6)
    print("  Loading Table S6 for interaction data...")
    metadata_s6 = pd.read_excel(excel_path, sheet_name='Table S6')
    metadata_s6 = metadata_s6.set_index(metadata_s6.columns[2])[metadata_s6.columns[1]].to_dict()
    
    # Create dataframe from Table S6 metadata
    metadata_df = pd.DataFrame(list(metadata_s6.items()), columns=['pair', 'interaction_class'])
    metadata_df['genome1'] = metadata_df['pair'].apply(lambda x: x.split(':')[0])
    metadata_df['genome2'] = metadata_df['pair'].apply(lambda x: x.split(':')[1])
    
    # Create symmetric pairs (both directions)
    metadata_df_rev = metadata_df.copy()
    metadata_df_rev = metadata_df_rev.rename(columns={'genome1': 'genome2', 'genome2': 'genome1'})
    metadata_df = pd.concat([metadata_df, metadata_df_rev])
    
    # Merge with genomic data
    print("  Merging data...")
    df = df.merge(metadata_df[['genome1', 'genome2', 'interaction_class']], 
                  on=['genome1', 'genome2'], 
                  how='inner')
    
    print(f"  Final dataset shape: {df.shape}")
    print(f"  Unique strain pairs: {df[['genome1', 'genome2']].drop_duplicates().shape[0]}")
    
    return df

def calculate_correlations(df, genomic_metric='RPS', experimental_metrics=None):
    """
    Calculate Spearman correlations between genomic and experimental metrics.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Dataframe containing both genomic and experimental metrics
    genomic_metric : str
        Name of genomic metric column (RPS)
    experimental_metrics : list
        List of experimental metric column names
    
    Returns:
    --------
    pandas.DataFrame
        Dataframe with correlation results
    """
    if experimental_metrics is None:
        # Get all interaction class columns (one-hot encoded)
        interaction_cols = [col for col in df.columns if col.startswith('interaction_')]
        if interaction_cols:
            experimental_metrics = interaction_cols
        else:
            # If not one-hot encoded, we'll create dummy variables
            experimental_metrics = []
    
    results = []
    
    print("\nCalculating correlations...")
    
    # First, check if we need to create one-hot encoding
    if 'interaction_class' in df.columns and not experimental_metrics:
        # Create one-hot encoding for interaction classes
        interaction_dummies = pd.get_dummies(df['interaction_class'], prefix='interaction')
        df = pd.concat([df, interaction_dummies], axis=1)
        experimental_metrics = interaction_dummies.columns.tolist()
    
    for exp_metric in experimental_metrics:
        if exp_metric in df.columns and genomic_metric in df.columns:
            # Remove NaN values for correlation
            valid_data = df[[genomic_metric, exp_metric]].dropna()
            
            if len(valid_data) > 3:  # Minimum data points for correlation
                r, p = spearmanr(valid_data[genomic_metric], valid_data[exp_metric])
                results.append({
                    'genomic_metric': genomic_metric,
                    'experimental_metric': exp_metric,
                    'rho': r,
                    'p': p,
                    'n': len(valid_data)
                })
    
    if results:
        results_df = pd.DataFrame(results)
        
        # Apply FDR correction
        _, q_values, _, _ = multipletests(results_df['p'], method='fdr_bh')
        results_df['q'] = q_values
        
        # Sort by absolute correlation strength
        results_df['abs_rho'] = results_df['rho'].abs()
        results_df = results_df.sort_values('abs_rho', ascending=False)
        
        return results_df, df
    else:
        print("No valid correlations found.")
        return pd.DataFrame(), df

def prepare_pairwise_data(df, value_column='RPS'):
    """
    Create summary statistics for each unique strain pair.
    Correctly handles 8 strains in 28 unique pairs.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe
    value_column : str
        Column name for the value to analyze (RPS)
    
    Returns:
    --------
    pandas.DataFrame
        Pairwise summary statistics
    """
    print(f"Preparing pairwise data using {value_column}...")
    
    # Create a unique pair ID (order doesn't matter)
    df['pair_id'] = df.apply(
        lambda row: '_'.join(sorted([str(row['genome1']), str(row['genome2'])])), 
        axis=1
    )
    
    # Calculate summary statistics per pair
    pair_summary = []
    
    for pair_id, group in df.groupby('pair_id'):
        # Calculate median RPS for this pair
        median_value = group[value_column].median()
        
        # Get the most frequent interaction class (handle ties)
        interaction_counts = group['interaction_class'].value_counts()
        if not interaction_counts.empty:
            # Get the mode(s)
            max_count = interaction_counts.max()
            modes = interaction_counts[interaction_counts == max_count].index.tolist()
            modal_interaction = modes[0]  # Take first if multiple
            
            # Calculate proportions for each interaction type
            interaction_types = ['Competition', 'Exploitation', 'Neutrality', 'Mutualism']
            proportions = {}
            for itype in interaction_types:
                proportions[f'{itype.lower()}_proportion'] = (group['interaction_class'] == itype).mean()
            
            pair_summary.append({
                'pair_id': pair_id,
                value_column: median_value,
                'modal_interaction': modal_interaction,
                'n_observations': len(group),
                **proportions
            })
    
    pair_df = pd.DataFrame(pair_summary)
    print(f"  Found {len(pair_df)} unique strain pairs")
    
    return pair_df

def permutation_anova_by_interaction_class(df, value_column='RPS', n_permutations=9999):
    """
    Permutation-based ANOVA for interaction classes using unique pairs.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe
    value_column : str
        Column name for the value to analyze
    n_permutations : int
        Number of permutations
    
    Returns:
    --------
    tuple: (observed_f, p_value, interaction_classes)
    """
    print(f"\nRunning permutation ANOVA for {value_column}...")
    
    # Prepare pair summary
    pair_summary = prepare_pairwise_data(df, value_column)
    
    print(f"Number of unique strain pairs: {len(pair_summary)}")
    
    # Create groups by modal interaction
    interaction_classes = sorted(pair_summary['modal_interaction'].unique())
    groups = []
    
    print("\nGroup sizes by interaction class:")
    for ic in interaction_classes:
        group_data = pair_summary[pair_summary['modal_interaction'] == ic][value_column].values
        if len(group_data) > 0:
            groups.append(group_data)
            print(f"  {ic}: {len(group_data)} pairs")
    
    if len(groups) < 2:
        print("Need at least 2 groups with data for ANOVA")
        return None, None, interaction_classes
    
    def anova_f_statistic(groups):
        """Calculate F-statistic for one-way ANOVA."""
        # Remove empty groups
        groups = [g for g in groups if len(g) > 0]
        
        if len(groups) < 2:
            return 0.0
        
        # Flatten all data
        all_data = np.concatenate(groups)
        grand_mean = np.mean(all_data)
        
        # Calculate between-group variability (SSB)
        ssb = 0
        for group in groups:
            n = len(group)
            if n > 0:
                group_mean = np.mean(group)
                ssb += n * (group_mean - grand_mean) ** 2
        
        # Calculate within-group variability (SSW)
        ssw = 0
        for group in groups:
            if len(group) > 0:
                group_mean = np.mean(group)
                ssw += np.sum((group - group_mean) ** 2)
        
        # Degrees of freedom
        k = len(groups)  # number of groups
        n_total = len(all_data)  # total observations
        
        df_between = k - 1
        df_within = n_total - k
        
        # Avoid division by zero
        if df_within == 0 or ssw == 0:
            return 0.0
        
        # Mean squares
        ms_between = ssb / df_between
        ms_within = ssw / df_within
        
        # F-statistic
        return ms_between / ms_within
    
    # Observed F-statistic
    f_observed = anova_f_statistic(groups)
    print(f"\nObserved F-statistic: {f_observed:.4f}")
    
    # Permutation test
    f_permuted = []
    
    # Get all values and group sizes
    all_values = np.concatenate(groups)
    group_sizes = [len(g) for g in groups]
    
    print(f"\nRunning {n_permutations} permutations...")
    for i in range(n_permutations):
        if i % 1000 == 0 and i > 0:
            print(f"  Completed {i} permutations")
        
        # Permute values across all pairs
        permuted_values = np.random.permutation(all_values)
        
        # Split into groups maintaining original sizes
        permuted_groups = []
        start_idx = 0
        for size in group_sizes:
            permuted_groups.append(permuted_values[start_idx:start_idx+size])
            start_idx += size
        
        # Calculate F-statistic for permuted data
        f_perm = anova_f_statistic(permuted_groups)
        f_permuted.append(f_perm)
    
    # Calculate p-value
    f_permuted = np.array(f_permuted)
    p_value = (np.sum(f_permuted >= f_observed) + 1) / (n_permutations + 1)
    
    return f_observed, p_value, interaction_classes

def random_forest_analysis(df, target_column='RPS', feature_columns=None):
    """
    Perform Random Forest analysis to identify important features.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe
    target_column : str
        Target variable column
    feature_columns : list
        Feature columns to use
    
    Returns:
    --------
    pandas.Series
        Feature importances
    """
    print("\nRunning Random Forest analysis...")
    
    if feature_columns is None:
        # Use interaction class dummy variables as features
        if 'interaction_class' in df.columns:
            # Create one-hot encoding
            interaction_dummies = pd.get_dummies(df['interaction_class'], prefix='interaction')
            df_rf = pd.concat([df[[target_column]], interaction_dummies], axis=1)
            feature_columns = interaction_dummies.columns.tolist()
        else:
            print("No feature columns specified and cannot create from interaction_class")
            return pd.Series()
    
    # Prepare data
    X = df_rf[feature_columns]
    y = df_rf[target_column]
    
    # Remove rows with NaN
    valid_mask = ~(X.isna().any(axis=1) | y.isna())
    X = X[valid_mask]
    y = y[valid_mask]
    
    if len(X) < 10:
        print(f"Not enough data for Random Forest (n={len(X)})")
        return pd.Series()
    
    # Train Random Forest
    regr = RandomForestRegressor(max_depth=2, random_state=42, n_estimators=100)
    regr.fit(X, y)
    
    # Get feature importances
    varimp = pd.Series(
        regr.feature_importances_,
        index=feature_columns,
        name='variable_importance'
    ).sort_values(ascending=False)
    
    print("\nRandom Forest feature importances:")
    for feature, importance in varimp.items():
        print(f"  {feature}: {importance:.4f}")
    
    return varimp

def create_interaction_heatmaps(df, value_column='RPS', excel_path=None):
    """
    Create heatmaps showing RPS by carbon source count and interaction class.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe
    value_column : str
        Column name for the value to plot
    excel_path : str
        Path to Excel file for additional metadata
    
    Returns:
    --------
    tuple: (mean_df, median_df)
        DataFrames with summary statistics
    """
    print("\nCreating interaction heatmaps...")
    
    if excel_path is not None and 'carbonsource_count' not in df.columns:
        # Load carbonsource_count from Table S8
        try:
            print("  Loading Table S8 for carbon source count...")
            metadata_s8 = pd.read_excel(excel_path, sheet_name='Table S8')
            metadata_s8 = metadata_s8[['coculture_species_GFP', 'coculture_species_mScarlet', 
                                      'carbonsource_count', 'interaction_class']]
            
            # Create symmetric pairs
            a = metadata_s8.rename(
                {'coculture_species_GFP': 'genome1', 
                 'coculture_species_mScarlet': 'genome2'}, 
                axis=1
            )
            b = metadata_s8.rename(
                {'coculture_species_GFP': 'genome2', 
                 'coculture_species_mScarlet': 'genome1'}, 
                axis=1
            )
            metadata_s8 = pd.concat([a, b])
            
            # Merge with RPS data
            df = df.merge(
                metadata_s8[['genome1', 'genome2', 'carbonsource_count', 'interaction_class']],
                on=['genome1', 'genome2'],
                how='inner'
            )
            print(f"  After merging: {df.shape[0]} observations")
        except Exception as e:
            print(f"  Could not load carbon source count: {e}")
    
    if 'carbonsource_count' not in df.columns:
        print("  No carbon source count data available")
        return None, None
    
    # Calculate summary statistics
    mean_df = df.groupby(['carbonsource_count', 'interaction_class'])[value_column].agg(['mean', 'count']).reset_index()
    median_df = df.groupby(['carbonsource_count', 'interaction_class'])[value_column].agg(['median', 'count']).reset_index()
    
    # Save summary statistics
    mean_df.to_csv('csource_RPS_interaction_mean_count.tsv', sep='\t', index=False)
    median_df.to_csv('csource_RPS_interaction_median_count.tsv', sep='\t', index=False)
    
    # Create heatmaps
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Mean heatmap
    pivot_mean = mean_df.pivot(
        index='carbonsource_count',
        columns='interaction_class',
        values='mean'
    )
    
    sns.heatmap(pivot_mean, ax=axes[0], cmap='coolwarm', center=0,
                cbar_kws={'label': f'Mean {value_column}'})
    axes[0].set_title(f'Mean {value_column} by Carbon Source Count and Interaction Class')
    axes[0].set_xlabel('Interaction Class')
    axes[0].set_ylabel('Carbon Source Count')
    
    # Median heatmap
    pivot_median = median_df.pivot(
        index='carbonsource_count',
        columns='interaction_class',
        values='median'
    )
    
    sns.heatmap(pivot_median, ax=axes[1], cmap='coolwarm', center=0,
                cbar_kws={'label': f'Median {value_column}'})
    axes[1].set_title(f'Median {value_column} by Carbon Source Count and Interaction Class')
    axes[1].set_xlabel('Interaction Class')
    axes[1].set_ylabel('Carbon Source Count')
    
    plt.tight_layout()
    plt.savefig('heatmap_RPS_csources_interaction.png', dpi=300, bbox_inches='tight')
    plt.savefig('heatmap_RPS_csources_interaction.svg', bbox_inches='tight')
    plt.close()
    
    print("  Heatmaps saved as 'heatmap_RPS_csources_interaction.png' and '.svg'")
    
    return mean_df, median_df

def run_comprehensive_analysis(df):
    """
    Run comprehensive analysis of RPS vs interaction data.
    
    Parameters:
    -----------
    df : pandas.DataFrame
    
    Returns:
    --------
    dict
        Dictionary with analysis results
    """
    print("=" * 60)
    print("COMPREHENSIVE ANALYSIS OF RPS VS INTERACTION DATA")
    print("=" * 60)
    
    # Ensure RPS column exists (rename if EIT exists)
    if 'EIT' in df.columns and 'RPS' not in df.columns:
        print("Renaming 'EIT' column to 'RPS'...")
        df['RPS'] = df['EIT']
    
    if 'RPS' not in df.columns:
        print("ERROR: No RPS column found in data")
        return {}
    
    print("\n1. BASIC STATISTICS")
    print("-" * 40)
    
    # Count unique pairs
    unique_pairs = df.apply(
        lambda row: tuple(sorted([str(row['genome1']), str(row['genome2'])])), 
        axis=1
    ).unique()
    print(f"Total unique strain pairs: {len(unique_pairs)}")
    print(f"Total observations: {len(df)}")
    
    # Calculate overall statistics by interaction class
    if 'interaction_class' in df.columns:
        overall_stats = df.groupby('interaction_class')['RPS'].agg(['mean', 'std', 'median', 'count'])
        print("\nOverall RPS by interaction class:")
        for idx, row in overall_stats.iterrows():
            print(f"  {idx}: mean={row['mean']:.4f}, median={row['median']:.4f}, "
                  f"std={row['std']:.4f}, n={int(row['count'])}")
    
    print("\n2. CORRELATION ANALYSIS")
    print("-" * 40)
    
    # Calculate correlations
    corr_results, df = calculate_correlations(df, genomic_metric='RPS')
    
    if not corr_results.empty:
        print("\nCorrelation results:")
        print(corr_results.to_string())
        corr_results.to_csv('correlations_RPS_interaction.tsv', sep='\t', index=False)
    
    print("\n3. PERMUTATION ANOVA")
    print("-" * 40)
    
    f_obs, p_value, classes = permutation_anova_by_interaction_class(
        df, value_column='RPS', n_permutations=5000
    )
    
    if f_obs is not None:
        print(f"\nPermutation ANOVA Result:")
        print(f"  F = {f_obs:.4f}, p = {p_value:.4f}")
        if p_value < 0.05:
            print(f"  -> Significant difference between interaction classes (p < 0.05)")
        else:
            print(f"  -> No significant difference between interaction classes")
    
    print("\n4. RANDOM FOREST ANALYSIS")
    print("-" * 40)
    
    varimp = random_forest_analysis(df, target_column='RPS')
    
    if not varimp.empty:
        varimp.to_csv('random_forest_importances.tsv', sep='\t')
    
    print("\n5. SUMMARY")
    print("-" * 40)
    
    # Create pair summary
    pair_summary = prepare_pairwise_data(df, value_column='RPS')
    
    if not pair_summary.empty:
        print(f"\nPair-level summary (n={len(pair_summary)}):")
        print(pair_summary[['pair_id', 'modal_interaction', 'RPS', 'n_observations']].to_string())
        pair_summary.to_csv('pair_summary_RPS.tsv', sep='\t', index=False)
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    
    return {
        'correlations': corr_results,
        'pair_summary': pair_summary,
        'anova_result': {'f': f_obs, 'p': p_value, 'classes': classes},
        'feature_importances': varimp
    }

def main():
    """Main execution function."""
    # Configuration
    CARBONCOMP_PATH = 'carboncomp_output.tsv.xz'  # output from running CaCo in the genomes retrieved from Ono et al.
    EXCEL_PATH = 'supplementary-tables_wraf224.xlsx'  # Ono et al. supplementary information table
    
    print("=" * 60)
    print("ANALYSIS OF RPS VS INTERACTION DATA (Ono et al.)")
    print("=" * 60)
    
    # Load and prepare data
    df = load_and_prepare_data(CARBONCOMP_PATH, EXCEL_PATH)
    
    # Run comprehensive analysis
    results = run_comprehensive_analysis(df)
    
    # Create heatmaps if carbon source data is available
    create_interaction_heatmaps(df, value_column='RPS', excel_path=EXCEL_PATH)
    
    print("\nAll analyses completed successfully!")
    print(f"Output files have been saved to current directory.")
    
    return results

if __name__ == "__main__":
    main()
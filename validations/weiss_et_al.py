
#!/usr/bin/env python3
"""
Script for correlating CaCo-derived genomic metrics with experimental competition data.
Compares niche overlap and competition scores with experimentally measured interaction metrics.
"""

import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests

def load_and_prepare_data(carboncomp_path, excel_path):
    """
    Load and merge genomic and experimental data.
    
    Parameters:
    -----------
    carboncomp_path : str
        Path to carboncomp_output.tsv.xz file
    excel_path : str
        Path to Excel file with experimental data
    
    Returns:
    --------
    pandas.DataFrame
        Merged dataframe with genomic and experimental metrics
    """
    # Load genomic data
    print("Loading genomic data...")
    df = pd.read_table(carboncomp_path)
    
    # Clean genome identifiers
    df.genome1 = df.genome1.apply(lambda x: x.split('_')[-1].replace('.fasta', ''))
    df.genome2 = df.genome2.apply(lambda x: x.split('_')[-1].replace('.fasta', ''))
    
    # Load and merge experimental data from different sheets
    print("Loading experimental data...")
    
    # Sheet 1: SM depletion overlap
    df2 = pd.read_excel(excel_path, 
                       sheet_name='SM depletion overlap', 
                       index_col=0)
    df2 = df2.stack().reset_index()
    df2.columns = ['genome1', 'genome2', 'sm_depletion']
    df = df.merge(on=['genome1', 'genome2'], right=df2)
    
    # Sheet 2: dAUC SM
    df2 = pd.read_excel(excel_path, 
                       sheet_name='dAUC SM', 
                       index_col=0)
    df2 = df2.stack().reset_index()
    df2.level_1 = df2.level_1.apply(lambda x: x.replace('sm', ''))
    df2.columns = ['genome1', 'genome2', 'AUC_sm_depletion']
    df = df.merge(on=['genome1', 'genome2'], right=df2)
    
    # Sheet 3: dGR SM
    df2 = pd.read_excel(excel_path, 
                       sheet_name='dGR SM', 
                       index_col=0)
    df2 = df2.stack().reset_index()
    df2.level_1 = df2.level_1.apply(lambda x: x.replace('sm', ''))
    df2.columns = ['genome1', 'genome2', 'GR_sm_depletion']
    df = df.merge(on=['genome1', 'genome2'], right=df2)
    
    print(f"Final dataframe shape: {df.shape}")
    return df

def calculate_correlations(df, genomic_metrics=['RPS', 'competition'], 
                          experimental_metrics=['sm_depletion', 'AUC_sm_depletion', 'GR_sm_depletion']):
    """
    Calculate Spearman correlations between genomic and experimental metrics.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        Dataframe containing both genomic and experimental metrics
    genomic_metrics : list
        List of genomic metric column names
    experimental_metrics : list
        List of experimental metric column names
    
    Returns:
    --------
    pandas.DataFrame
        Dataframe with correlation results (r, p, q values)
    """
    results = []
    
    print("Calculating correlations...")
    for genomic in genomic_metrics:
        for experimental in experimental_metrics:
            # Remove NaN values for correlation
            valid_data = df[[genomic, experimental]].dropna()
            
            if len(valid_data) > 3:  # Minimum data points for correlation
                r, p = stats.spearmanr(valid_data[genomic], valid_data[experimental])
                results.append({
                    'genomic_metric': genomic,
                    'experimental_metric': experimental,
                    'r': r,
                    'p': p,
                    'n': len(valid_data)
                })
    
    results_df = pd.DataFrame(results)
    
    # Apply FDR correction
    _, q_values, _, _ = multipletests(results_df['p'], method='fdr_bh')
    results_df['q'] = q_values
    
    return results_df

def plot_correlations(results_df, output_path='correlation_results.png'):
    """
    Create visualization of correlation results.
    
    Parameters:
    -----------
    results_df : pandas.DataFrame
        Dataframe with correlation results
    output_path : str
        Path to save the output plot
    """
    plt.figure(figsize=(10, 6))
    
    # Create heatmap-style visualization
    pivot_table = results_df.pivot(index='genomic_metric', 
                                   columns='experimental_metric', 
                                   values='r')
    
    # Create annotation with r and q values
    annot_data = []
    for genomic in pivot_table.index:
        row = []
        for exp in pivot_table.columns:
            mask = (results_df['genomic_metric'] == genomic) & (results_df['experimental_metric'] == exp)
            if mask.any():
                r_val = results_df.loc[mask, 'r'].values[0]
                q_val = results_df.loc[mask, 'q'].values[0]
                row.append(f"r={r_val:.2f}\nq={q_val:.3f}")
            else:
                row.append('')
        annot_data.append(row)
    
    sns.heatmap(pivot_table, annot=annot_data, fmt='', 
                cmap='coolwarm', center=0, 
                cbar_kws={'label': 'Spearman ρ'})
    
    plt.title('Correlation between Genomic and Experimental Competition Metrics')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Plot saved to {output_path}")

def save_results(results_df, output_path='correlation_results.csv'):
    """
    Save correlation results to CSV file.
    
    Parameters:
    -----------
    results_df : pandas.DataFrame
        Dataframe with correlation results
    output_path : str
        Path to save the results
    """
    results_df.to_csv(output_path, index=False)
    print(f"Results saved to {output_path}")

def main():
    """Main execution function."""
    # Configuration - update these paths for your system
    CARBONCOMP_PATH = 'carboncomp_output.tsv.xz'  # output from running CaCo default in the genomes downloaded from Weiss et al.
    EXCEL_PATH = '41396_2021_1153_MOESM2_ESM.xlsx' # supplementary materials from Weiss et al.
    
    # Load data
    df = load_and_prepare_data(CARBONCOMP_PATH, EXCEL_PATH)
    
    # Calculate correlations
    results = calculate_correlations(df)
    
    # Display results
    print("\n=== Correlation Results ===")
    print(results.to_string())
    
    # Save results
    save_results(results, 'correlation_results.csv')
    
    # Create visualization
    plot_correlations(results, 'correlation_heatmap.png')
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()
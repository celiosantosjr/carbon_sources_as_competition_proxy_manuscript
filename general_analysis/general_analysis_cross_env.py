"""
Carbon Metabolism Analysis Pipeline
====================================
This pipeline analyzes carbon metabolism in microbial genomes from different environments.
It includes data preparation, CaCo analysis, and downstream statistical/visualization steps.

Pipeline Steps:
1. Data preparation and filtering
2. Downloading and processing MAGs
3. Running CaCo analysis for carbon competition
4. Protein feature calculations
5. Statistical analysis and visualization
6. Taxonomic analysis

Outputs are organized in ./output/ directory
"""

# ============================================
# CONFIGURATION
# ============================================
import os
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from glob import glob
from tqdm import tqdm
from subprocess import run
from itertools import chain, combinations
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from matplotlib.lines import Line2D
import tarfile
import shutil

# Create output directories
os.makedirs('output/plots', exist_ok=True)
os.makedirs('output/tables', exist_ok=True)
os.makedirs('genomes', exist_ok=True)
os.makedirs('tmp', exist_ok=True)

# Color maps for environments
ENV_COLORS = {
    'Ocean Water': '#d95f02',    # Orange
    'Freshwater': '#1b9e77',     # Teal
    'Soil': '#7570b3',           # Purple
    'Human Gut': '#e7298a'       # Pink
}

# ============================================
# STEP 1: DATA PREPARATION
# ============================================

def prepare_genome_metadata(input_file='data/mOTUs3.genome_metadata.tsv.gz',
                           output_file='sample_genomes_to_work.tsv'):
    """
    Filter and prepare genome metadata for analysis.
    
    Parameters:
    -----------
    input_file : str
        Path to input metadata file
    output_file : str
        Path to output filtered metadata file
    
    Returns:
    --------
    pd.DataFrame: Filtered genome metadata
    """
    print("Step 1: Preparing genome metadata...")
    
    df = pd.read_table(input_file, sep='\t', low_memory=False)
    
    # Quality filtering
    df['OVERALL_QUAL'] = df.CHECKM_COMPLETENESS - (5 * df.CHECKM_CONTAMINATION)
    df = df[df.OVERALL_QUAL >= 50]
    df = df[(df.CHECKM_COMPLETENESS >= 90) & (df.CHECKM_CONTAMINATION <= 5)]
    df = df[df['#GENOME'].apply(lambda x: '.' not in x)]
    df = df[df['GUNC_PROGENOMES2.1'] == True]
    df = df[~df.ENVIRONMENT.isna()]
    
    # Extract sample information
    df['SAMPLE'] = df['#GENOME'].apply(lambda x: x.split('_')[1])
    df = df[~df.SAMPLE.str.contains('GENOME')]
    
    # Filter for samples with at least 2 MAGs
    sample_counts = df.SAMPLE.value_counts()
    valid_samples = sample_counts[sample_counts >= 2].index
    df = df[df.SAMPLE.isin(valid_samples)]
    
    # Filter for target environments
    target_envs = ['Human Gut', 'Ocean Water', 'Soil', 'Freshwater']
    df = df[df.ENVIRONMENT.isin(target_envs)]
    
    print(f"Found {len(df)} genomes after filtering")
    print("Environment distribution:")
    print(df.ENVIRONMENT.value_counts())
    
    # Parse GTDB taxonomy
    df['GTDB-R207'] = df['GTDB-R207'].apply(lambda x: x.split(';'))
    taxonomic_levels = ['DOMAIN', 'PHYLUM', 'CLASS', 'ORDER', 'FAMILY', 'GENUS', 'SPECIES']
    
    for i, level in enumerate(taxonomic_levels):
        df[level] = df['GTDB-R207'].apply(lambda x: x[i].split('__')[1] if i < len(x) else '')
    
    df = df.drop('GTDB-R207', axis=1)
    
    # Save filtered data
    df.to_csv(output_file, sep='\t', index=False)
    
    # Plot sample distribution
    sample_env_counts = df[['SAMPLE', 'ENVIRONMENT']].drop_duplicates()
    env_counts = sample_env_counts.ENVIRONMENT.value_counts()
    
    plt.figure(figsize=(10, 6))
    env_counts.plot(kind='bar', color=[ENV_COLORS.get(env, 'gray') for env in env_counts.index])
    plt.title('Number of Samples per Environment')
    plt.xlabel('Environment')
    plt.ylabel('Number of Samples')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('output/plots/number_of_samples_per_environment.svg')
    plt.close()
    
    return df

# ============================================
# STEP 2: DOWNLOAD MAGs
# ============================================

def download_mags(metadata_file='sample_genomes_to_work.tsv'):
    """
    Download MAGs from provided links.
    
    Parameters:
    -----------
    metadata_file : str
        Path to metadata file with download links
    """
    print("\nStep 2: Downloading MAGs...")
    
    def check_downloaded(link, checker_df):
        """Check which files from a download link are already downloaded."""
        done_files = [name.replace('.fa.gz', '').replace('genomes/', '') 
                     for name in glob('genomes/*.fa.gz')]
        files = checker_df.loc[checker_df['DOWNLOAD_LINK'] == link, 'link'].tolist()[0]
        n = len(set(files).intersection(set(done_files)))
        return n == len(files)
    
    df = pd.read_table(metadata_file)
    download_groups = df.groupby('DOWNLOAD_LINK').apply(lambda x: x['#GENOME'].tolist())
    download_groups = download_groups.reset_index().rename({0: 'link'}, axis=1)
    download_groups['L'] = download_groups.link.apply(lambda x: len(x))
    download_groups = download_groups.sort_values(by='L', ascending=False).drop('L', axis=1)
    checker = download_groups.copy()
    
    # Download each group
    for _, link, genomes in tqdm(download_groups.itertuples(), total=len(download_groups)):
        if not check_downloaded(link, checker):
            print(f"Downloading from: {link}")
            run(['wget', link, '-O', 'temp.tar'])
            
            # Extract tar file
            with tarfile.open('temp.tar') as tar:
                tar.extractall()
            
            # Move genome files
            for genome in genomes:
                tar_dir = 'temp.tar'.replace('.tar', '')
                src_path = f"{tar_dir}/{genome}.fa.gz"
                if os.path.exists(src_path):
                    os.rename(src_path, f'genomes/{genome}.fa.gz')
            
            # Cleanup
            if os.path.exists(tar_dir):
                shutil.rmtree(tar_dir)
            os.remove('temp.tar')
    
    print(f"Downloaded {len(glob('genomes/*.fa.gz'))} genome files")

# ============================================
# STEP 3: GENOME DESCRIPTION
# ============================================

def describe_genomes(genome_dir='genomes', output_file='output/tables/genomes_description.tsv'):
    """
    Calculate basic statistics for each genome.
    
    Parameters:
    -----------
    genome_dir : str
        Directory containing genome FASTA files
    output_file : str
        Path to output description file
    """
    print("\nStep 3: Calculating genome statistics...")
    
    results = []
    for infile in tqdm(glob(f'{genome_dir}/*.fa*')):
        genome = os.path.basename(infile).replace('.fa', '').replace('.gz', '')
        length, g_count, c_count, contigs = 0, 0, 0, 0
        
        # Handle gzipped files
        open_func = gzip.open if infile.endswith('.gz') else open
        mode = 'rt' if infile.endswith('.gz') else 'r'
        
        with open_func(infile, mode) as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                seq = str(record.seq).upper()
                length += len(seq)
                g_count += seq.count('G')
                c_count += seq.count('C')
                contigs += 1
        
        gc_content = (g_count + c_count) / length if length > 0 else 0
        results.append((genome, length, gc_content, contigs))
    
    df = pd.DataFrame(results, columns=['genome', 'length', 'GC%', 'contigs'])
    df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Processed {len(df)} genomes")
    return df

# ============================================
# STEP 4: CaCo ANALYSIS
# ============================================

def run_caco_analysis(metadata_file='sample_genomes_to_work.tsv'):
    """
    Run CaCo analysis in batches for carbon competition prediction.
    
    Parameters:
    -----------
    metadata_file : str
        Path to metadata file
    """
    print("\nStep 4: Running CaCo analysis...")
    
    samples = pd.read_table(metadata_file)
    sample_groups = samples.groupby('SAMPLE').apply(lambda x: x['#GENOME'].tolist())
    
    # Create batch files
    for idx, genomes in enumerate(sample_groups):
        with open(f'tmp/caco_gl_batch_{idx}.txt', 'w') as f:
            for genome in genomes:
                f.write(f'genomes/{genome}.fa\n')
    
    # Create CaCo script
    batch_files = glob('tmp/caco_gl_batch_*.txt')
    with open('tmp/calling_CaCo.sh', 'w') as f:
        for infile in batch_files:
            batch = os.path.basename(infile).replace('caco_gl_', '').replace('.txt', '')
            f.write(f'python3 CaCo.py -m from_nucleotides -gl {infile}\n')
            f.write(f'mv allfams.tsv output/tables/allfams_{batch}.tsv\n')
            f.write(f'mv allsubs.tsv output/tables/allsubs_{batch}.tsv\n')
            f.write(f'mv carboncomp_output.tsv.xz output/tables/carboncomp_output_{batch}.tsv.xz\n')
            f.write(f'rm -rf tmp/\n\n')
        
        f.write(f'cat output/tables/allfams_batch*.tsv > output/tables/allfams.tsv\n')
        f.write(f'rm -rf output/tables/allfams_batch*.tsv\n\n')
        f.write(f'cat output/tables/allsubs_batch*.tsv > output/tables/allsubs.tsv\n')
        f.write(f'rm -rf output/tables/allsubs_batch*.tsv\n\n')
        f.write(f'xzcat output/tables/carboncomp_output_*.tsv.xz | xz > output/tables/carboncomp_output.tsv.xz\n')
        f.write(f'rm -rf output/tables/carboncomp_output_*.tsv.xz\n')
    
    # Run CaCo (commented out - uncomment to run)
    # print("Running CaCo analysis (this may take a while)...")
    # run(['bash', 'tmp/calling_CaCo.sh'])
    
    # Cleanup
    for infile in batch_files:
        os.remove(infile)
    os.remove('tmp/calling_CaCo.sh')
    
    print("CaCo analysis setup complete")

# ============================================
# STEP 5: DATA CLEANING
# ============================================

def clean_caco_outputs():
    """
    Clean and fix CaCo output tables.
    """
    print("\nStep 5: Cleaning CaCo outputs...")
    
    # Clean carboncomp_output
    df = pd.read_table('output/tables/carboncomp_output.tsv.xz', sep='\t')
    df = df[df.genome1 != 'genome1']  # Remove header rows
    df.rename(columns={'EIT': 'RPS'}, inplace=True)  # Rename EIT to RPS
    df.to_csv('output/tables/carboncomp_output.tsv.xz', sep='\t', index=False, compression='xz')
    
    # Clean allfams
    df = pd.read_table('output/tables/allfams.tsv', sep='\t')
    df = df[df.genome != 'genome']
    df.to_csv('output/tables/allfams.tsv', sep='\t', index=False)
    
    # Clean allsubs
    df = pd.read_table('output/tables/allsubs.tsv', sep='\t')
    df = df[df.genome != 'genome']
    df.to_csv('output/tables/allsubs.tsv', sep='\t', index=False)
    
    print("Output cleaning complete")

# ============================================
# STEP 6: COMPETITION ANALYSIS
# ============================================

def analyze_competition_by_environment():
    """
    Analyze competition scores by environment with statistical testing.
    """
    print("\nStep 6: Analyzing competition by environment...")
    
    df = pd.read_table('output/tables/carboncomp_output.tsv.xz')
    samples = pd.read_table('sample_genomes_to_work.tsv')
    env_dict = samples.set_index('#GENOME')['ENVIRONMENT'].to_dict()
    
    df['env'] = df.genome1.apply(lambda x: env_dict.get(x))
    df = df[df.prob < 0.05]  # Significant interactions only
    
    # Order environments by median competition
    env_order = df.groupby('env').apply(lambda x: x['competition'].median())
    env_order = env_order.sort_values().index
    
    # Create boxplot with swarmplot
    plt.figure(figsize=(12, 8))
    sns.boxplot(data=df, x='env', y='competition', showfliers=False, 
                color='white', order=env_order)
    
    # Sample for swarmplot to avoid overplotting
    swarm_sample = df.groupby('env').apply(lambda x: x.sample(min(2000, len(x))))
    swarm_sample = swarm_sample.reset_index(drop=True)
    sns.swarmplot(data=swarm_sample, x='env', y='competition', 
                  s=1.5, order=env_order, alpha=0.3)
    
    plt.xlabel('Environment', fontsize=12)
    plt.ylabel('Competition Score', fontsize=12)
    plt.title('Competition Scores by Environment', fontsize=14)
    plt.tight_layout()
    plt.savefig('output/plots/competition_by_environment.svg', dpi=300)
    plt.close()
    
    # Statistical testing between environments
    test_results = []
    for env1, env2 in combinations(env_order, 2):
        u_stat, p_val = mannwhitneyu(df[df.env == env1]['competition'],
                                    df[df.env == env2]['competition'])
        test_results.append((env1, env2, u_stat, p_val))
    
    tests_df = pd.DataFrame(test_results, columns=['env1', 'env2', 'U', 'p_value'])
    _, p_adj, _, _ = multipletests(tests_df['p_value'])
    tests_df['p_adj'] = p_adj
    
    tests_df.to_csv('output/tables/environment_comparison_competition.tsv', 
                   sep='\t', index=False)
    
    print("Competition analysis complete")
    return df

# ============================================
# STEP 7: DIVERSITY AND COMPETITION ANALYSIS
# ============================================

def analyze_diversity_competition():
    """
    Analyze relationship between diversity and competition.
    """
    print("\nStep 7: Analyzing diversity and competition relationships...")
    
    df = pd.read_table('output/tables/carboncomp_output.tsv.xz')
    samples = pd.read_table('sample_genomes_to_work.tsv')
    
    # Convert RPS column (formerly EIT) to binary
    df['RPS'] = df['RPS'].apply(lambda x: x >= 0)
    
    # Map sample and environment info
    sample_dict = samples.set_index('#GENOME')['SAMPLE'].to_dict()
    env_dict = samples.set_index('#GENOME')['ENVIRONMENT'].to_dict()
    
    df['sample'] = df.genome1.apply(lambda x: sample_dict.get(x))
    df['env'] = df.genome1.apply(lambda x: env_dict.get(x))
    
    # Calculate positive interactions per sample
    pos_interactions = df.groupby('sample').apply(lambda x: x['RPS'].sum())
    total_interactions = df.groupby('sample').size()
    
    diversity_df = pd.concat([pos_interactions, total_interactions], axis=1)
    diversity_df.columns = ['positive_interactions', 'total_interactions']
    diversity_df = diversity_df.reset_index()
    
    # Add environment info
    sample_env_dict = df.set_index('sample')['env'].to_dict()
    diversity_df['env'] = diversity_df['sample'].apply(lambda x: sample_env_dict.get(x))
    
    # Plot relationship
    plt.figure(figsize=(10, 8))
    for env, color in ENV_COLORS.items():
        env_data = diversity_df[diversity_df.env == env]
        plt.scatter(env_data['total_interactions'], env_data['positive_interactions'],
                   c=color, label=env, alpha=0.6, s=50)
        
        # Add trend line
        if len(env_data) > 1:
            z = np.polyfit(env_data['total_interactions'], 
                          env_data['positive_interactions'], 1)
            p = np.poly1d(z)
            plt.plot(env_data['total_interactions'], p(env_data['total_interactions']),
                    color=color, linestyle='--', linewidth=2)
    
    plt.xlabel('Total Interactions per Sample', fontsize=12)
    plt.ylabel('Positive Interactions (RPS ≥ 0)', fontsize=12)
    plt.title('Diversity-Competition Relationship', fontsize=14)
    plt.legend(title='Environment')
    plt.tight_layout()
    plt.savefig('output/plots/diversity_competition.svg', dpi=300)
    plt.close()
    
    # Calculate correlation for each environment
    corr_results = []
    for env in diversity_df['env'].unique():
        env_data = diversity_df[diversity_df.env == env]
        if len(env_data) > 2:
            r, p = spearmanr(env_data['total_interactions'],
                            env_data['positive_interactions'])
            corr_results.append((env, r, p))
    
    corr_df = pd.DataFrame(corr_results, columns=['environment', 'spearman_r', 'p_value'])
    corr_df.to_csv('output/tables/diversity_competition_correlation.tsv',
                  sep='\t', index=False)
    
    print("Diversity analysis complete")
    return diversity_df

# ============================================
# STEP 8: SUBSTRATE ANALYSIS
# ============================================

def analyze_substrates():
    """
    Analyze substrate utilization patterns.
    """
    print("\nStep 8: Analyzing substrate utilization...")
    
    subs_df = pd.read_table('output/tables/allsubs.tsv')
    samples = pd.read_table('sample_genomes_to_work.tsv')
    env_dict = samples.set_index('#GENOME')['ENVIRONMENT'].to_dict()
    
    # Count substrates per genome
    subs_df['substrate_count'] = subs_df['substrates'].apply(
        lambda x: len(str(x).split(', ')) if pd.notna(x) and str(x) != 'nan' else 0
    )
    subs_df['env'] = subs_df['genome'].apply(lambda x: env_dict.get(x))
    
    # Plot substrate distribution by environment
    plt.figure(figsize=(12, 8))
    for env, color in ENV_COLORS.items():
        env_data = subs_df[subs_df.env == env]
        if len(env_data) > 0:
            sns.kdeplot(data=env_data, x='substrate_count',
                       color=color, label=env, alpha=0.7, linewidth=2)
    
    plt.xlim(0, 60)
    plt.xlabel('Number of Substrates per MAG', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.title('Substrate Utilization Distribution', fontsize=14)
    plt.legend(title='Environment')
    plt.tight_layout()
    plt.savefig('output/plots/substrate_distribution.svg', dpi=300)
    plt.close()
    
    # Statistical testing between environments
    test_results = []
    envs = list(ENV_COLORS.keys())
    for env1, env2 in combinations(envs, 2):
        data1 = subs_df[subs_df.env == env1]['substrate_count']
        data2 = subs_df[subs_df.env == env2]['substrate_count']
        if len(data1) > 0 and len(data2) > 0:
            u_stat, p_val = mannwhitneyu(data1, data2)
            test_results.append((env1, env2, u_stat, p_val))
    
    if test_results:
        test_df = pd.DataFrame(test_results, columns=['env1', 'env2', 'U', 'p_value'])
        test_df.to_csv('output/tables/substrate_count_comparison.tsv',
                      sep='\t', index=False)
    
    # Analyze substrate-genome length relationship
    genome_info = pd.read_table('output/tables/genomes_description.tsv')
    merged_df = subs_df.merge(genome_info, on='genome', how='left')
    merged_df['length_mbp'] = merged_df['length'] / 1e6
    merged_df['gc_percent'] = merged_df['GC%'] * 100
    
    # Plot substrate count vs genome length
    plt.figure(figsize=(10, 8))
    for env, color in ENV_COLORS.items():
        env_data = merged_df[merged_df.env == env]
        if len(env_data) > 0:
            plt.scatter(env_data['substrate_count'], env_data['length_mbp'],
                       c=color, label=env, alpha=0.5, s=30)
    
    plt.xlabel('Number of Substrates', fontsize=12)
    plt.ylabel('Genome Length (Mbp)', fontsize=12)
    plt.title('Substrate Count vs Genome Length', fontsize=14)
    plt.legend(title='Environment')
    plt.tight_layout()
    plt.savefig('output/plots/substrates_vs_length.svg', dpi=300)
    plt.close()
    
    print("Substrate analysis complete")
    return merged_df

# ============================================
# STEP 9: NMDS AND PCA ANALYSIS
# ============================================

def perform_nmds_analysis():
    """
    Perform NMDS/PCA analysis on substrate data.
    """
    print("\nStep 9: Performing dimensionality reduction analysis...")
    
    data = pd.read_table('output/tables/allsubs.tsv')
    samples = pd.read_table('sample_genomes_to_work.tsv')
    
    # Parse substrates
    data['substrate_list'] = data['substrates'].apply(
        lambda x: str(x).split(', ') if pd.notna(x) and str(x) != 'nan' else []
    )
    
    # Create binary matrix
    all_substrates = set(chain.from_iterable(data['substrate_list']))
    binary_matrix = []
    
    for _, row in data.iterrows():
        binary_row = [1 if sub in row['substrate_list'] else 0 for sub in all_substrates]
        binary_row.insert(0, row['genome'])
        binary_matrix.append(binary_row)
    
    columns = ['genome'] + list(all_substrates)
    binary_df = pd.DataFrame(binary_matrix, columns=columns)
    binary_df = binary_df.set_index('genome').T
    
    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(binary_df.T)
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'])
    pca_df.index = binary_df.columns
    
    # Add environment information
    pca_df = pca_df.reset_index().rename(columns={'genome': '#GENOME'})
    pca_df = pca_df.merge(samples[['#GENOME', 'ENVIRONMENT']], on='#GENOME', how='left')
    
    # Plot PCA
    plt.figure(figsize=(12, 10))
    for env, color in ENV_COLORS.items():
        env_data = pca_df[pca_df.ENVIRONMENT == env]
        plt.scatter(env_data['PC1'], env_data['PC2'],
                   c=color, label=env, alpha=0.6, s=50)
    
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})', fontsize=12)
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})', fontsize=12)
    plt.title('PCA of Substrate Utilization Patterns', fontsize=14)
    plt.legend(title='Environment')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('output/plots/pca_substrates.svg', dpi=300)
    plt.close()
    
    print(f"PCA explained variance: PC1={pca.explained_variance_ratio_[0]:.3f}, "
          f"PC2={pca.explained_variance_ratio_[1]:.3f}")
    
    return pca_df

# ============================================
# MAIN PIPELINE FUNCTION
# ============================================

def run_pipeline(config=None):
    """
    Main pipeline function to run all analysis steps.
    
    Parameters:
    -----------
    config : dict, optional
        Configuration dictionary
    """
    print("=" * 60)
    print("CARBON METABOLISM ANALYSIS PIPELINE")
    print("=" * 60)
    
    # Default configuration
    if config is None:
        config = {
            'run_data_preparation': True,
            'run_download': False,  # Set to True to download MAGs
            'run_genome_description': True,
            'run_caco': False,  # Set to True to run CaCo
            'run_analysis': True,
            'run_protein_features': False
        }
    
    # Step 1: Data preparation
    if config['run_data_preparation']:
        metadata = prepare_genome_metadata()
    
    # Step 2: Download MAGs (optional)
    if config['run_download']:
        download_mags()
    
    # Step 3: Genome description
    if config['run_genome_description']:
        genome_stats = describe_genomes()
    
    # Step 4: CaCo analysis (optional)
    if config['run_caco']:
        run_caco_analysis()
        clean_caco_outputs()
    
    # Step 5-9: Analysis steps
    if config['run_analysis']:
        comp_df = analyze_competition_by_environment()
        div_df = analyze_diversity_competition()
        subs_df = analyze_substrates()
        pca_df = perform_nmds_analysis()
    
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETED SUCCESSFULLY")
    print("=" * 60)
    print(f"Outputs saved in:")
    print(f"  - output/plots/ (visualizations)")
    print(f"  - output/tables/ (data tables)")
    print(f"  - genomes/ (MAG sequences)")
    print("=" * 60)

# ============================================
# RUN PIPELINE
# ============================================

if __name__ == "__main__":
    # Configuration - modify based on needs
    config = {
        'run_data_preparation': True,
        'run_download': False,  # Requires wget and large downloads
        'run_genome_description': True,
        'run_caco': False,  # Requires CaCo installation
        'run_analysis': True,
        'run_protein_features': False
    }
    
    try:
        run_pipeline(config)
    except Exception as e:
        print(f"\nError in pipeline execution: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
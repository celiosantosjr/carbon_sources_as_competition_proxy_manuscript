# Robustness Analysis for CaCo Predictions

## **OVERVIEW**

This comprehensive robustness validation pipeline performs **six independent analyses** to test the robustness of CaCo (Carbon Competition) predictions against methodological biases and technical artifacts. The script validates that observed environmental patterns in carbon metabolism are not driven by sampling inequality, genome quality issues, annotation biases, or phylogenetic artifacts.

## **SIX VALIDATION ANALYSES PERFORMED**

1. **Bootstrap Subsampling** - Tests robustness to sampling inequality across environments
2. **Specialist/Generalist Cutoff Sensitivity** - Validates classification stability across different substrate count thresholds
3. **Genome Quality Correlation** - Tests for technical artifacts from genome assembly quality
4. **Substrate Prevalence Threshold Sensitivity** - Validates community metabolic profile clustering
5. **Annotation Bias Validation (4-part)** - Comprehensive test to rule out technical artifacts:
   - Heteroscedasticity-robust inference (HC3)
   - High-confidence genome filtering
   - Phylogenetic downsampling
   - Simulated annotation error
6. **Phylogenetic Pattern Robustness** - Tests Darwin's limiting similarity hypothesis

## **REQUIRED INPUT FILES:**

### **Core Input Files (from main pipeline):**
- **`carboncomp_output.tsv.xz`** - CaCo competition predictions (with `EIT` column renamed to `RPS`)
- **`allsubs.tsv`** - Substrate utilization predictions per genome
- **`allfams.tsv`** - Enzyme family predictions per genome
- **`sample_genomes_to_work.tsv`** - Filtered metadata with environment/taxonomy info
- **`genomes_description.tsv`** - Genome statistics (GC%, length, contigs)

### **File Naming Convention:**
The script automatically renames `EIT` → `RPS` (Relative Pairwise Score) and `relEIT` → `relRPS` for clarity.

## **DEPENDENCIES:**

### **Python Packages:**
```bash
# Core data analysis
pip install pandas numpy scipy statsmodels

# Visualization
pip install matplotlib seaborn

# Dimensionality reduction
pip install scikit-learn

# Progress bars
pip install tqdm

# Statistics
pip install statsmodels
```

### **System Dependencies:**
```bash
# Linux/Ubuntu
sudo apt-get install xz-utils

# macOS
brew install xz
```

## **FOLDER STRUCTURE:**
```
project/
├── robustness_analysis.py           # This script
├── data/                            # Input data directory
│   ├── carboncomp_output.tsv.xz     # CaCo competition predictions
│   ├── allsubs.tsv                  # Substrate utilization
│   ├── allfams.tsv                  # Enzyme families
│   ├── sample_genomes_to_work.tsv   # Genome metadata
│   └── genomes_description.tsv      # Genome statistics
├── robustness_outputs/              # Auto-created output directory
│   ├── plots/                       # All visualization files (.svg)
│   └── tables/                      # All data tables (.csv)
└── README.md                        # This file
```

## **CONFIGURATION:**

Edit the `CONFIG` dictionary at the top of the script:

```python
CONFIG = {
    # File paths (relative to script location)
    'competition_data': 'data/carboncomp_output.tsv.xz',
    'sample_metadata': 'data/sample_genomes_to_work.tsv',
    'genome_info': 'data/genomes_description.tsv',
    'substrate_data': 'data/allsubs.tsv',
    'family_data': 'data/allfams.tsv',
    'taxonomy_data': 'data/sample_genomes_to_work.tsv',
    
    # Analysis parameters
    'n_bootstrap_iterations': 1000,      # Reduce to 100 for quick testing
    'specialist_cutoffs': [20, 23, 25, 27, 30, 33],
    'thresholds': [0.5, 0.6, 0.75, 0.8, 0.9],
    'annotation_noise_levels': [0.1, 0.2, 0.3],
    
    # Output settings
    'output_dir': './robustness_outputs/',
    'seed': 42                           # For reproducibility
}
```

## **RUNNING THE ANALYSIS:**

### **Basic Run:**
```bash
python robustness_analysis.py
```

### **Run with Custom Configuration:**
```python
# Option 1: Edit CONFIG in the script
# Option 2: Run with different paths
import sys
sys.path.append('.')
from robustness_analysis import main

# Create custom configuration
custom_config = {
    'competition_data': 'path/to/your/carboncomp_output.tsv.xz',
    'sample_metadata': 'path/to/your/sample_genomes_to_work.tsv',
    'output_dir': './custom_outputs/',
    'n_bootstrap_iterations': 500  # Faster run
}

# Update config and run
import robustness_analysis
robustness_analysis.CONFIG.update(custom_config)
robustness_analysis.main()
```

## **EXPECTED OUTPUTS:**

After successful run, you'll find in `robustness_outputs/`:

### **Visualizations (SVG format):**
- **`bootstrap_subsampling_results.svg`** - Bootstrap analysis results (4 panels)
- **`cutoff_sensitivity_analysis.svg`** - Specialist cutoff sensitivity (4 panels)
- **`genome_quality_vs_substrates.svg`** - Quality-substrate correlations (2 panels)
- **`threshold_sensitivity_pca.svg`** - PCA across prevalence thresholds (6 panels)
- **`phylogenetic_patterns.svg`** - Phylogenetic structuring analysis (2 panels)

### **Data Tables (CSV format):**
- **`phylogenetic_analysis_results.csv`** - Statistical tests for phylogenetic patterns
- **`validation_summary.csv`** - Comprehensive summary statistics by environment

### **Console Output Summary:**
The script prints detailed progress and a final summary including:
1. Data loading statistics
2. Results for each of the six analyses
3. Key validation findings
4. Overall conclusions for manuscript

## **ANALYSIS DETAILS:**

### **1. Bootstrap Subsampling Analysis**
- **Purpose**: Test if environmental patterns persist after correcting for unequal sampling
- **Method**: 1000 bootstrap iterations, balanced resampling per environment
- **Output**: Median substrate counts and specialist proportions with 95% CIs

### **2. Specialist/Generalist Cutoff Sensitivity**
- **Purpose**: Validate that environmental rankings are stable across different specialist definitions
- **Method**: Test 6 cutoff values (20-33 substrates)
- **Output**: Spearman correlations and rank preservation statistics

### **3. Genome Quality Correlation**
- **Purpose**: Test if substrate predictions correlate with genome quality metrics
- **Method**: Quality score = Completeness - (5 × Contamination)
- **Output**: Spearman correlations by environment

### **4. Substrate Prevalence Threshold Sensitivity**
- **Purpose**: Validate community metabolic profile clustering stability
- **Method**: PCA of sample profiles at different prevalence thresholds (50-90%)
- **Output**: PCA plots and centroid movement metrics

### **5. Annotation Bias Validation (4-part)**
- **Purpose**: Rule out technical artifacts in phylogenetic patterns
- **Methods**:
  1. **HC3 robust standard errors** - Account for heteroscedasticity
  2. **High-confidence filtering** - Remove extreme CAZyme densities
  3. **Phylogenetic downsampling** - Balance phylum representation
  4. **Simulated annotation error** - Test robustness to 10-30% noise
- **Output**: Stability metrics and correlation coefficients

### **6. Phylogenetic Pattern Robustness**
- **Purpose**: Test Darwin's limiting similarity hypothesis
- **Method**: Intra-taxon vs inter-taxon niche overlap comparisons
- **Output**: Mann-Whitney U tests across taxonomic ranks (Phylum to Genus)

## **INTERPRETATION OF RESULTS:**

### **Key Metrics to Check:**
1. **Bootstrap CIs** - Narrow confidence intervals indicate robust patterns
2. **Spearman ρ > 0.8** - High correlation indicates cutoff stability
3. **Quality correlations < 0.3** - Low correlation suggests no technical artifacts
4. **Centroid movement < 0.5** - Low movement indicates threshold stability
5. **Phylogenetic p < 0.05** - Significant differences validate biological patterns

### **Success Criteria:**
- ✅ All six analyses complete without errors
- ✅ Bootstrap CIs do not cross zero for key comparisons
- ✅ Specialist proportions remain ranked similarly across cutoffs
- ✅ Genome quality correlations are weak (|r| < 0.3)
- ✅ PCA clusters remain stable across thresholds
- ✅ Phylogenetic patterns show significant structuring

## **TROUBLESHOOTING:**

### **Common Issues:**

1. **Missing Input Files:**
   ```
   FileNotFoundError: [Errno 2] No such file or directory: 'data/carboncomp_output.tsv.xz'
   ```
   **Solution**: Update file paths in `CONFIG` dictionary

2. **Memory Issues:**
   ```
   MemoryError: Unable to allocate array with shape...
   ```
   **Solution**: Reduce `n_bootstrap_iterations` to 100 or 500

3. **Missing Dependencies:**
   ```
   ModuleNotFoundError: No module named 'statsmodels'
   ```
   **Solution**: Install missing package: `pip install statsmodels`

4. **Empty Results:**
   ```
   No valid environments found for analysis
   ```
   **Solution**: Check that environment names match (`Ocean Water`, `Human Gut`, `Soil`, `Freshwater`)

5. **Column Renaming Warnings:**
   ```
   Column 'EIT' not found, using existing 'RPS' column
   ```
   **Solution**: This is informational - script handles both column names

### **Performance Tips:**
- For quick testing: Set `n_bootstrap_iterations = 100`
- Disable specific analyses by commenting out function calls in `main()`
- Use `thresholds = [0.75]` to test only the optimal threshold
- Process environments separately by modifying `valid_envs` list

## **CUSTOMIZATION OPTIONS:**

### **Modify Analysis Parameters:**
```python
# In the CONFIG dictionary:
'specialist_cutoffs': [15, 20, 25, 30, 35],  # Custom cutoffs
'thresholds': [0.7, 0.75, 0.8],              # Fewer thresholds
'annotation_noise_levels': [0.05, 0.1, 0.15], # Lower noise levels
'output_dir': './my_validation_results/',    # Custom output directory
```

### **Add New Environments:**
```python
'environment_colors': {
    'Freshwater': '#1b9e77',
    'Ocean Water': '#d95f02',
    'Soil': '#7570b3',
    'Human Gut': '#e7298a',
    'Sediment': '#66a61e',      # Add new environment
    'Wastewater': '#e6ab02'     # Add new environment
}
```

### **Run Specific Analyses Only:**
```python
# Comment out analyses you don't need in main() function:
# bootstrap_results = bootstrap_subsampling_analysis(...)  # Skip bootstrap
cutoff_results = analyze_cutoff_sensitivity(...)
quality_results = analyze_genome_quality_correlation(...)
# threshold_results = analyze_threshold_sensitivity(...)   # Skip threshold
annotation_results = perform_annotation_bias_validation(...)
phylogenetic_results = analyze_phylogenetic_patterns(...)
```

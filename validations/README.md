# Comprehensive Validation Suite for CaCo (Carbon Competition Predictions)

## **OVERVIEW**

This repository contains a comprehensive validation suite for the **CaCo (Carbon Competition)** tool, which predicts carbon source competition and niche partitioning from genomic data. The suite includes **four independent validation pipelines** that compare genomic predictions with experimental data across multiple studies and methodologies.

## **VALIDATION PIPELINES**

### **1. BIOLOG Phenotype Validation** (`caco_biolog_validation_pipeline.py`)
**Validates:** Genome-predicted carbon utilization against experimental BIOLOG™ PM1 phenotype data  
**Dataset:** 10 bacterial isolates from Colombian Amazon (river, soil, sediment)  
**Key Analysis:** Compares substrate counts and niche overlap between genomic predictions and empirical measurements

### **2. Coello-Camba et al. Perturbation Experiment** (`coello-camba_et_all.py`)
**Validates:** RPS (Resource Partitioning Score) across nutrient perturbation treatments  
**Dataset:** Marine mesocosm experiments with N, P, Si additions  
**Key Analysis:** Tests how nutrient availability affects predicted competition dynamics

### **3. Ono et al. Interaction Validation** (`ono_et_al.py`)
**Validates:** RPS predictions against experimentally measured interaction types  
**Dataset:** 8 bacterial strains with defined interaction classes (Competition, Mutualism, etc.)  
**Key Analysis:** Correlation between genomic competition scores and experimental interaction outcomes

### **4. Weiss et al. Competition Validation** (`weiss_et_al.py`)
**Validates:** Niche overlap and competition scores against growth measurements  
**Dataset:** Co-culture experiments with growth rate and AUC measurements  
**Key Analysis:** Spearman correlations between genomic metrics and experimental competition outcomes

## **QUICK START**

### **1. Install Dependencies:**

```
# Core packages
pip install pandas numpy scipy statsmodels scikit-learn matplotlib seaborn tqdm

# Additional for specific scripts
pip install openpyxl lxml  # For Excel file reading
pip install xlrd           # Legacy Excel support
```

### **2. Run All Validations (Complete Pipeline):**

```
# 1. First run CaCo on your genomes
python CaCo.py -m from_nucleotides -gl genome_list.txt -o carboncomp_output.tsv.xz

# 2. Run BIOLOG validation
python caco_biolog_validation_pipeline.py -g genome_list.txt -b csource.tsv -o biolog_validation

# 3. Run perturbation experiment validation (if data available)
python coello-camba_et_all.py

# 4. Run interaction validation (if data available)
python ono_et_al.py

# 5. Run competition validation (if data available)
python weiss_et_al.py
```

## **DETAILED SCRIPT DOCUMENTATION**

### **1. `caco_biolog_validation_pipeline.py`**

#### **Purpose:**
Validates CaCo predictions against experimental BIOLOG phenotype data from environmental isolates.

#### **Input Files:**
- `genome_list.txt` - List of genome FASTA file paths
- `csource.tsv` - BIOLOG PM1 phenotype data (binary matrix)
- (Optional) `data/substrate_key.json` - Mapping of substrates to CAZy families  # It should be used within CaCo installation

#### **Outputs:**
- `caco_biolog_validation_figure.png` - Four-panel validation figure
- `validation_summary_report.txt` - Detailed statistics and interpretation
- `genome_predictions/` - CaCo results from genome analysis
- `biolog_predictions/` - RPS calculations from BIOLOG data
- `comparison_results/` - Side-by-side comparison data

#### **Figure Panels:**
1. **Panel A:** Genome-predicted vs phenotype substrate counts (Spearman correlation)
2. **Panel B:** CaCo-predicted vs phenotypic niche overlap
3. **Panel C:** Niche overlap distributions from genomes (specialist-specialist, generalist-generalist, mixed)
4. **Panel D:** Niche overlap distributions from phenotypes

#### **Usage:**

```
python caco_biolog_validation_pipeline.py -g genome_list.txt -b csource.tsv -o output_dir
```

#### **Key Features:**

- Automatic isolate name matching between genome files and BIOLOG data
- Complete EIT/RPS calculation on both datasets using same formulas
- Proper terminology (RPS instead of EIT)
- Robust error handling with fallback to mock data
- Comprehensive statistics and publication-ready figures

---

### **2. `coello-camba_et_all.py`**

#### **Purpose:**

Analyzes CaCo-derived RPS across nutrient perturbation experiments (N, P, Si additions).

#### **Input Files:**

- `carboncomp_output.tsv.xz` - CaCo predictions from Coello-Camba et al. genomes
- `data/refactored_metadata.tsv` - Metadata for filtered MAGs

#### **Outputs:**

- `RPS_summary_statistics_per_condition.tsv` - Summary statistics by treatment
- `pairwise_mannwhitney_results.tsv` - Statistical comparisons between conditions
- `environmental_correlations.tsv` - Correlations with environmental variables
- `random_forest_feature_importances.tsv` - Important environmental predictors
- Multiple visualization files (KDE, box-swarm, violin plots)

#### **Statistical Tests:**

1. **Kruskal-Wallis** - Differences across all conditions
2. **Pairwise Mann-Whitney U** - Specific condition comparisons (FDR-corrected)
3. **Levene's test** - Homogeneity of variances
4. **Spearman correlations** - RPS vs environmental variables
5. **Random Forest** - Identify important environmental predictors

#### **Usage:**

```
python coello-camba_et_all.py
```

#### **Experimental Conditions Analyzed:**

- N,P Day 1 single addition
- N,P 2 week continuous addition
- Control mesocosm
- N,P,Si Day 1 single addition
- N,P,Si 2 week continuous addition

---

### **3. `ono_et_al.py`**

#### **Purpose:**

Correlates CaCo-predicted RPS with experimentally measured interaction types.

#### **Input Files:**

- `carboncomp_output.tsv.xz` - CaCo predictions from Ono et al. genomes
- `supplementary-tables_wraf224.xlsx` - Ono et al. supplementary tables

#### **Outputs:**

- `correlations_RPS_interaction.tsv` - Correlation results
- `pair_summary_RPS.tsv` - Summary by unique strain pairs
- `heatmap_RPS_csources_interaction.png` - Heatmap visualization
- `random_forest_importances.tsv` - Feature importances

#### **Interaction Classes:**

- Competition
- Exploitation
- Neutrality
- Mutualism

#### **Key Analyses:**

1. **Spearman correlations** - RPS vs interaction classes
2. **Permutation ANOVA** - Statistical differences between interaction classes
3. **Random Forest** - Identify predictive features
4. **Heatmap visualization** - RPS by carbon source count and interaction class

#### **Usage:**

```
python ono_et_al.py
```

#### **Dataset Characteristics:**

- 8 bacterial strains
- 28 unique pairwise combinations
- Experimentally determined interaction types
- Carbon source count data

---

### **4. `weiss_et_al.py`**

#### **Purpose:**

Correlates CaCo metrics with experimental competition measurements from co-culture studies.

#### **Input Files:**

- `carboncomp_output.tsv.xz` - CaCo predictions from Weiss et al. genomes
- `41396_2021_1153_MOESM2_ESM.xlsx` - Supplementary materials from Weiss et al.

#### **Outputs:**

- `correlation_results.csv` - Correlation results table
- `correlation_heatmap.png` - Visualization of correlations

#### **Experimental Metrics:**

- `sm_depletion` - Substrate depletion overlap
- `AUC_sm_depletion` - Area under curve for substrate depletion
- `GR_sm_depletion` - Growth rate changes

#### **Genomic Metrics:**

- `RPS` - Resource Partitioning Score (formerly EIT)
- `competition` - Direct competition score

#### **Statistical Analysis:**

- **Spearman correlations** between all genomic and experimental metrics
- **FDR correction** for multiple testing
- **Heatmap visualization** of correlation matrix

#### **Usage:**

```
python weiss_et_al.py
```

## **COMMON FILE FORMATS**

### **Input File Specifications:**

#### **1. Genome List (`genome_list.txt`):**

```
genomes/20C_canu_polished.fa
genomes/22C_flye_polished.fa
genomes/23C_canu_polished.fa
...
```

#### **2. BIOLOG Data (`csource.tsv`):**

Tab-separated matrix with substrates as rows and isolates as columns:

```
Csource	7C	20C	TP21	22C	...
D-Serina	0	1	0	0	...
D-Glucosa-6-fosfato	0	1	0	1	...
...
```

#### **3. CaCo Output (`carboncomp_output.tsv.xz`):**

Compressed tab-separated file with pairwise competition metrics:

```
genome1	genome2	set1	set2	intersection	competition	relcomp	prob	RPS	relRPS
...
```

#### **4. Substrate Key (`data/substrate_key.json`):**

JSON mapping of substrates to CAZy families:

```json
{
  "D-Glucose": ["GH1", "GH3", "GH5"],
  "L-Arabinose": ["GH43", "GH51"],
  ...
}
```

## **INSTALLATION AND SETUP**

### **Complete Installation:**

```
# Clone or download the repository
git clone <repository-url>
cd caco-validation-suite

# Create virtual environment (optional but recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install all dependencies
pip install -r requirements.txt

# Install additional dependencies if needed
pip install openpyxl lxml xlrd
```

### **Requirements File (`requirements.txt`):**

```
# Core data analysis
pandas>=1.5.0
numpy>=1.23.0
scipy>=1.9.0
statsmodels>=0.13.0
scikit-learn>=1.1.0

# Visualization
matplotlib>=3.6.0
seaborn>=0.12.0

# Progress and utilities
tqdm>=4.64.0
lzma  # Usually included with Python

# Excel file handling
openpyxl>=3.0.0
xlrd>=2.0.0
lxml>=4.9.0
```

## **WORKFLOW EXAMPLES**

### **Example 1: Complete Validation from Scratch**

```
# Step 1: Run CaCo on your genomes
python CaCo.py -m from_nucleotides -gl genome_list.txt -o carboncomp_output.tsv.xz

# Step 2: Validate against BIOLOG data
python caco_biolog_validation_pipeline.py -g genome_list.txt -b csource.tsv -o validation_results

# Step 3: Check outputs
ls validation_results/
# Should see: caco_biolog_validation_figure.png, validation_summary_report.txt, etc.
```

### **Example 2: Using Existing CaCo Results**

```
# If you already have CaCo results, just run validation:
python caco_biolog_validation_pipeline.py -g genome_list.txt -b csource.tsv -o validation_results --skip_caco
```

### **Example 3: Batch Processing Multiple Studies**

```
# Create a batch script (run_all_validations.sh)
#!/bin/bash

echo "=== Running BIOLOG Validation ==="
python caco_biolog_validation_pipeline.py -g genome_list.txt -b csource.tsv -o biolog_results

echo "=== Running Coello-Camba Analysis ==="
python coello-camba_et_all.py

echo "=== Running Ono et al. Analysis ==="
python ono_et_al.py

echo "=== Running Weiss et al. Analysis ==="
python weiss_et_al.py

echo "=== All validations complete ==="
```

## **EXPECTED RESULTS AND INTERPRETATION**

### **Successful Validation Indicators:**

1. **BIOLOG Validation:**
   - **Substrate count correlation:** Spearman ρ > 0.7, P < 0.05
   - **Niche overlap correlation:** Spearman ρ > 0.6, P < 0.05
   - **Consistent specialist-generalist patterns** in both genomic and phenotypic data

2. **Coello-Camba Analysis:**
   - **Significant differences** in RPS across nutrient treatments
   - **Interpretable environmental correlations** (e.g., RPS vs nutrient concentrations)

3. **Ono et al. Analysis:**
   - **Significant correlation** between RPS and interaction classes
   - **Higher RPS** for competitive pairs, lower for mutualistic pairs

4. **Weiss et al. Analysis:**
   - **Negative correlation** between RPS and growth overlap (expected: more competition = less coexistence)
   - **Consistent patterns** across different experimental metrics

### **Troubleshooting Common Issues:**

1. **"Column 'EIT' not found":**
   - **Solution:** Scripts automatically rename EIT → RPS, but ensure your CaCo output has either column

2. **Missing input files:**
   - **Solution:** Update file paths in script configuration or use command-line arguments

3. **Memory issues with large datasets:**
   - **Solution:** Reduce bootstrap iterations or use subset of data for testing

4. **Excel file reading errors:**
   - **Solution:** Ensure `openpyxl` and `xlrd` are installed, check Excel file format

## **CUSTOMIZATION AND EXTENSION**

### **Customizing Analysis Parameters:**

Each script has configurable parameters at the top. For example, in `coello-camba_et_all.py`:

```
# Configuration section
CARBONCOMP_PATH = 'carboncomp_output.tsv.xz'
METADATA_PATH = 'data/refactored_metadata.tsv'
```

### **Adding New Datasets:**

To validate CaCo with new experimental data:

1. **Format your data** to match existing templates
2. **Create a new script** following the structure of existing validations
3. **Add to the pipeline** by including in batch scripts

### **Modifying Visualization:**

All scripts use `matplotlib` and `seaborn` for plotting. Customize:
- Color schemes in `environment_colors` dictionaries
- Figure sizes and DPI settings
- Plot styles and annotations
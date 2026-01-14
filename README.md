# CaCo Validation and Analysis Suite

## **Overview**

This repository contains a comprehensive suite for validating and analyzing **CaCo (Carbon Competition)** predictions — a tool for predicting carbon source competition and niche partitioning from genomic data. The suite is organized into three main modules, each addressing different aspects of validation and analysis.

## **Repository Structure**

```
caco-validation-suite/
├── general_analysis/        # Main pipeline for environmental pattern analysis
├── robustness_analysis/     # Tests robustness to methodological biases  
├── validations/             # Validation against experimental datasets
├── README.md               # This file
└── requirements.txt        # Dependencies for all modules
```

## **1. General Analysis (`general_analysis/`)**

**Purpose:** Main pipeline for analyzing carbon metabolism patterns across environments using mOTUs database.

**Workflow:**
1. Data preparation from mOTUs metadata
2. MAG downloading and processing
3. CaCo execution on environmental genomes
4. Statistical analysis of competition patterns
5. Visualization of environmental differences

**Input:** Requires `mOTUs3.genome_metadata.tsv.gz` from the mOTUs database.

**Usage:**
```bash
cd general_analysis
python carbon_metabolism_pipeline.py
```

## **2. Robustness Analysis (`robustness_analysis/`)**

**Purpose:** Tests CaCo predictions against methodological biases and technical artifacts.

**Key Analyses:**
- Bootstrap subsampling for sampling inequality
- Specialist/generalist cutoff sensitivity
- Genome quality correlation tests
- Annotation bias validation (4-part)
- Phylogenetic pattern robustness

**Input:** Requires CaCo outputs (`carboncomp_output.tsv.xz`, `allsubs.tsv`, `allfams.tsv`) and genome metadata.

**Usage:**
```bash
cd robustness_analysis
python robustness_analysis.py
```

## **3. Validations (`validations/`)**

**Purpose:** Validates CaCo predictions against multiple experimental datasets.

### **Three Experimental Validations:**

1. **BIOLOG Phenotype Validation** - Compares predictions with experimental BIOLOG™ PM1 data
2. **Nutrient Perturbation Experiments** (Coello-Camba et al.) - Tests RPS across nutrient treatments
3. **Strain Interaction Validation** (Ono et al.) - Correlates RPS with experimental interaction types
4. **Competition Validation** (Weiss et al.) - Compares with growth competition measurements

**Input:** Each validation requires study-specific data files (see individual READMEs).

**Usage:**
```bash
cd validations
# Run specific validation
python caco_biolog_validation_pipeline.py -g genome_list.txt -b csource.tsv -o output_dir
python coello-camba_et_all.py
python ono_et_al.py
python weiss_et_al.py
```

## **Quick Start**

### **Installation:**
```bash
# Clone repository
git clone <repository-url>
cd caco-validation-suite

# Install dependencies
pip install -r requirements.txt

# Install CaCo (required for full analysis)
git clone https://github.com/amanbasu/CaCo
cd CaCo
# Follow CaCo installation instructions
```

### **Basic Usage Examples:**

1. **Start with general analysis:**
```bash
cd general_analysis
python carbon_metabolism_pipeline.py
```

2. **Validate with experimental data:**
```bash
cd validations
python caco_biolog_validation_pipeline.py -g ../data/genome_list.txt -b ../data/csource.tsv -o results
```

3. **Test robustness of results:**
```bash
cd robustness_analysis
python robustness_analysis.py
```

## **Dependencies**

**Core Requirements:**
```bash
pip install pandas numpy scipy statsmodels scikit-learn matplotlib seaborn tqdm
```

**Additional for specific modules:**
```bash
# For validations with Excel files
pip install openpyxl lxml xlrd

# For general analysis
pip install biopython gdown wget

# For robustness analysis
pip install statsmodels
```

## **Expected Outputs**

### **Each module generates:**
- **Publication-ready figures** (SVG/PNG formats)
- **Statistical summary tables** (TSV/CSV formats)
- **Detailed reports** with interpretation

### **Key Validation Metrics:**
- **Substrate count correlation:** ρ > 0.7 indicates good prediction accuracy
- **Niche overlap correlation:** ρ > 0.6 validates competition predictions
- **Significant phylogenetic patterns:** p < 0.05 supports biological interpretation
- **Robustness to methodology:** Stable results across analyses confirm reliability

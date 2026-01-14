The pipeline requires the following inputs and dependencies to run:

## **REQUIRED INPUT FILES:**

### **1. Core Input Files:**
- **`mOTUs3.genome_metadata.tsv.gz`** - Metadata file containing genome information with columns:
  - `#GENOME`: Genome identifier
  - `CHECKM_COMPLETENESS`: Completeness percentage
  - `CHECKM_CONTAMINATION`: Contamination percentage
  - `GUNC_PROGENOMES2.1`: Boolean for GUNC filtering
  - `ENVIRONMENT`: Environment classification
  - `DOWNLOAD_LINK`: URL to download genome
  - `GTDB-R207`: Taxonomy string

### **2. Optional Input Files:**
- **If skipping download step**: Genome FASTA files in `genomes/` directory named as `{GENOME_ID}.fa` or `{GENOME_ID}.fa.gz`

## **DEPENDENCIES:**

### **Python Packages:**
```bash
# Core data analysis
pip install pandas numpy scipy statsmodels

# Visualization
pip install matplotlib seaborn

# Bioinformatics
pip install biopython

# Utilities
pip install tqdm scikit-learn

# File handling
pip install gdown wget
```

### **System Dependencies:**
```bash
# Linux/Ubuntu
sudo apt-get install wget tar xz-utils

# macOS
brew install wget gnu-tar xz
```

### **External Tools:**
- **CaCo (Carbon Competition tool)**: Required if running CaCo analysis
  - Clone from: `https://github.com/amanbasu/CaCo`
  - Install dependencies as per CaCo documentation

## **FOLDER STRUCTURE:**
```
project/
├── mOTUs3.genome_metadata.tsv.gz    # REQUIRED INPUT
├── carbon_metabolism_pipeline.py     # This script
├── genomes/                          # Auto-created for downloaded genomes
├── output/                           # Auto-created for results
│   ├── plots/                        # All visualizations
│   └── tables/                       # All data tables
└── tmp/                              # Temporary files
```

## **CONFIGURATION OPTIONS:**

The pipeline can be run with different configurations:

```python
# Minimal analysis (no downloads, no CaCo)
config = {
    'run_data_preparation': True,      # REQUIRES: mOTUs3.genome_metadata.tsv.gz
    'run_download': False,             # Skip MAG download
    'run_genome_description': False,   # Skip if genomes already described
    'run_caco': False,                 # Skip CaCo analysis
    'run_analysis': True,              # Run statistical analysis
}

# Full analysis (requires internet and CaCo)
config = {
    'run_data_preparation': True,
    'run_download': True,              # REQUIRES: Internet connection, 50-100GB disk
    'run_genome_description': True,
    'run_caco': True,                  # REQUIRES: CaCo installation
    'run_analysis': True,
}
```

## **PREPARATION STEPS:**

1. **Obtain the metadata file:**
   - Download `mOTUs3.genome_metadata.tsv.gz` from the mOTUs database
   - Or use the provided example if available

2. **Set up the environment:**
   ```bash
   # Create virtual environment
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   
   # Install dependencies
   pip install -r requirements.txt
   ```

3. **Prepare configuration:**
   - Edit the `config` dictionary at the bottom of the script
   - Set paths if different from defaults

## **RUNNING THE PIPELINE:**

```bash
# Basic run (requires only metadata file)
python carbon_metabolism_pipeline.py

# Or with custom configuration
python -c "
import sys
sys.path.append('.')
from carbon_metabolism_pipeline import run_pipeline
config = {'run_data_preparation': True, 'run_analysis': True}
run_pipeline(config)
"
```

## **EXPECTED OUTPUTS:**

After successful run, you'll find in `output/`:

### **Plots:**
- `competition_by_environment.svg` - Boxplots of competition scores
- `diversity_competition.svg` - Scatter plots of diversity vs competition
- `substrate_distribution.svg` - KDE plots of substrate counts
- `pca_substrates.svg` - PCA plot of substrate utilization
- `number_of_samples_per_environment.svg` - Bar plot of sample distribution

### **Tables:**
- `sample_genomes_to_work.tsv` - Filtered genome metadata
- `genomes_description.tsv` - Genome statistics (GC%, length, contigs)
- `environment_comparison_competition.tsv` - Statistical test results
- `diversity_competition_correlation.tsv` - Correlation results

## **TROUBLESHOOTING:**

1. **Missing metadata file:**
   ```
   FileNotFoundError: [Errno 2] No such file or directory: 'mOTUs3.genome_metadata.tsv.gz'
   ```
   **Solution:** Download from mOTUs or adjust path in `prepare_genome_metadata()` call.

2. **Missing dependencies:**
   ```
   ModuleNotFoundError: No module named 'pandas'
   ```
   **Solution:** Install missing packages: `pip install pandas`

3. **Disk space issues during download:**
   ```
   OSError: [Errno 28] No space left on device
   ```
   **Solution:** Need 50-100GB free space for genomes, or skip download step.

4. **CaCo not found:**
   ```
   FileNotFoundError: [Errno 2] No such file or directory: 'CaCo.py'
   ```
   **Solution:** Install CaCo or set `run_caco: False` in config.

The pipeline is designed to be modular - you can run individual steps by importing and calling specific functions if you don't want to run the entire pipeline.


## Workflow summarized


Here's the pipeline workflow in table format:

| Step | Input Files | Output Files | Dependencies/Tools | Notes |
|------|-------------|--------------|-------------------|-------|
| **1. Data Preparation** | `mOTUs3.genome_metadata.tsv.gz` (raw metadata from mOTUs) | `sample_genomes_to_work.tsv` (filtered metadata) | Python, pandas | Filters by quality, environment, and sample size |
| **2. Download MAGs** | `sample_genomes_to_work.tsv` (from Step 1, contains DOWNLOAD_LINK column) | `genomes/*.fa.gz` (downloaded genome files) | wget, tar, internet connection | Downloads MAGs from provided URLs, requires ~50-100GB disk space |
| **3. Genome Description** | `genomes/*.fa*` (from Step 2) | `genomes_description.tsv` (genome statistics) | BioPython (SeqIO), gzip (for compressed files) | Calculates GC%, length, contig counts for each genome |
| **4. CaCo Analysis** | 1. `sample_genomes_to_work.tsv` (from Step 1)<br>2. `genomes/*.fa` (decompressed from Step 2) | 1. `carboncomp_output.tsv.xz`<br>2. `allfams.tsv`<br>3. `allsubs.tsv` | CaCo.py tool, Python 3 | Requires CaCo installation; generates carbon competition predictions |
| **5. Data Cleaning** | 1. `carboncomp_output.tsv.xz`<br>2. `allfams.tsv`<br>3. `allsubs.tsv` (all from Step 4) | Cleaned versions of same files (with `EIT` renamed to `RPS`) | Python, pandas | Removes header rows, renames EIT column to RPS |
| **6. Competition Analysis** | 1. `carboncomp_output.tsv.xz` (cleaned, from Step 5)<br>2. `sample_genomes_to_work.tsv` (from Step 1) | 1. `competition_by_environment.svg`<br>2. `environment_comparison_competition.tsv` | Python, seaborn, matplotlib, scipy | Statistical tests comparing competition between environments |
| **7. Diversity & Competition Analysis** | 1. `carboncomp_output.tsv.xz` (cleaned, from Step 5)<br>2. `sample_genomes_to_work.tsv` (from Step 1) | 1. `diversity_competition.svg`<br>2. `diversity_competition_correlation.tsv` | Python, seaborn, scipy | Analyzes relationship between diversity and positive interactions |
| **8. Substrate Analysis** | 1. `allsubs.tsv` (cleaned, from Step 5)<br>2. `sample_genomes_to_work.tsv` (from Step 1)<br>3. `genomes_description.tsv` (from Step 3) | 1. `substrate_distribution.svg`<br>2. `substrate_count_comparison.tsv`<br>3. `substrates_vs_length.svg` | Python, seaborn, scipy | Analyzes substrate utilization patterns and correlations |
| **9. NMDS/PCA Analysis** | 1. `allsubs.tsv` (cleaned, from Step 5)<br>2. `sample_genomes_to_work.tsv` (from Step 1) | 1. `pca_substrates.svg`<br>2. `pca_substrates.tsv` | Python, scikit-learn, seaborn | Dimensionality reduction on substrate utilization data |

## **Pipeline Dependency Flow:**

```
mOTUs3.genome_metadata.tsv.gz
        ↓ (Step 1)
sample_genomes_to_work.tsv
        ├─────→ (Step 2) → genomes/*.fa* → (Step 3) → genomes_description.tsv
        │                         ↓ (decompress)
        │                   genomes/*.fa
        │                         ↓
        │                   (Step 4: CaCo) 
        ↓                   (requires both inputs)
carboncomp_output.tsv.xz    allfams.tsv    allsubs.tsv
        ↓ (Step 5)              ↓ (Step 5)       ↓ (Step 5)
 cleaned versions        cleaned versions   cleaned versions
        │                         │                │
        └─────→ (Step 6 & 7)     │                │
        │                │        │                │
        │                │        └─────→ (Step 8) │
        │                │                │        │
        │                │                │        └─────→ (Step 9)
        ↓                ↓                ↓                ↓
   competition    diversity        substrate       dimensionality
   analysis       analysis         analysis        reduction
```

## **Optional Steps Configuration:**

| Step | Required? | Can Skip If... |
|------|-----------|----------------|
| Step 2 (Download) | Optional | You already have genomes in `genomes/` directory |
| Step 3 (Description) | Optional | You already have `genomes_description.tsv` |
| Step 4 (CaCo) | Optional | You already have CaCo output files |
| Steps 5-9 | Required | For full analysis |

## **Minimal Input Set for Partial Runs:**

| To Run These Steps | Need These Files |
|-------------------|------------------|
| Steps 1, 6, 7 | 1. `mOTUs3.genome_metadata.tsv.gz`<br>2. `carboncomp_output.tsv.xz` (from elsewhere) |
| Steps 8, 9 | 1. `allsubs.tsv`<br>2. `sample_genomes_to_work.tsv`<br>3. `genomes_description.tsv` |
| Complete pipeline | Only `mOTUs3.genome_metadata.tsv.gz` (will download/process everything else) |

The pipeline is modular, allowing you to start from different points depending on what data you already have available.
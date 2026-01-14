## Key Features of This Pipeline:

### 1. **Complete EIT/RPS Calculation on Both Datasets:**
- **For genomes:** Runs the full CaCo pipeline (or custom implementation)
- **For BIOLOG:** Processes the csource.tsv file to generate `allsubs_biolog.tsv` and calculates RPS (formerly EIT) using the same formulas

### 2. **Proper Terminology:**
- Uses **RPS** (Relative Pairwise Score) instead of EIT
- `relRPS` instead of `relEIT`
- All outputs use consistent terminology

### 3. **Generates All Required Outputs:**
- Genome predictions (allsubs.tsv, carboncomp_output.tsv.xz)
- BIOLOG predictions (allsubs_biolog.tsv, caco_biolog_results.tsv.xz)
- Comparison data (substrate counts, overlap comparisons)
- Validation figures with all 4 panels

### 4. **Figure Panels Exactly as Requested:**
- **Panel A:** Genome-predicted substrate counts vs phenotype data (Spearman correlation)
- **Panel B:** CaCo-predicted niche overlap vs phenotypic measurements
- **Panel C:** Niche overlap distributions from genomes (specialist-specialist, generalist-generalist, mixed)
- **Panel D:** Niche overlap distributions from phenotypes (specialist-specialist, generalist-generalist, mixed)

### 5. **Robust Error Handling:**
- Falls back to mock data if CaCo fails
- Handles missing files gracefully
- Provides detailed error messages

### 6. **Comprehensive Output:**
- PNG and PDF versions of the figure
- Detailed summary report
- All intermediate data files
- Structured directory organization

## To Use:

1. **Save all files** in the same directory:
   - `caco_biolog_validation_pipeline.py` (main pipeline)
   - `run_validation_pipeline.py` (convenience runner)
   - `genome_list.txt` (your genome paths -- available in data/ folder)
   - `csource.tsv` (BIOLOG data)
   - `CaCo.py` (original CaCo script)

2. **Install dependencies:**

```
pip install pandas numpy matplotlib scipy seaborn tqdm
```

3. **Run the pipeline:**

```
python run_validation_pipeline.py
```
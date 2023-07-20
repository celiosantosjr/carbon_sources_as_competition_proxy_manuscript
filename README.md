# Carbon Sources as competition proxy manuscript

This folder contains the scripts needed for reproducing the work 'Carbon sources
as a proxy for microbial competition in the oceans'.

### Requirements

All scripts were written either in bash/shell or python3 language.
To run these scripts you will need the following packages:

- tqdm
- pandas
- numpy
- scipy
- geopandas
- macrel
- seaborn
- matplotlib

You will need to install in your path:

- prodigal
- mmseqs
- CaCo

### Info

We need to run the first and most crucial step of these scripts the program
[CaCo](github.com/celiosantosjr/CaCo) that we designed to calculate the carbon substrates competition 
scores. You can find more information on its github repository, just follow
the [link](github.com/celiosantosjr/CaCo).

### How to run

Scripts were written in their correct order, so run the scripts from 00 until 13th 
using the classical:

```
python3 <script>
```

There are some exceptions:

- After script `00_download_files.py` you should run the shell script it generates 
  before running the script `01_cleaning_files.py`

- The script `10_get_AMPs.sh` is a shell script and should be run in the Macrel
  installed environment.

### Distribution

In this folder it is also distributed an enlarged version of the previously
AMPRG database annotated by [Kintses et al. (2018)](https://pubmed.ncbi.nlm.nih.gov/30559406/).

### Organization by analysis

| **Figure/Analysis** | **Script** |
| :---: | :---: | 
| figure 1A-C | 02_graphs_fig1.py |
| figure 1D | 03_niche_overlap.py |
| figure 2A | 'fluxogram designed on Canvas.com' |
| figure 2B-C | 04_substrates_per_depth__ocean.py |
| figure 3 | 09_analysis_by_species.py |
| figure 4A | 08_test_most_competitive.py |
| figure 4B | 07_result_intertaxa.py |
| figure 5 | 06_competition_vs_longhurst.py |
| figure S1A-D | 02_graphs_fig1.py |
| figure S2 | 05_correlate_prevalence_in_MAGs.py |


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

df = pd.read_table('output/substrates_by_site.tsv')

# correlation with latitude
cors = []
for col in df.columns[6:]:
    r, p = spearmanr(df.latitude, df[col])
    cors.append(('latitude', col, r, p))

cors = pd.DataFrame(cors, columns=['x', 'y', 'spearman_rho', 'p-value'])
_, cors['padj'], _, _ = multipletests(cors['p-value'])
cors.to_csv('output/correlations_latitude.tsv', sep='\t', header=True, index=None)

# correlation with longitude
cors = []
for col in df.columns[6:]:
    r, p = spearmanr(df.longitude, df[col])
    cors.append(('longitude', col, r, p))

cors = pd.DataFrame(cors, columns=['x', 'y', 'spearman_rho', 'p-value'])
_, cors['padj'], _, _ = multipletests(cors['p-value'])
cors.to_csv('output/correlations_longitude.tsv', sep='\t', header=True, index=None)

# correlation with ocean
cors = []
for col in df.columns[6:]:
    r, p = spearmanr(df.ocean, df[col])
    cors.append(('ocean', col, r, p))

cors = pd.DataFrame(cors, columns=['x', 'y', 'spearman_rho', 'p-value'])
_, cors['padj'], _, _ = multipletests(cors['p-value'])
cors.to_csv('output/correlations_ocean.tsv', sep='\t', header=True, index=None)

# correlation with depth
cors = []
for col in df.columns[6:]:
    r, p = spearmanr(df.depthlayer, df[col])
    cors.append(('depthlayer', col, r, p))

cors = pd.DataFrame(cors, columns=['x', 'y', 'spearman_rho', 'p-value'])
_, cors['padj'], _, _ = multipletests(cors['p-value'])
cors.to_csv('output/correlations_depth_layer.tsv', sep='\t', header=True, index=None)

# plotting
# latitude
lat = pd.read_table('output/correlations_latitude.tsv')
lat = lat[lat.padj < 0.05]
for _, x, y, _, _, _ in lat.itertuples():
    sns.scatterplot(data=df, x=x, y=y, s=0.5, alpha=0.25, color='gray')
    sns.regplot(data=df, x=x, y=y, color='gray')
    plt.title(y)
    plt.ylabel('Prevalence in MAGs')
    plt.xlabel('Latitude')
    plt.tight_layout()
    plt.savefig(f'output/{x}__{y}.svg')
    plt.close()

# longitude
lon = pd.read_table('output/correlations_longitude.tsv')
lon = lon[lon.padj < 0.05]
for _, x, y, _, _, _ in lon.itertuples():
    sns.scatterplot(data=df, x=x, y=y, s=0.5, alpha=0.25, color='gray')
    sns.regplot(data=df, x=x, y=y, color='gray')
    plt.title(y)
    plt.ylabel('Prevalence in MAGs')
    plt.xlabel('Longitude')
    plt.tight_layout()
    plt.savefig(f'output/{x}__{y}.svg')
    plt.close()

# ocean
oc = pd.read_table('output/correlations_ocean.tsv')
oc = oc[oc.padj < 0.05]
for _, x, y, _, _, _ in oc.itertuples():
    sns.boxplot(data=df, x=x, y=y, showfliers=False, color='white', width=0.4)
    sns.swarmplot(data=df, x=x, y=y, color='gray')
    plt.title(y)
    plt.ylabel('Prevalence in MAGs')
    plt.xlabel('Ocean')
    plt.tight_layout()
    plt.savefig(f'output/{x}__{y}.svg')
    plt.close()

# depth
dl = pd.read_table('output/correlations_depth_layer.tsv')
dl = dl[dl.padj < 0.05]
for _, x, y, _, _, _ in dl.itertuples():
    sns.boxplot(data=df, x=x, y=y, showfliers=False, color='white', width=0.4)
    sns.swarmplot(data=df, x=x, y=y, color='gray')
    plt.title(y)
    plt.ylabel('Prevalence in MAGs')
    plt.xlabel('Depth Layer')
    plt.tight_layout()
    plt.savefig(f'output/{x}__{y}.svg')
    plt.close()
    

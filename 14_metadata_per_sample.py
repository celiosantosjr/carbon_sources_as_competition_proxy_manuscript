import pandas as pd
import seaborn as sns
import geopandas as gpd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from matplotlib import pyplot as plt

smeta_tara = pd.read_excel('sampling_sites_metadata.xlsx',
                          sheet_name='1.TARA_CONTEXT_95_SEQ_STATIONS_')


#meta_malaspina = pd.read_excel('sampling_sites_metadata.xlsx',
#                                sheet_name='MALASPINA_ATLANTECO_WP2_amplicon_samples_20210713')


genomes = pd.read_table('data/summarized_hq_genomes.tsv')
genomes['BioSamples'] = genomes.Genome.apply(lambda x: x.split('_')[1])

eit = pd.read_table('data/carboncomp_output.tsv.xz')
eit['s1'] = eit.genome1.apply(lambda x: x.split("_")[1])
eit['s2'] = eit.genome2.apply(lambda x: x.split("_")[1])
eit = eit[eit.s1 == eit.s2]

a = genomes.merge(on='BioSamples', right=meta_tara)

eit = eit.drop('s2', axis=1)
eit = eit.rename({'s1': 'BioSamples'}, axis=1)
eit = eit.groupby('BioSamples')['competition'].agg('mean')
eit = eit.reset_index()

eit = eit.merge(on='BioSamples', right=meta_tara)

eit = eit.drop(['Sample.ID', 'Sample.mat', 'Sample',
                'Station.label', 'Campaign.label',
                'Event.label', 'Event.device',
                'Event.comment', 'Event.date',
                'lower.size.fraction', 'upper.size.fraction',
                'Method', 'Marine.biome',
                'Ocean.region', 'Biogeographical.province',
                'Depth.nominal'], axis=1)

res = []
for col in eit.columns[2:]:
     r, p = spearmanr(eit['competition'], eit[col], nan_policy='omit')
     res.append(('competition', col, r, p))

res = pd.DataFrame(res, columns=['x', 'y', 'rho', 'p-value'])
_, res['padj'], _, _ = multipletests(res['p-value'])

print(res[res['padj'] < 0.05])

sns.scatterplot(data=eit, x='CO3', y='competition')
plt.ylabel('Average competition per sample')
plt.tight_layout()
plt.savefig("output/sample_competition.svg")
plt.close()

sns.scatterplot(data=eit, x='FC.picoeukaryotes.cells.mL', y='competition')
plt.ylabel('Average competition per sample')
plt.tight_layout()
plt.savefig("output/sample_competition_FCpicoeuk.svg")
plt.close()

fig, ax = plt.subplots()
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
world.plot(ax=ax, color='lightgray') 

sns.scatterplot(ax=ax, data=eit, x='Longitude', y='Latitude', hue='competition', s=5, palette='YlOrBr')
plt.savefig('output/sample_competition_global.svg')
plt.savefig('output/sample_competition_global.png', dpi=300)
plt.close()



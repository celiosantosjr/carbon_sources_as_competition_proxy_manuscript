import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from tqdm import tqdm
from scipy.stats import kruskal

print('Loading data') 
genomes = pd.read_table('data/summarized_hq_genomes.tsv')
genomes = genomes[['Genome','dRep Species Cluster', 'GTDB Taxonomy']]
genomes['GTDB Taxonomy'] = genomes['GTDB Taxonomy'].apply(lambda x: x.split(";")[1])
genomes = genomes.rename({'GTDB Taxonomy': 'phylum', 'Genome': 'genome'}, axis=1)

sps = genomes.set_index('genome')['dRep Species Cluster'].to_dict()
phyla = genomes[['dRep Species Cluster', 'phylum']].drop_duplicates().set_index('dRep Species Cluster')['phylum'].to_dict()

eit = pd.read_table('data/carboncomp_output.tsv.xz')
eit = eit[eit.prob < 0.05] 
eit = eit[eit.genome1.isin(genomes.genome) & eit.genome2.isin(genomes.genome)]
eit = eit.drop(['set1', 'set2', 'intersection',
                'relcomp', 'prob', 'relEIT'], axis=1)

print('Attributting species')
eit['dRep1'] = eit.genome1.apply(lambda x: sps.get(x))
eit['dRep2'] = eit.genome2.apply(lambda x: sps.get(x))
eit['ph1'] = eit.dRep1.apply(lambda x: phyla.get(str(x)))
eit['ph2'] = eit.dRep2.apply(lambda x: phyla.get(str(x))) 
eit = eit.dropna()

print('Processing genome species')
g = []
for i in tqdm(set(sps.values()), desc='Processing clusters'):
    x = eit[(eit.dRep1 == i) | (eit.dRep2 == i)]['competition']
    n1 = len(x)
    n2 = x.mean()
    n3 = x.std()
    g.append([i, n1, n2, n3])
    
g = pd.DataFrame(g, columns=['dRep', 'total interactions', 'average comp', 'std'])

g['error'] = 1.96 * np.sqrt((g['average comp'] * (1-g['average comp']) / g['total interactions']))

g['phylum'] = g.dRep.apply(lambda x: phyla.get(x))

g.to_csv('output/phylum_test_result.tsv', sep='\t', header=True, index=None)

g = g.dropna()
k = g['phylum'].value_counts()
k = k[k >= 10]
k = k.index
g = g[g.phylum.isin(k)]

print('Testing groups with Kruskal-Wallis')
ktest = {p: np.array(g[g.phylum == p]['average comp'].tolist()) for p in k if len(g[g.phylum == p]) >= 10}

print(kruskal(*list(ktest.values())))
#KruskalResult(statistic=152.61709530129886, pvalue=1.0783532366668832e-27)

print('Plotting differences')

# defining order
order = g.groupby('phylum').apply(lambda x: x['average comp'].quantile(0.5)).sort_values().index
sns.boxplot(data=g, x='phylum', y='average comp', showfliers=False, width=0.4, color='white', order=order)
sns.swarmplot(data=g, x='phylum', y='average comp', s=2, order=order)
plt.ylabel('Average competition')
plt.xlabel('')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('output/phylum_comp.svg')



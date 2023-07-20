import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests

print('Preparing oceans db')
data = pd.read_table('provinces_output.tsv', names=['latitude', 'longitude', 'longhurst', 'desc'])

print('Genomes data')
genomes = pd.read_table('data/summarized_hq_genomes.tsv')
genomes = genomes.merge(on=['latitude', 'longitude'], right=data)
provdict = genomes[['Genome', 'longhurst']].set_index('Genome').to_dict()['longhurst']

print('Loading EIT')
eit = pd.read_table('data/carboncomp_output.tsv.xz')
eit = eit[eit.prob < 0.05]
eit = eit[eit.genome1.isin(genomes.Genome) & eit.genome2.isin(genomes.Genome)]

eit['a'] = eit.genome1.apply(lambda x: provdict.get(x))
eit['b'] = eit.genome2.apply(lambda x: provdict.get(x))
eit = eit[eit['a'] == eit['b']]
eit = eit.rename({'a': 'longhurst'}, axis=1).drop('b', axis=1)

k = eit.longhurst.value_counts()
k = k[k>=10].index
eit = eit[eit.longhurst.isin(k)]

print('testing')
s, p = kruskal(*[eit[eit.longhurst == x]['competition'].tolist() for x in k])
print(f'Kruskal had a significant result (p={p:.2E})')
# Kruskal had a significant result (p=4.02E-70)

test = []
for k1, k2 in combinations(k, 2):
    _, p = mannwhitneyu(eit[eit.longhurst == k1]['competition'].tolist(), eit[eit.longhurst == k2]['competition'].tolist())
    test.append([k1, k2, p])

test = pd.DataFrame(test, columns=['x1', 'x2', 'pval'])
_, test['padj'], _, _ = multipletests(test['pval'])
test = test[test.padj < 0.05]
test.to_csv('output/province_comp_test.tsv', sep='\t', header=True, index=None)

print('plotting')
q50 = {x: eit[eit.longhurst == x]['competition'].quantile(0.5) for x in k}
q50 = pd.DataFrame.from_dict(q50, orient='index').sort_values(by=0).index
sns.boxplot(data=eit, y='competition', x='longhurst', order=q50, color='white', showfliers=False)
sns.stripplot(data=eit, x='longhurst', y='competition', order=q50, s=1.5, palette='Dark2')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('output/provinces.svg')
plt.savefig('output/provinces.png', dpi=300)


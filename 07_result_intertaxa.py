import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

print('Load data')
genomes = pd.read_table('data/summarized_hq_genomes.tsv')
eit = pd.read_table('data/carboncomp_output.tsv.xz')
eit = eit[eit.prob < 0.05]

print('Prepare taxonomy keys')
genomes = genomes[['Genome', 'GTDB Taxonomy']]
genomes['GTDB Taxonomy'] = genomes['GTDB Taxonomy'].apply(lambda x: x.split(';'))

genomes['p'] = genomes['GTDB Taxonomy'].apply(lambda x: x[1].replace('p__', '0'))
phylum = genomes[['Genome', 'p']].set_index('Genome').to_dict()['p']

genomes['c'] = genomes['GTDB Taxonomy'].apply(lambda x: x[2].replace('c__', '0'))
clas = genomes[['Genome', 'c']].set_index('Genome').to_dict()['c']

genomes['o'] = genomes['GTDB Taxonomy'].apply(lambda x: x[3].replace('o__', '0'))
order = genomes[['Genome', 'o']].set_index('Genome').to_dict()['o']

genomes['f'] = genomes['GTDB Taxonomy'].apply(lambda x: x[4].replace('f__', '0'))
family = genomes[['Genome', 'f']].set_index('Genome').to_dict()['f']

genomes['g'] = genomes['GTDB Taxonomy'].apply(lambda x: x[5].replace('g__', '0'))
genus = genomes[['Genome', 'g']].set_index('Genome').to_dict()['g']

print('Prepare phyla')
eit['p1'] = eit.genome1.apply(lambda x: phylum.get(x, '0'))
eit['p2'] = eit.genome2.apply(lambda x: phylum.get(x, '0'))

print('Prepare class')
eit['c1'] = eit.genome1.apply(lambda x: clas.get(x, '0'))
eit['c2'] = eit.genome2.apply(lambda x: clas.get(x, '0'))

print('Prepare order')
eit['o1'] = eit.genome1.apply(lambda x: order.get(x, '0'))
eit['o2'] = eit.genome2.apply(lambda x: order.get(x, '0'))

print('Prepare family')
eit['f1'] = eit.genome1.apply(lambda x: family.get(x, '0'))
eit['f2'] = eit.genome2.apply(lambda x: family.get(x, '0'))

print('Prepare genus')
eit['g1'] = eit.genome1.apply(lambda x: genus.get(x, '0'))
eit['g2'] = eit.genome2.apply(lambda x: genus.get(x, '0'))

print('Fix string replacement')
eit = eit.replace('0', np.nan)

# all results
print('Plot distributions')
def plot():
    fig, ([ax1, ax2, ax3, ax4, ax5]) = plt.subplots(5)
    sns.kdeplot(eit[eit.p1 == eit.p2]['EIT'], label='Intrataxa', ax=ax1)
    sns.kdeplot(eit[eit.p1 != eit.p2]['EIT'], label='Intertaxa', ax=ax1)
    sns.kdeplot(eit[eit.c1 == eit.c2]['EIT'], ax=ax2)
    sns.kdeplot(eit[eit.c1 != eit.c2]['EIT'], ax=ax2)
    sns.kdeplot(eit[eit.o1 == eit.o2]['EIT'], ax=ax3)
    sns.kdeplot(eit[eit.o1 != eit.o2]['EIT'], ax=ax3)
    sns.kdeplot(eit[eit.f1 == eit.f2]['EIT'], ax=ax4)
    sns.kdeplot(eit[eit.f1 != eit.f2]['EIT'], ax=ax4)
    sns.kdeplot(eit[eit.g1 == eit.g2]['EIT'], ax=ax5)
    sns.kdeplot(eit[eit.g1 != eit.g2]['EIT'], ax=ax5)
    labels = ['Phylum', 'Class', 'Order', 'Family', 'Genus']
    for idx, x in enumerate([ax1, ax2, ax3, ax4, ax5]):
        if idx < 4:
            x.set_xticks([])
            x.set_xlabel(np.nan)
        x.set_title(labels[idx])
    plt.subplots_adjust(right=0.85)
    fig.legend(loc='center right')
    plt.savefig('output/different_taxonomy_levels.svg')

plot()

print('Test statistics')
_, p = mannwhitneyu(eit[eit.p1 == eit.p2]['EIT'],
                    eit.loc[(eit.p1 != eit.p2), ['p1', 'p2', 'EIT']].dropna()['EIT'])

print(f'Phylum median two-sided difference of EIT: {p:.2E}')
#Phylum median two-sided difference of EIT: 0.00E+00

_, p = mannwhitneyu(eit[eit.c1 == eit.c2]['EIT'],
                    eit.loc[(eit.c1 != eit.c2), ['c1', 'c2', 'EIT']].dropna()['EIT'])

print(f'Class median two-sided difference of EIT: {p:.2E}')
#Class median two-sided difference of EIT: 1.22E-242

_, p = mannwhitneyu(eit[eit.o1 == eit.o2]['EIT'],
                    eit.loc[(eit.o1 != eit.o2), ['o1', 'o2', 'EIT']].dropna()['EIT'])

print(f'Order median two-sided difference of EIT: {p:.2E}')
#Order median two-sided difference of EIT: 2.08E-211

_, p = mannwhitneyu(eit[eit.f1 == eit.f2]['EIT'],
                    eit.loc[(eit.f1 != eit.f2), ['f1', 'f2', 'EIT']].dropna()['EIT'])

print(f'Family median two-sided difference of EIT: {p:.2E}')
#Family median two-sided difference of EIT: 1.32E-176

_, p = mannwhitneyu(eit[eit.g1 == eit.g2]['EIT'],
                    eit.loc[(eit.g1 != eit.g2), ['g1', 'g2', 'EIT']].dropna()['EIT'])

print(f'Genus median two-sided difference of EIT: {p:.2E}')
#Genus median two-sided difference of EIT: 1.55E-123


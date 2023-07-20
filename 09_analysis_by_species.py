import pandas as pd
from tqdm import tqdm
from scipy.stats import spearmanr
import seaborn as sns
import matplotlib.pyplot as plt

print('Loading data') 
genomes = pd.read_table('data/summarized_hq_genomes.tsv')
genomes = genomes[['Genome','dRep Species Cluster']]
genomes = genomes.set_index('Genome').to_dict()['dRep Species Cluster']

eit = pd.read_table('data/carboncomp_output.tsv.xz')
eit = eit[eit.prob < 0.05] 

print('Attributting species')
eit['p1'] = eit.genome1.apply(lambda x: genomes.get(x))
eit['p2'] = eit.genome2.apply(lambda x: genomes.get(x))
eit = eit.dropna()

print('Processing genome species')
g = []
for i in tqdm(set(genomes.values()), desc='Processing clusters'):
    x = eit[(eit.p1 == i) | (eit.p2 == i)]['EIT']
    n1 = len(x)
    n2 = sum(x >= 0)
    n3 = x.mean()
    g.append([i, n1, n2, n3])
    
g = pd.DataFrame(g,
                 columns=['species',
                          'total interactions',
                          'positive interactions',
                          'average EIT'])

g = g.sort_values(by=['total interactions',
                      'positive interactions',
                      'average EIT'],
                  ascending=[False, False, False])

g = g.dropna()
g.to_csv('output/results_per_drep_species.tsv.gz', sep='\t', header=True, index=None)

print('Testing correlation of total interactions and average EIT')
print(spearmanr(g['total interactions'], g['average EIT']))
# SignificanceResult(statistic=0.04695979308344871, pvalue=0.15694169533691085)

print('Testing correlation of total interactions and positive interactions')
print(spearmanr(g['total interactions'], g['positive interactions']))
# SignificanceResult(statistic=0.9298547209964025, pvalue=0.0)

sns.regplot(data=g, y='total interactions', x='positive interactions')
plt.xlabel('Total interactions per species')
plt.ylabel('Number of positive or zeroed EITs')
plt.savefig('output/Totalint_vs_positive.svg')
plt.close()

sns.jointplot(data=g, y='total interactions', x='average EIT')
plt.xlabel('Average EIT per species')
plt.ylabel('Total interactions per species')
plt.savefig('output/totalint_vs_EIT.svg')
plt.close()


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

print('Genomes data')
genomes = pd.read_table('data/summarized_hq_genomes.tsv')

print('selecting most abundant species')
k = genomes.groupby('dRep Species Cluster').agg('size').sort_values()
k = k[k >= 20].index

print('getting taxonomy key')
taxkey = genomes[genomes['dRep Species Cluster'].isin(['3085_1', '4138_1', '546_1', '6203_1'])]
taxkey = taxkey[['dRep Species Cluster', 'GTDB Taxonomy']].drop_duplicates()
taxkey = taxkey.set_index('dRep Species Cluster').to_dict()
taxkey = taxkey['GTDB Taxonomy']

# Manually editted to:
# {'6203_1': 'Amylibacter_A spp.', '3085_1': 'TCS55 sp002715035', '4138_1': 'UBA10364 spp.', '546_1': 'GCA-2725285 spp.'}

print('getting niche')
x = genomes[genomes['dRep Species Cluster'].isin(k)]
x = x[['dRep Species Cluster', 'temperature (°C)', 'oxygen (µmol/kg)']]
x = x.sort_values(by='dRep Species Cluster')
x['taxonomy'] = x['dRep Species Cluster'].apply(lambda x: taxkey.get(x))

print('plotting niche overlap')
sns.jointplot(data=x,
              x='temperature (°C)',
              y='oxygen (µmol/kg)',
              hue='taxonomy',
              kind='kde',
              palette='Dark2')

plt.tight_layout()
plt.savefig('output/niche_overlap.svg')
plt.savefig('output/niche_overlap.png', dpi=300)
plt.close()


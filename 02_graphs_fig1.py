import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import geopandas as gpd
from glob import glob
from Bio import SeqIO

genomes = pd.read_table('data/summarized_hq_genomes.tsv')
genomes['phylum'] = genomes['GTDB Taxonomy'].apply(lambda x: x.split(";")[1])

# plot fig1C
sns.histplot(genomes.groupby(['latitude', 'longitude', 'depth']).agg('size'))
plt.xlabel('MAGs per site')
plt.tight_layout()
plt.savefig('output/MAGs_per_site.svg')
plt.savefig('output/MAGs_per_site.png', dpi=300)
plt.close()

# plot fig1B
x = genomes.groupby(['depth layer', 'phylum']).agg('size').reset_index().pivot(index='phylum', columns='depth layer', values=0).fillna(0)
x = x.loc[x.sum(axis=1).sort_values().index]
x.loc['others'] = x.iloc[0:21].sum(axis=0)
x = x.iloc[21:] / x.iloc[21:].sum(axis=0)
x = x.loc[x.sum(axis=1).sort_values().index]

sns.heatmap(x)
plt.tight_layout()
plt.savefig('output/taxonomy.svg')
plt.savefig('output/taxonomy.png', dpi=300)
plt.close()

# plot fig1A
fig, ax = plt.subplots()
world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
world.plot(ax=ax, color='lightgray') 

sns.scatterplot(ax=ax, data=genomes, x='longitude', y='latitude', hue='depth layer', s=5, palette='Dark2')
plt.savefig('output/genome_distribution.svg')
plt.savefig('output/genome_distribution.png', dpi=300)
plt.close()


# plotting figure S1A-D
sns.scatterplot(data=genomes, x='CheckM Completeness', y='CheckM Contamination', alpha=0.5)
plt.tight_layout()
plt.savefig('output/S1A.svg')
plt.close()

prots = pd.read_table('protein_per_genome.tsv')
sns.histplot(prots)
plt.xlabel('Predicted proteins per MAG')
plt.tight_layout()
plt.savefig('output/S1D.svg')
plt.close()

gfolder = 'data/genomes-fasta/'

seqinfo = []
for fname in genomes.Genome:
    L, g, c, contigs = 0, 0, 0, 0
    for record in SeqIO.parse(f'{gfolder}/{fname}.fa', 'fasta'):
        contigs += 1
        g += str(record.seq).upper().count('G')
        c += str(record.seq).upper().count('C')
        L += len(record.seq)
    seqinfo.append((fname, L, (g+c)/L, contigs))
    print(fname, L, (g+c)/L, contigs)

seqinfo = pd.DataFrame(seqinfo, columns=['genome', 'length', 'GC', 'contigs'])
seqinfo.to_csv('output/genomes_contigs_length_gc.tsv', sep='\t', header=True, index=None)
seqinfo['length'] = seqinfo['length'] / 1e6

sns.scatterplot(data=seqinfo, x='length', y='GC', alpha=0.5)
plt.xlabel('Genome Length (Mbp)')
plt.ylabel('GC content')
plt.tight_layout()
plt.savefig('output/S1B.svg')
plt.savefig('output/S1B.png', dpi=300)
plt.close()

sns.histplot(data=seqinfo, x='contigs')
plt.xlabel("Contigs per MAG")
plt.tight_layout()
plt.savefig('output/S1C.svg')
plt.savefig('output/S1C.png', dpi=300)
plt.close()

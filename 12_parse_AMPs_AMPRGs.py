import json
import pandas as pd
from tqdm import tqdm
from glob import glob

print('Load table')
df = []
for infile in tqdm(glob('data/*_macrel.tsv.gz')):
    name = infile.split('/')[-1].replace('_macrel.tsv.gz', '')
    ndf = pd.read_table(infile, header='infer', comment='#')
    ndf['Genome'] = name
    df.append(ndf)

print('Form db')
df = pd.concat(df)

print('Export')
df.to_csv('output/macrel_all_output.tsv.gz',
          sep='\t',
          header=True,
          index=None)

print('Load dicts')
genes = json.load(open('AMPRG_database/access2genes.json'))
targets = json.load(open('AMPRG_database/genes2AMPs.json'))
mechanism = json.load(open('AMPRG_database/genes2mechanisms.json'))

print('Load table')
df = []
for infile in tqdm(glob("data/*_AMPRGs.m8")):
    name = infile.split("/")[-1].replace('_AMPRGs.m8', '')
    try:
        ndf = pd.read_table(infile, sep='\t', header=None)
        ndf['Genome'] = name
        df.append(ndf)
    except:
        pass  # in this case the file is empty

df = pd.concat(df)

print('Reduce db')
df = df.groupby(0)
df = df.apply(lambda x: x.sort_values(by=[11, 10, 2],
                                      ascending=[False, True, False]).head(1))

df = df.reset_index(drop=True)

df.columns=['query', 'target', 'identity', 'alnlen',
            'mismatch', 'gaps', 'qstart', 'qend',
            'tstart', 'tend', 'evalue', 'bitscore',
            'Genome']

print('Convert results')
df['gene'] = df['target'].apply(lambda x: genes.get(x))
df['amp_res'] = df.gene.apply(lambda x: targets.get(x))
df['mechanism'] = df.gene.apply(lambda x: mechanism.get(x))

print('Export')
df.to_csv('output/AMPRGs.tsv.gz',
          sep='\t',
          header=True,
          index=None)

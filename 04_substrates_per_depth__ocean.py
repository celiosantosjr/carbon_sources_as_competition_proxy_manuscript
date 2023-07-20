import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import chain
from collections import Counter

# import tables
df = pd.read_table('data/summarized_hq_genomes.tsv')
subs = pd.read_table('data/allsubs.tsv', sep='\t')

# merge info
subs.rename({'genome':'Genome'}, axis=1, inplace=True)
df = df.merge(on='Genome', right=subs)
df['substrates'] = df['substrates'].apply(lambda x: x.split(', ') if x == x else [''])

# select substrates per sample
x = df.groupby(['latitude',
                'longitude',
                'size fraction',
                'depth',
                'depth layer'])

x = x.apply(lambda x: dict(Counter(chain.from_iterable(x['substrates']))))

w = df.groupby(['latitude',
                'longitude',
                'size fraction',
                'depth',
                'depth layer']).agg('size')

x = pd.concat([x, w], axis=1).rename({0: 'subs', 1: 'size'}, axis=1)
x['threshold'] = x['size'] * 0.75
x['threshold'] = x.threshold.apply(lambda w: 2 if w == 1 else int(w))
x = x[x['size'] > 1]


# determine substrates per depth
x = x.reset_index()
layers = x['depth layer'].value_counts().to_dict()

epi, mes, bat, aby = [],[],[],[]
for idx, latitude, longitude, sizefraction, depth, depthlayer, subs, t0, t in x.itertuples():
    if t == 1: t = t0
    subs = [k for k, v in subs.items() if v >= t]
    if depthlayer == 'EPI': epi.append(subs)
    if depthlayer == 'MES': mes.append(subs)
    if depthlayer == 'BAT': bat.append(subs)
    if depthlayer == 'ABY': aby.append(subs)

epi = dict(Counter(chain.from_iterable(epi)))
mes = dict(Counter(chain.from_iterable(mes)))
bat = dict(Counter(chain.from_iterable(bat)))
aby = dict(Counter(chain.from_iterable(aby)))

epi = {k: v/layers.get('EPI') for k, v in epi.items()}
mes = {k: v/layers.get('MES') for k, v in mes.items()}
bat = {k: v/layers.get('BAT') for k, v in bat.items()}
aby = {k: v/layers.get('ABY') for k, v in aby.items()}

## joining keys
subs = [list(w.keys()) for w in [epi, bat, mes, aby]]
subs = set(chain.from_iterable(subs))
subs = {k for k in subs if len(k) > 0} # eliminating empty keys

print(f'It was found {len(subs)} substrates')

## preparing dataframe
L = len(subs)
ndf = np.zeros((4, L))

for idx, w in enumerate(subs):
    ndf[0, idx] = epi.get(w, 0)
    ndf[1, idx] = mes.get(w, 0)
    ndf[2, idx] = bat.get(w, 0)
    ndf[3, idx] = aby.get(w, 0)

ndf = pd.DataFrame(ndf, columns=subs, index=['EPI', 'MES', 'BAT', 'ABY']).T
ndf = ndf.drop('ABY', axis=1)

# plot substrates per depth
sns.clustermap(data=ndf, yticklabels=1)
plt.tight_layout()
plt.savefig('output/fig2B.svg')
plt.close()

# select substrates per ocean
long_oceans = {'ANTA': 'Antarctic',
               'APLR': 'Antarctic',
               'ARAB': 'Indian',
               'ARCH': 'Pacific', 
               'ARCT': 'Atlantic', 
               'AUSE': 'Pacific', 
               'AUSW': 'Indian', 
               'BENG': 'Atlantic', 
               'BPLR': 'Arctic', 
               'CAMR': 'Pacific', 
               'CARB': 'Atlantic',
               'CCAL': 'Pacific', 
               'CHIL': 'Pacific',  
               'CNRY': 'Atlantic',
               'EAFR': 'Indian', 
               'FKLD': 'Atlantic',
               'GFST': 'Atlantic',
               'GUIA': 'Atlantic',
               'ISSG': 'Indian', 
               'MEDI': 'Atlantic',
               'MONS': 'Indian', 
               'NADR': 'Atlantic',
               'NASE': 'Atlantic',
               'NASW': 'Atlantic',
               'NATR': 'Atlantic',
               'NPTG': 'Pacific',
               'NWCS': 'Atlantic',
               'PEQD': 'Pacific',
               'PNEC': 'Pacific',
               'REDS': 'Indian',
               'SANT': 'Antarctic',
               'SARC': 'Atlantic',
               'SATL': 'Atlantic',
               'SPSG': 'Pacific',
               'SSTC': 'Antarctic',
               'WTRA': 'Atlantic',
               }

long = pd.read_table('provinces_output.tsv', names=['latitude', 'longitude', 'longhurst', 'definition'])
long['oceans'] = long.longhurst.apply(lambda x: long_oceans.get(x))

ocean = x.merge(on=['latitude', 'longitude'], right=long).drop(['definition', 'longhurst'], axis=1)

# determine substrates per depth
antarctic, atlantic, arctic, indian, pacific = [], [], [], [], []
for idx, latitude, longitude, sizefraction, depth, depthlayer, subs, t0, t, sea in ocean.itertuples():
    if t == 1: t = t0
    subs = [k for k, v in subs.items() if v >= t]
    if sea == 'Antarctic': antarctic.append(subs)
    if sea == 'Atlantic': atlantic.append(subs)
    if sea == 'Arctic': arctic.append(subs)
    if sea == 'Indian': indian.append(subs)
    if sea == 'Pacific': pacific.append(subs)

antarctic = dict(Counter(chain.from_iterable(antarctic)))
atlantic = dict(Counter(chain.from_iterable(atlantic)))
arctic = dict(Counter(chain.from_iterable(arctic)))
indian = dict(Counter(chain.from_iterable(indian)))
pacific = dict(Counter(chain.from_iterable(pacific)))

layers = ocean['oceans'].value_counts().to_dict()

antarctic = {k: v/layers.get('Antarctic') for k, v in antarctic.items()}
atlantic = {k: v/layers.get('Atlantic') for k, v in atlantic.items()}
arctic = {k: v/layers.get('Arctic') for k, v in arctic.items()}
indian = {k: v/layers.get('Indian') for k, v in indian.items()}
pacific = {k: v/layers.get('Pacific') for k, v in pacific.items()}

## joining keys
subs = [list(w.keys()) for w in [antarctic, atlantic, arctic, indian, pacific]]
subs = set(chain.from_iterable(subs))
subs = {k for k in subs if len(k) > 0} # eliminating empty keys

print(f'It was found {len(subs)} substrates')

## preparing dataframe
L = len(subs)
ndf = np.zeros((5, L))

for idx, w in enumerate(subs):
    ndf[0, idx] = antarctic.get(w, 0)
    ndf[1, idx] = atlantic.get(w, 0)
    ndf[2, idx] = arctic.get(w, 0)
    ndf[3, idx] = indian.get(w, 0)
    ndf[4, idx] = pacific.get(w, 0)

ndf = pd.DataFrame(ndf, columns=subs, index=['Antarctic', 'Atlantic', 'Arctic', 'Indian', 'Pacific']).T

# plot substrates per depth
sns.clustermap(data=ndf, yticklabels=1)
plt.tight_layout()
plt.savefig('output/fig2C.svg')
plt.close()

# create table of substrates per site
keys = list(set(chain.from_iterable([k.keys() for k in ocean.subs])))[1:]
ndf = []
for idx, latitude, longitude, sizefraction, depth, depthlayer, subs, t0, t, sea in ocean.itertuples():
    sdf = [latitude, longitude, sizefraction, depth, depthlayer, sea]
    for substrate in keys:
        sdf.append(subs.get(substrate, 0) / t0)
    ndf.append(sdf)  

cols = ['latitude', 'longitude', 'sizefraction', 'depth', 'depthlayer', 'ocean'] + keys  
ndf = pd.DataFrame(ndf, columns=cols)
ndf.to_csv('output/substrates_by_site.tsv', sep='\t', header=True, index=None)


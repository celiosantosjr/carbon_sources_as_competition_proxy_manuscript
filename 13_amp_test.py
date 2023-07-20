import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

print('Load data')
hq = pd.read_table('data/summarized_hq_genomes.tsv')
hq = set(hq.Genome)

print('Load data')
rg = pd.read_table('output/AMPRGs.tsv.gz')
rg = rg[rg.Genome.isin(hq)]
rg = set(rg.Genome)

print('Load data')
amps = pd.read_table('output/macrel_all_output.tsv.gz')
amps = amps[amps.Genome.isin(hq)]
amps = set(amps.Genome)

print('Load data')
subs = pd.read_table('data/carboncomp_output.tsv.xz')
subs = subs[subs.genome1.isin(hq) & subs.genome2.isin(hq) & (subs.prob < 0.05)]

print('Generate none set')
nothing = hq - rg - amps

print('Subsetting')
amp_eit = subs[subs.genome1.isin(amps) & subs.genome2.isin(amps)]['EIT']
rg_eit = subs[subs.genome1.isin(rg) & subs.genome2.isin(rg)]['EIT']
amp_rg_eit = subs[(subs.genome1.isin(amps) & subs.genome2.isin(rg)) | (subs.genome1.isin(rg) & subs.genome2.isin(amps))]['EIT']
nothing_amp_eit = subs[(subs.genome1.isin(amps) & subs.genome2.isin(nothing)) | (subs.genome1.isin(nothing) & subs.genome2.isin(amps))]['EIT']
nothing_rg_eit = subs[(subs.genome1.isin(rg) & subs.genome2.isin(nothing)) | (subs.genome1.isin(nothing) & subs.genome2.isin(rg))]['EIT']

print('Preparing structure')
all_eit = pd.DataFrame([], columns=['All', 'AMPvsAMP', 'AMPRGvsAMPRG', 'AMPvsAMPRG', 'AMPvsNone', 'AMPRGvsNone'])
all_eit['All'] = subs['EIT']
all_eit['AMPvsAMP'] = amp_eit
all_eit['AMPRGvsAMPRG'] = rg_eit
all_eit['AMPvsAMPRG'] = amp_rg_eit
all_eit['AMPvsNone'] = nothing_amp_eit
all_eit['AMPRGvsNone'] = nothing_rg_eit
print(all_eit.head())

print('Testing variables')
test = []
for x, y in combinations(all_eit.columns, 2):
    _, p = mannwhitneyu(all_eit[x].dropna(), all_eit[y].dropna())
    test.append((x, y, p))

test = pd.DataFrame(test, columns=['group1', 'group2', 'p-value'])
test.to_csv("output/AMP_AMPRGs_groups_comparison.tsv", sep='\t', header=True, index=None)
 
print('Plotting')
sns.boxplot(data=all_eit, showfliers=False, width=0.4, color='white', order=['All', 'AMPvsAMP', 'AMPvsAMPRG', 'AMPvsNone', 'AMPRGvsAMPRG', 'AMPRGvsNone'])
sns.swarmplot(data=all_eit.sample(2000), s=1.5, order=['All', 'AMPvsAMP', 'AMPvsAMPRG', 'AMPvsNone', 'AMPRGvsAMPRG', 'AMPRGvsNone'])
plt.ylabel('Ecological Interaction Type (EIT)')
plt.xlabel('Group tested')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('output/AMPeffect2.svg')
plt.savefig('output/AMPeffect2.png', dpi=300)
plt.show()


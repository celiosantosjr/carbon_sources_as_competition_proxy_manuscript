import os
import pandas as pd
from subprocess import run
from glob import glob

print('Creating directories')
os.makedirs('data', exist_ok=True)
os.makedirs('output', exist_ok=True)

print('Downloading resource files')
genomes = 'https://sunagawalab.ethz.ch/share/microbiomics/ocean/suppl_data/genomes-fasta.tar.gz'
metadata = 'https://sunagawalab.ethz.ch/share/microbiomics/ocean/suppl_data/genomes-summary.csv.gz'

run(['wget', '-O', 'data/genomes-fasta.tar.gz', genomes])
run(['wget', '-O', 'data/genomes-summary.csv.gz', metadata])

print('Filtering data')
metadata = pd.read_table('data/genomes-summary.csv.gz', sep=',')
metadata['Quality'] = metadata['CheckM Completeness'] - 5*metadata['CheckM Contamination']

f1 = metadata['CheckM Completeness'] >= 90
f2 = metadata['CheckM Contamination'] <= 5
f3 = metadata['Quality'] >= 50
f4 = (metadata['# scaffolds'] <= 500) & (metadata['N50 (scaffolds)'] >= 10000)
f5 = metadata['# tRNAs (AA)'] >= 15

metadata = metadata[f1 & f2 & f3 & f4 & f5]

metadata = metadata[~metadata['size fraction'].fillna('-').str.contains('<')]

f6 = metadata.groupby(['latitude', 'longitude', 'depth layer', 'size fraction']).agg('size')
f6 = f6[f6>1].index

metadata = pd.concat([metadata[(metadata.latitude == x[0]) & (metadata.longitude == x[1]) & (metadata['depth layer'] == x[2]) & (metadata['size fraction'] == x[3])] for x in f6])

metadata = metadata.reset_index(drop=True)

metadata = metadata[['Genome', 'CheckM Completeness', 'CheckM Contamination',
                     '# tRNAs (AA)', '16S rRNA', 'GTDB Taxonomy',
                     'dRep Species Cluster', 'station', 'depth', 'depth layer',
                     'latitude', 'longitude', 'size fraction', 'temperature (°C)',
                     'oxygen (µmol/kg)']]

metadata.to_csv('data/summarized_hq_genomes.tsv', sep='\t', header=True, index=None)

run(['rm', '-rf', 'data/genomes-summary.csv.gz'])

print('Uncompressing genomes')
run(['tar', '-zxvf', 'data/genomes-fasta.tar.gz', '-C', 'data/'])
run(['rm', '-rf', 'data/genomes-fasta.tar.gz'])

for infile in glob("data/genomes-fasta/*"):
    sample = infile.split("/")[-1].replace('.fa', '')
    if sample in metadata.Genome.tolist():
        pass
    else:
        run(['rm', '-rf', infile])
        
print('Creating lists for CaCo')
glist = metadata.groupby(['latitude', 'longitude', 'depth']).apply(lambda x: x.Genome.tolist())

for idx, gs in enumerate(glist):
    with open(f'data/caco_gl_batch_{idx}.txt', 'w') as ofile:
        for g in gs:
            ofile.write(f'genomes-fasta/{g}.fa\n')

print('Creating CaCo script')
with open('data/calling_CaCo.sh', 'w') as ofile:
    for infile in glob('data/caco_gl_batch*.txt'):
        batch = infile.split("/")[-1].replace('caco_gl_', '').replace('.txt', '')
        infile = infile.split('/')[-1]
        ofile.write(f'python3 CaCo.py -m from_nucleotides -gl {infile}\n')
        ofile.write(f'mv allfams.tsv allfams_{batch}.tsv\n')
        ofile.write(f'mv allsubs.tsv allsubs_{batch}.tsv\n')
        ofile.write(f'mv carboncomp_output.tsv.xz carboncomp_output_{batch}.tsv.xz\n')
        ofile.write(f'rm -rf tmp/\n\n')
    ofile.write(f'cat allfams_batch*.tsv > allfams.tsv\n')
    ofile.write(f'rm -rf allfams_batch*.tsv\n\n')
    ofile.write(f'cat allsubs_batch*.tsv > allsubs.tsv\n')
    ofile.write(f'rm -rf allsubs_batch*.tsv\n\n')
    ofile.write(f'xzcat carboncomp_output_*.tsv.xz | xz > carboncomp_output.tsv.xz\n')
    ofile.write(f'rm -rf carboncomp_output_*.tsv.xz\n')
    
# now you have all batch files and the script to run CaCo
# after running the script calling_CaCo.sh you should have
# the files for the dbCAN families, substrates and competition/EIT tables.
    

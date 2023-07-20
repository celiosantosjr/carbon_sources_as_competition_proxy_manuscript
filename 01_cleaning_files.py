import os
import pandas as pd
from subprocess import run
from glob import glob

print('Cleaning GLs')
for infile in glob('data/caco_gl_batch_*.txt'):
    run(['rm', '-rf', infile])
    
run(['rm', '-rf', 'data/calling_CaCo.sh'])

print('Fixing tables')

df = pd.read_table('data/carboncomp_output.tsv.xz', sep='\t', header='infer')
df = df[df.genome1 != 'genome1']
df.to_csv('data/carboncomp_output.tsv.xz', sep='\t', header=True, index=None)

df = pd.read_table('data/allfams.tsv', sep='\t', header='infer')
df = df[df.genome != 'genome']
df.to_csv('data/allfams.tsv', sep='\t', header=True, index=None)

df = pd.read_table('data/allsubs.tsv', sep='\t', header='infer')
df = df[df.genome != 'genome']
df.to_csv('data/allsubs.tsv', sep='\t', header=True, index=None)


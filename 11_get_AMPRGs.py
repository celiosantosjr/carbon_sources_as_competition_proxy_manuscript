from glob import glob
from subprocess import run

for i in glob('data/genomes-fasta/*.fa'):
    sample = i.split("/")[-1].replace('.fa', '')
    # run gene prediction
    run(['prodigal',
         '-a', 'tmp.faa',
         '-c',
         '-i', i,
         '-m',
         '-p', 'single'])
    # search AMPRGs
    run(['mmseqs',
         'easy-search',
         '-e', '1e-10',
         '-c', '0.6',
         '--cov-mode', '2',
         '--min-seq-id', '0.7',
         'tmp.faa',
         'AMPRG_database/AMPRGs_db.faa',
         f'data/{sample}_AMPRGs.m8', 'tmp/'])
    # clean
    run(['rm', '-rf', 'tmp/'])
    run(['rm', '-rf', 'tmp.faa'])
    

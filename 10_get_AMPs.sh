for i in $(ls data/genomes-fasta/);
    do
        macrel contigs --fasta data/genomes-fasta/$i --out outtest
        mv outtest/macrel.out.prediction.gz data/${i/.fa/_macrel.tsv.gz}
        rm -rvf outtest
    done

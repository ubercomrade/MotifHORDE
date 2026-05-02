uv run motifhorde \
 ./PEAKS035260_GATA2_O09100_MACS2.fa \
 ./PEAKS035260_GATA2_O09100_MACS2_gb.fa \
/home/anton/Documents/genomes/mm/mm10.ensembl.promoters.fa \
 output/ -t bamm -l 8-20-4 -c overlap --cj-metric co --cj-perm 100 --cj-filter score --cj-threshold 0.8 -v -b peaks --fpr 0.5 -m pauROC

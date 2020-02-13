#!/bin/bash -xe

### NOTE: make sure that bowtie (v1) and genmap are up-to-date and installed in PATH.

# download an ecoli reference and unzip
curl -o - ftp://ftp.ensemblgenomes.org/pub/bacteria/release-46/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.toplevel.fa.gz | gunzip > ecoli.fa
# create a bowtie1 index
bowtie-build ecoli.fa ecoli.bowtie1
# create a fasta file containing a k-mer for every position in the genomic (ordered)
awk '!/^>/ { seq = seq $0 } END { for (i = 1; i <= length(seq) - 20 + 1; ++i) { print ">" i; print substr(seq, i, 20) } }' ecoli.fa > ecoli.kmers.fa
# build index for genmap
genmap index -F ecoli.fa -I ecoli.gm.index
# use genmap to compute mappability
genmap map -I ecoli.gm.index -O . -K 20 -E 1 -fl -bg
# search kmers with bowtie1 (in its all-mapping mode)
bowtie -f -n 1 -l 20 -a -S -p 24 --sam-nohead ecoli.bowtie1 ecoli.kmers.fa > ecoli.sam
# sort sam and convert to bed file
cut -f1 ecoli.sam | sort -k1,1 -n | uniq -c | awk -v OFS='\t' 'BEGIN{ startPos=1 } { if (lastPos && $1 != lastOccCount) { print "Chromosome", (startPos - 1), ($2 - 1), lastOccCount; startPos=$2 } lastOccCount=$1; lastPos=$2 } END { print "Chromosome", (startPos - 1), lastPos, lastOccCount }' > ecoli.bowtie.bedgraph
# compare both mappability files
diff -s ecoli.bowtie.bedgraph ecoli.genmap.bedgraph

rm -r ecoli.bowtie.bedgraph ecoli.genmap.bedgraph ecoli.bowtie1* ecoli.fa ecoli.gm.index ecoli.sam ecoli.kmers.fa

#!/bin/bash

    clear

    date


## ~~ Clustering at desired threshold~~ ##

#For some analysis, it can be interesting to cluster all
#the sequences at a desired threshold. This script takes a fasta file
#and generates two documents: A centroid file with all the fastas of the threshold clusters,
#and a correspondence .tsv to be able to know the ASV/OTU relationships.

# Usearch path (where does your happy usearch live?)
usearch=scripts/programs/usearch

#Parameters!
id_job="16s_blanes_otu2asv99"
input='data/abund-tax-raw/asv_Blanes_16S.fasta' # the one with only 16S Blanes seqs of interestB
output=data/abund-tax-raw


# This is an example, any threshold can be applied
$usearch -cluster_smallmem ${input} \
         -id 0.99 \
         -centroids ${output}/otus_99.fa \
         -sortedby size \
         -uc ${output}/${id_job}_uc_clusters.uc


# A small script to create the tsv!
python3 scripts/preprocessing/clustering-otu/3_asv2otu.py ${input} \
                     ${output}/${id_job}_uc_clusters.uc > $output/${id_job}_correspondence.tsv

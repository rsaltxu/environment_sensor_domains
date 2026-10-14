#!/usr/bin/env R

#R script to test correlation between cophenetic distances of GTDB representatives and their jaccard distances
# on sensor domain composition

library(vegan)


###Load precomputed matrices
gtdb_coph=read.csv('/home/roman.sanz/data/diamond_pblast/gtdb_reps/GTDB_famrep_coph_dist_v2.csv',
 header = TRUE,row.names = 1)
jacc_HK=read.csv('/home/roman.sanz/data/diamond_pblast/gtdb_reps/GTDB_famrep_HK_jacc_dist.csv',
 header = TRUE,row.names = 1)

gtdb_coph=as.matrix(gtdb_coph)
jacc_HK=as.matrix(jacc_HK)

coph_rows=rownames(gtdb_coph)

jacc_rows=rownames(jacc_HK)

#Make rows and columns go in the same order
jacc_HK <- jacc_HK[coph_rows, coph_rows]

jacc_HK[is.na(jacc_HK)] <- 1

#Perform test
mantel(gtdb_coph,jacc_HK,method='spearman')
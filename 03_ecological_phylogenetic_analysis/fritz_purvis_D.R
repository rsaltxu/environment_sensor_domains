#!/usr/bin/env R

#R script to test the phylogenetic signal on ligand binding domain absence/presence

library(caper)
library(ape)
library(dplyr)
library(readr)


#Load data and keep one representative per family on the top 11 phyla
df <- read.csv("/home/roman.sanz/data/diamond_pblast/gtdb_reps/family_rep_HK_abundance_matrix_prot_90cov_50score.csv",
                 header = TRUE, row.names = 1)
keep_rows=read.table("/home/roman.sanz/data/diamond_pblast/gtdb_reps/family_reps_top_phyla_list.txt", 
                        header = FALSE, 
                        stringsAsFactors = FALSE)$V1
df=df[(rownames(df) %in% keep_rows),,drop=FALSE]

#load annotations for sensor domains-ligand
hk_annotation=read.csv('/home/roman.sanz/data/HK_related/outer_hk_fam_ligand_annot_v2.csv',header=TRUE, row.names = 1)

#load GTDB tree
tree=read.tree(file="/home/roman.sanz/data/phylomodel/family_reps_top_phyla_gtdb.tree")

###Process data to binary table
trait_groups <- list(
  Aminoacids = rownames(hk_annotation[hk_annotation$Ligand=='Aminoacids',]),
  Hydroxy = rownames(hk_annotation[hk_annotation$Ligand=='4-hydroxycinnamic acid',]),
  Biogenic_amine = rownames(hk_annotation[hk_annotation$Ligand=='Biogenic amine',]),
  FAD = rownames(hk_annotation[hk_annotation$Ligand=='FAD',]),
  Dival_cat = rownames(hk_annotation[hk_annotation$Ligand=='Divalent cation',]),
  b_type_heme = rownames(hk_annotation[hk_annotation$Ligand=='b-type heme binding',]),
  c_type_heme = rownames(hk_annotation[hk_annotation$Ligand=='c-type heme binding',]),
  
  Carboxylate = rownames(hk_annotation[hk_annotation$Ligand=='Carboxilate',])
)

binary_traits <- data.frame(Genome = rownames(df))

for(trait in names(trait_groups)) {
    print(trait)
    cols <- trait_groups[[trait]]
    print(cols)
    binary_traits[[trait]] <-as.integer(rowSums(df[, (colnames(df) %in% cols)] > 0, na.rm = TRUE) > 0)
}


#Load comparative object
GTDB_ligands <- comparative.data(tree,binary_traits, 'Genome')


###Tests for every ligand

res <- phylo.d(GTDB_ligands, binvar = Aminoacids)
print(res)
print('---')
res <- phylo.d(GTDB_ligands, binvar = Hydroxy)
print(res)
print('---')
res <- phylo.d(GTDB_ligands, binvar = Biogenic_amine)
print(res)
print('---')
res <- phylo.d(GTDB_ligands, binvar = Dival_cat)
print(res)
print('---')
res <- phylo.d(GTDB_ligands, binvar = b_type_heme)
print(res)
print('---')
res <- phylo.d(GTDB_ligands, binvar = FAD)
print(res)
print('---')
res <- phylo.d(GTDB_ligands, binvar = Carboxylate)
print(res)
print('---')
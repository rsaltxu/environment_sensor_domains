#!/usr/bin/env R

library(ggplot2)
library(ComplexUpset)
library(readr)

df <- read.csv("/home/roman.sanz/data/diamond_pblast/gtdb_reps/HK_abspres_matrix_prot_90cov_50score.csv",
                 header = TRUE, row.names = 1)

keep_cols <- read.table("/home/roman.sanz/data/tmhmm/outer_HKs.txt", 
                        header = FALSE, 
                        stringsAsFactors = FALSE)$V1

# Filter your data matrix to keep only these columns
df <- df[, (colnames(df) %in% keep_cols), drop = FALSE]                 
phyla = rownames(df)

phyla = c('Bacteroidota','Chloroflexota','Bacillota','Gammaproteobacteria','Alphaproteobacteria','Actinomycetota',
'Verrucomicrobiota','Acidobacteriota','Desulfobacterota','Planctomycetota','Cyanobacteriota')
df_bool <- df[, -ncol(df)]
df_bool <- t(df_bool)
df_bool<- df_bool == 1

sets=colnames(df_bool)

df_bool <- as.data.frame(df_bool) 

png("/home/roman.sanz/data/figures/upset_plot_gtdb_HK_nosetsize_v2_outer_hks.png", 
    width = 2100, height = 1600, res = 300)

upset(
  df_bool,
  intersect = phyla,
  base_annotations=list(
      'Number of shared sensory domains'=intersection_size(counts=FALSE)
  ),
  name = "Phyla",
  width_ratio = 0.1,
    min_size=80,
    min_degree=1, set_sizes=FALSE)+ theme(
    axis.text.y = element_text(size = 12))
dev.off()


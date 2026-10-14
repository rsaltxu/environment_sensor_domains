#!/usr/bin/env R

library(grid)
library(pheatmap)

hk_to_keep= read.table("/home/roman.sanz/data/tmhmm/outer_HKs.txt", 
                        header = FALSE, 
                        stringsAsFactors = FALSE)$V1#we must filter to just get external HKS

original=read.csv("/home/roman.sanz/data/singlem/australia/australia_custom_abundances.csv",header=TRUE,
row.names=1)

signif_corr2 = read.csv("/home/roman.sanz/data/australia_significant_correlations_singlem.csv",header = TRUE, row.names = 1)
all_corr = read.csv("/home/roman.sanz/data/australia_all_correlations_singlem.csv",header = TRUE, row.names = 1)

signif_corr2 <- signif_corr2[(rownames(signif_corr2) %in% hk_to_keep), , drop = FALSE]# in this case, phyla are pre-filtered

all_corr <- all_corr[(rownames(all_corr) %in% hk_to_keep), , drop = FALSE]

##CHANGE LABELS
signif_corr2 =as.matrix(signif_corr2)
all_corr=as.matrix(all_corr)

signif_corr2 <- signif_corr2[rowSums(!is.na(signif_corr2)) > 1,
                             colSums(!is.na(signif_corr2)) > 1]
mat <- signif_corr2
mat[is.na(mat)] <- 0

row_dist <- dist(mat, method = "euclidean")
col_dist <- dist(t(mat), method = "euclidean")

Aridity_prevalence <- read.csv("/home/roman.sanz/data/singlem/australia/australia_custom_WAridity.csv", header=TRUE, row.names = 1)

Aridity_prevalence=Aridity_prevalence[c('W')]

png("/home/roman.sanz/data/figures/spercorr_r_singlem_zoom_v2_W.png", 
    width = 2500, height = 2500, res = 300)
pheatmap(mat,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_method = "complete",
            annotation_col = Aridity_prevalence,
        annotation_colors = list(W = colorRampPalette(c('#666666',"#fecc5c","#e31a1c"))(50)),
         na_col = "grey90",
        show_rownames = FALSE,
        fontsize_col=15,
        angle_col = 45,
         clustering_distance_rows =row_dist,
         clustering_distance_cols = col_dist
  )
dev.off()
# ----------------------
# Write column names to separate files
writeLines(cluster1_rows, "/home/roman.sanz/data/HK_related/spear_cor_arid_1201.txt")
writeLines(cluster2_rows, "/home/roman.sanz/data/HK_related/spear_cor_nonarid_1213.txt")

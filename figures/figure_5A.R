#!/usr/bin/env R

library(pheatmap)
library(gridExtra)
library(grid)
library(vegan)

# --- Load data ---
data <- read.csv("/home/roman.sanz/data/diamond_pblast/australia/HK_abundance_matrix_prot_50cov_50score.csv",
                 header = TRUE, row.names = 1)

# Read the list of allowed column names
keep_cols <- read.table("/home/roman.sanz/data/tmhmm/outer_HKs.txt", 
                        header = FALSE, 
                        stringsAsFactors = FALSE)$V1

# Filter your data matrix to keep only these columns
data <- data[, (colnames(data) %in% keep_cols), drop = FALSE]

data_matrix <- as.matrix(data)
data_matrix <- data_matrix[, apply(data_matrix, 2, function(x) any(x != 0))]

threshold <- 10   # for example: require at least 2 non-zero entries

nonzero_counts <- colSums(data_matrix != 0)
# Keep columns where non-zero occurrences meet the threshold
data_matrix <- data_matrix[, nonzero_counts >= threshold]
cat("Number of columns:", ncol(data_matrix), "\n")
data_matrix <- log2(1e6*(data_matrix) + 1)


# ---------------------------
s_metadata <- read.csv("/home/roman.sanz/scripts/data/mapfile_base_subset_Aus.csv",
                       header = TRUE, row.names = 1)

metagen_metadata <- read.csv("/home/roman.sanz/scripts/data/Metadata_metagenomes_selected_255_Aus.csv",
                       header = TRUE, row.names = 2)
taxprof_data <- read.csv("/home/roman.sanz/data/singlem/australia/taxprof_singlem.csv",header=TRUE, row.names = 1)

HK_annotation <- read.csv("/home/roman.sanz/data/HK_related/outer_hk_fam_ligand_annot_v2.csv", header=TRUE, row.names = 1)

HK_annotation <- HK_annotation[colnames(data_matrix), ]
HK_annotation <- HK_annotation[order(HK_annotation$Family), ]

data_matrix <- data_matrix[,rownames(HK_annotation)]


clean_names <- gsub("^X", "", rownames(s_metadata))
rownames(s_metadata) <- clean_names

# --- Keep only the single metric column ---
metric_col <- "Aridity.index...annual.mean"

# Subset to only rows present in the heatmap
aridity_vector <- s_metadata[rownames(data_matrix), "Aridity.index...annual.mean"]
aridity_vector <- cut(
  aridity_vector,
  breaks = c(-Inf, 0.2, 0.65, Inf),
  labels = c("Arid", "Semiarid", "Humid"),
  right = TRUE,
  include.lowest = TRUE
)
aridity_vector <- s_metadata[rownames(data_matrix), "Aridity.index...annual.mean"]

names(aridity_vector) <- rownames(data_matrix)
npp_vector <- s_metadata[rownames(data_matrix), "NPP.Mean"]
names(npp_vector) <- rownames(data_matrix)
conductivity_vector <- metagen_metadata[rownames(data_matrix), "conductivity"]
names(conductivity_vector) <- rownames(data_matrix)
general_env <- metagen_metadata[rownames(data_matrix), "general_env_feature"]
names(general_env) <- rownames(data_matrix)
carbon <- metagen_metadata[rownames(data_matrix), "organic_carbon"]
names(carbon) <- rownames(data_matrix)
ph <- metagen_metadata[rownames(data_matrix), "ph"]
names(ph) <- rownames(data_matrix)
gammaprot <- taxprof_data[rownames(data_matrix),"Gammaproteobacteria"]
names(gammaprot) <- rownames(data_matrix)

####COLUMN ANNOTATION

HK_domain <- HK_annotation[colnames(data_matrix),"Family"]
names(HK_domain) <- colnames(data_matrix)

names(HK_domain) <- colnames(data_matrix)
Ligand <- HK_annotation[colnames(data_matrix),"Ligand"]
names(Ligand) <- colnames(data_matrix)

# ---------------------------
annotation_df <- data.frame(Aridity = aridity_vector)
annotation_df <- data.frame(
  Aridity = aridity_vector,
  Habitat = general_env,
  pH = ph,
  NPP= npp_vector,
  Conductivity = conductivity_vector,
  Organic_C=carbon
)

rownames(annotation_df) <- names(general_env)

annotation_cols <- data.frame(HK_Family=HK_domain,
  Ligand=Ligand)
rownames(annotation_cols) <- colnames(data_matrix)

# Replace some name entries
annotation_df$Habitat[annotation_df$Habitat == "Coastal (coastal terrestrial including cliff; rocky shore)"] <- 
  "Coastal"
annotation_df$Habitat[annotation_df$Habitat == "Riverine (rivers and streams)"] <- 
  "Riverine"
  annotation_df$Habitat[annotation_df$Habitat == ""] <- 
  "Other"

annotation_df$Habitat <- factor(annotation_df$Habitat)

annotation_cols$HK_Family <- factor(annotation_cols$HK_Family)
annotation_cols$Ligand <- factor(annotation_cols$Ligand)

# ---------------------------
# Combine barplot + heatmap
# ---------------------------
png("/home/roman.sanz/data/figures/CPM_Austr_outer_HKs_ligand_known.png", 
    width = 2000, height = 1600, res = 300)+
pheatmap(
  data_matrix,
  annotation_row = annotation_df,
  annotation_col = annotation_cols,
  annotation_colors = list(
  Aridity = c(
    "Arid" = "#E31A1C",
    "Semiarid" = "#F89540",
    "Humid" = "#666666"
  ),
  Habitat = c(
    "Arid" = "#d95f02",  
    "Coastal" = "#7570b3", 
    "Montane" = "#e7298a",
    "Riverine" = "#1f78b4",
    "Temperate" = "#1b9e77",
    "Tropical (dry)" = "#e6ab02",
    "Tropical (wet)" = "#66a61e",
    "Other" = "white"
  ),
  HK_Family = c(
    "dCache" = "#d95f02", 
    "PAS" = "#1f78b4"
  ),
  Ligand = c(
    "4-hydroxycinnamic acid" = "#1b9e77",
    "Aminoacids" = "#d95f02",
    "Biogenic amine" = "#cab600",
    "b-type heme binding" = "#9e0202",
    "Carboxilate" = "#e7298a",
    "Divalent cation" ="#66a61e",
    "FAD" = "#1f78b4"
  )
  ),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = FALSE,
  show_colnames=FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  fontsize_col = 6
)
dev.off()
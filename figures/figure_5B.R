#!/usr/bin/env R

library(ggplot2)

ligands=c('AminoAc','Purine','dival_cat', 
'btype_heme','ctype_heme','FAD','Carboxilate', 
'4HydroCinn','BioAm')

### Now we follow the same procedure for every ligand type
### BioAm

df <- read.csv('/data/HK_related/BioAm_CPM_Aus.csv',header=TRUE,row.names=1)
df$Habitat <- factor(df$Habitat,levels = c("Arid", "Semiarid", "Humid"))
#print(df)

# Shapiro-Wilk normality test for each group
#by(df$BioAm_CPM, df$Habitat, shapiro.test)

# Define pairwise comparisons
comparisons <- list(
  c("Arid", "Semiarid"),
  c("Semiarid", "Humid"),
  c("Arid", "Humid")
)

# If not normal → Kruskal-Wallis
kruskal.test(BioAm_CPM ~ Habitat, data = df)

# Pairwise comparisons
#pairwise.wilcox.test(df$BioAm_CPM, df$Habitat,
#                     p.adjust.method = "bonferroni")
# Base boxplot with ggplot2
p <- ggplot(df, aes(x = Habitat, y = BioAm_CPM, fill = Habitat)) +
    geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(
    aes(colour = Habitat),
    width = 0.15,
    size = 2,
    alpha = 0.6
  ) +
  scale_color_manual(values = c(
    "Arid" = "#bb0508",
    "Semiarid" = "#d2701a",
    "Humid" = "#4c4b4b"
  ))+
  scale_fill_manual(values = c(
    "Arid" = "#E31A1C",
    "Semiarid" = "#F89540",
    "Humid" = "#666666"
  ))+
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18),  # center title
    #axis.title.x = element_text(size = 22),             # x-axis label
    #axis.title.y = element_text(size = 22),             # y-axis label
    axis.text.x = element_text(size = 20,angle=45,hjust=1,vjust=1),              # x tick labels
    axis.text.y = element_text(size = 20)               # y tick labels
    ) +
  labs(title = "",
       x = "",
       y = "")

# Add significance bars manually
p 
ggsave("/data/figures/Aus_BioAm_CPM_aridity_boxplot_annotfree.png",
       width = 4,
       height = 5,
       units = "in",
       dpi = 300)


### AminoAc

df <- read.csv('/data/HK_related/AminoAc_CPM_Aus.csv',header=TRUE,row.names=1)
df$Habitat <- factor(df$Habitat,levels = c("Arid", "Semiarid", "Humid"))

# Shapiro-Wilk normality test for each group
by(df$AminoAc_CPM, df$Habitat, function(x) {
  if(length(unique(x)) < 2) {
    return("all values identical")
  }
  
  if(length(x) < 3) {
    return("too few observations")
  }
  
  shapiro.test(x)
})

# If not normal → Kruskal-Wallis
kruskal.test(AminoAc_CPM ~ Habitat, data = df)

# Pairwise comparisons
#pairwise.wilcox.test(df$AminoAc_CPM, df$Habitat,
#                     p.adjust.method = "bonferroni")

# Base boxplot with ggplot2
p <- ggplot(df, aes(x = Habitat, y = AminoAc_CPM, fill = Habitat)) +
    geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(
    aes(colour = Habitat),
    width = 0.15,
    size = 2,
    alpha = 0.6
  ) +
  scale_color_manual(values = c(
    "Arid" = "#bb0508",
    "Semiarid" = "#d2701a",
    "Humid" = "#4c4b4b"
  ))+
  scale_fill_manual(values = c(
    "Arid" = "#E31A1C",
    "Semiarid" = "#F89540",
    "Humid" = "#666666"
  ))+
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.text.x = element_text(size = 20,angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size = 20)
    ) +
  labs(title = "",
       x = "",
       y = "")

# Add significance bars manually
p

ggsave("/data/figures/Aus_AminoAc_CPM_aridity_boxplot_annotfree.png",
       width = 4,
       height = 5,
       units = "in",
       dpi = 300)


### Purine

df <- read.csv('/data/HK_related/Purine_CPM_Aus.csv',header=TRUE,row.names=1)
df$Habitat <- factor(df$Habitat,levels = c("Arid", "Semiarid", "Humid"))

# Shapiro-Wilk normality test for each group
by(df$Purine_CPM, df$Habitat, function(x) {
  if(length(unique(x)) < 2) {
    return("all values identical")
  }
  
  if(length(x) < 3) {
    return("too few observations")
  }
  
  shapiro.test(x)
})

# If not normal → Kruskal-Wallis
kruskal.test(Purine_CPM ~ Habitat, data = df)

# Pairwise comparisons
#pairwise.wilcox.test(df$Purine_CPM, df$Habitat,
#                     p.adjust.method = "bonferroni")

# Base boxplot with ggplot2
p <- ggplot(df, aes(x = Habitat, y = Purine_CPM, fill = Habitat)) +
    geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(
    aes(colour = Habitat),
    width = 0.15,
    size = 2,
    alpha = 0.6
  ) +
  scale_color_manual(values = c(
    "Arid" = "#bb0508",
    "Semiarid" = "#d2701a",
    "Humid" = "#4c4b4b"
  ))+
  scale_fill_manual(values = c(
    "Arid" = "#E31A1C",
    "Semiarid" = "#F89540",
    "Humid" = "#666666"
  ))+
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.text.x = element_text(size = 20,angle=45,hjust=1,vjust=1), 
    axis.text.y = element_text(size = 20)
    ) +
  labs(title = "",
       x = "",
       y = "")

# Add significance bars manually
p 

ggsave("/data/figures/Aus_Purine_CPM_aridity_boxplot_annotfree.png",
       width = 4,
       height = 5,
       units = "in",
       dpi = 300)


### dival_cat

df <- read.csv('/data/HK_related/dival_cat_CPM_Aus.csv',header=TRUE,row.names=1)
df$Habitat <- factor(df$Habitat,levels = c("Arid", "Semiarid", "Humid"))

# Shapiro-Wilk normality test for each group
by(df$dival_cat_CPM, df$Habitat, function(x) {
  if(length(unique(x)) < 2) {
    return("all values identical")
  }
  
  if(length(x) < 3) {
    return("too few observations")
  }
  
  shapiro.test(x)
})

# If not normal → Kruskal-Wallis
kruskal.test(dival_cat_CPM ~ Habitat, data = df)

# Pairwise comparisons
#pairwise.wilcox.test(df$dival_cat_CPM, df$Habitat,
#                     p.adjust.method = "bonferroni")

# Base boxplot with ggplot2
p <- ggplot(df, aes(x = Habitat, y = dival_cat_CPM, fill = Habitat)) +
    geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(
    aes(colour = Habitat),
    width = 0.15,
    size = 2,
    alpha = 0.6
  ) +
  scale_color_manual(values = c(
    "Arid" = "#bb0508",
    "Semiarid" = "#d2701a",
    "Humid" = "#4c4b4b"
  ))+
  scale_fill_manual(values = c(
    "Arid" = "#E31A1C",
    "Semiarid" = "#F89540",
    "Humid" = "#666666"
  ))+
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18),  # center title
    axis.title.x = element_text(size = 22),             # x-axis label
    axis.title.y = element_text(size = 22),             # y-axis label
    axis.text.x = element_text(size = 20,angle=45,hjust=1,vjust=1),              # x tick labels
    axis.text.y = element_text(size = 20)               # y tick labels
  ) +
  labs(title = "",
       x = "",
       y = "")

# Add significance bars manually
p

ggsave("/data/figures/Aus_dival_cat_CPM_aridity_boxplot_annotfree.png",
       width = 4,
       height = 5,
       units = "in",
       dpi = 300)


### btype_heme

df <- read.csv('/data/HK_related/btype_heme_CPM_Aus.csv',header=TRUE,row.names=1)
df$Habitat <- factor(df$Habitat,levels = c("Arid", "Semiarid", "Humid"))

# Shapiro-Wilk normality test for each group
by(df$btype_heme_CPM, df$Habitat, function(x) {
  if(length(unique(x)) < 2) {
    return("all values identical")
  }
  
  if(length(x) < 3) {
    return("too few observations")
  }
  
  shapiro.test(x)
})

# If not normal → Kruskal-Wallis
kruskal.test(btype_heme_CPM ~ Habitat, data = df)

# Pairwise comparisons
#pairwise.wilcox.test(df$btype_heme_CPM, df$Habitat,
#                     p.adjust.method = "bonferroni")

# Base boxplot with ggplot2
p <- ggplot(df, aes(x = Habitat, y = btype_heme_CPM, fill = Habitat)) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(
    aes(colour = Habitat),
    width = 0.15,
    size = 2,
    alpha = 0.6
  ) +
  scale_color_manual(values = c(
    "Arid" = "#bb0508",
    "Semiarid" = "#d2701a",
    "Humid" = "#4c4b4b"
  ))+
  scale_fill_manual(values = c(
    "Arid" = "#E31A1C",
    "Semiarid" = "#F89540",
    "Humid" = "#666666"
  ))+
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18),  # center title
    axis.title.x = element_text(size = 22),             # x-axis label
    axis.title.y = element_text(size = 22),             # y-axis label
    axis.text.x = element_text(size = 20,angle=45,hjust=1,vjust=1),              # x tick labels
    axis.text.y = element_text(size = 20)               # y tick labels
    ) +
  labs(title = "",
       x = "",
       y = "")

# Add significance bars manually
p 

ggsave("/data/figures/Aus_btype_heme_CPM_aridity_boxplot_annotfree.png",
       width = 4,
       height = 5,
       units = "in",
       dpi = 300)


### ctype_heme

df <- read.csv('/data/HK_related/ctype_heme_CPM_Aus.csv',header=TRUE,row.names=1)
df$Habitat <- factor(df$Habitat,levels = c("Arid", "Semiarid", "Humid"))

# Shapiro-Wilk normality test for each group
by(df$ctype_heme_CPM, df$Habitat, function(x) {
  if(length(unique(x)) < 2) {
    return("all values identical")
  }
  
  if(length(x) < 3) {
    return("too few observations")
  }
  
  shapiro.test(x)
})

# If not normal → Kruskal-Wallis
kruskal.test(ctype_heme_CPM ~ Habitat, data = df)

# Pairwise comparisons
#pairwise.wilcox.test(df$ctype_heme_CPM, df$Habitat,
#                     p.adjust.method = "bonferroni")

# Base boxplot with ggplot2
p <- ggplot(df, aes(x = Habitat, y = ctype_heme_CPM, fill = Habitat)) +
    geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(
    aes(colour = Habitat),
    width = 0.15,
    size = 2,
    alpha = 0.6
  ) +
  scale_color_manual(values = c(
    "Arid" = "#bb0508",
    "Semiarid" = "#d2701a",
    "Humid" = "#4c4b4b"
  ))+
  scale_fill_manual(values = c(
    "Arid" = "#E31A1C",
    "Semiarid" = "#F89540",
    "Humid" = "#666666"
  ))+
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18),  # center title
    axis.title.x = element_text(size = 22),             # x-axis label
    axis.title.y = element_text(size = 22),             # y-axis label
    axis.text.x = element_text(size = 20,angle=45,hjust=1,vjust=1),              # x tick labels
    axis.text.y = element_text(size = 20)               # y tick labels
    ) +
  labs(title = "",
       x = "",
       y = "")

# Add significance bars manually
p

ggsave("/data/figures/Aus_ctype_heme_CPM_aridity_boxplot_annotfree.png",
       width = 4,
       height = 5,
       units = "in",
       dpi = 300)


### FAD

df <- read.csv('/data/HK_related/FAD_CPM_Aus.csv',header=TRUE,row.names=1)
df$Habitat <- factor(df$Habitat,levels = c("Arid", "Semiarid", "Humid"))

# Shapiro-Wilk normality test for each group
by(df$FAD_CPM, df$Habitat, function(x) {
  if(length(unique(x)) < 2) {
    return("all values identical")
  }
  
  if(length(x) < 3) {
    return("too few observations")
  }
  
  shapiro.test(x)
})

# If not normal → Kruskal-Wallis
#kruskal.test(FAD_CPM ~ Habitat, data = df)

# Pairwise comparisons
#pairwise.wilcox.test(df$FAD_CPM, df$Habitat,
#                     p.adjust.method = "bonferroni")

# Base boxplot with ggplot2
p <- ggplot(df, aes(x = Habitat, y = FAD_CPM, fill = Habitat)) +
    geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(
    aes(colour = Habitat),
    width = 0.15,
    size = 2,
    alpha = 0.6
  ) +
  scale_color_manual(values = c(
    "Arid" = "#bb0508",
    "Semiarid" = "#d2701a",
    "Humid" = "#4c4b4b"
  ))+
    scale_fill_manual(values = c(
    "Arid" = "#E31A1C",
    "Semiarid" = "#F89540",
    "Humid" = "#666666"
  ))+
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18),  # center title
    axis.title.x = element_text(size = 22),             # x-axis label
    axis.title.y = element_text(size = 22),             # y-axis label
    axis.text.x = element_text(size = 20,angle=45,hjust=1,vjust=1),              # x tick labels
    axis.text.y = element_text(size = 20)               # y tick labels
    ) +
  labs(title = "",
       x = "",
       y = "")

# Add significance bars manually
p 

ggsave("/data/figures/Aus_FAD_CPM_aridity_boxplot_annotfree.png",
       width = 4,
       height = 5,
       units = "in",
       dpi = 300)


### Carboxilate

df <- read.csv('/data/HK_related/Carboxilate_CPM_Aus.csv',header=TRUE,row.names=1)
df$Habitat <- factor(df$Habitat,levels = c("Arid", "Semiarid", "Humid"))
#print(df)

# Shapiro-Wilk normality test for each group
by(df$Carboxilate_CPM, df$Habitat, function(x) {
  if(length(unique(x)) < 2) {
    return("all values identical")
  }
  
  if(length(x) < 3) {
    return("too few observations")
  }
  
  shapiro.test(x)
})

# If not normal → Kruskal-Wallis
kruskal.test(Carboxilate_CPM ~ Habitat, data = df)

# Pairwise comparisons
pairwise.wilcox.test(df$Carboxilate_CPM, df$Habitat,
                     p.adjust.method = "bonferroni")

# Base boxplot with ggplot2
p <- ggplot(df, aes(x = Habitat, y = Carboxilate_CPM, fill = Habitat)) +
    geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(
    aes(colour = Habitat),
    width = 0.15,
    size = 2,
    alpha = 0.6
  ) +
  scale_color_manual(values = c(
    "Arid" = "#bb0508",
    "Semiarid" = "#d2701a",
    "Humid" = "#4c4b4b"
  ))+
  scale_fill_manual(values = c(
    "Arid" = "#E31A1C",
    "Semiarid" = "#F89540",
    "Humid" = "#666666"
  ))+
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18),  # center title
    axis.title.x = element_text(size = 22),             # x-axis label
    axis.title.y = element_text(size = 22),             # y-axis label
    axis.text.x = element_text(size = 20,angle=45,hjust=1,vjust=1),              # x tick labels
    axis.text.y = element_text(size = 20)               # y tick labels
    ) +
  labs(title = "",
       x = "",
       y = "")

# Add significance bars manually
p 

ggsave("/data/figures/Aus_Carboxilate_CPM_aridity_boxplot_annotfree.png",
       width = 4,
       height = 5,
       units = "in",
       dpi = 300)

### 4HydroCinn

df <- read.csv('/data/HK_related/HydroCinn_CPM_Aus.csv',header=TRUE,row.names=1)
df$Habitat <- factor(df$Habitat,levels = c("Arid", "Semiarid", "Humid"))

# Shapiro-Wilk normality test for each group
by(df$HydroCinn_CPM, df$Habitat, function(x) {
  if(length(unique(x)) < 2) {
    return("all values identical")
  }
  
  if(length(x) < 3) {
    return("too few observations")
  }
  
  shapiro.test(x)
})

# If not normal → Kruskal-Wallis
kruskal.test(HydroCinn_CPM ~ Habitat, data = df)

# Pairwise comparisons
#pairwise.wilcox.test(df$HydroCinn_CPM, df$Habitat,
#                     p.adjust.method = "bonferroni")

# Base boxplot with ggplot2
p <- ggplot(df, aes(x = Habitat, y = HydroCinn_CPM, fill = Habitat)) +
    geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(
    aes(colour = Habitat),
    width = 0.15,
    size = 2,
    alpha = 0.6
  ) +
  scale_color_manual(values = c(
    "Arid" = "#bb0508",
    "Semiarid" = "#d2701a",
    "Humid" = "#4c4b4b"
  ))+
  scale_fill_manual(values = c(
    "Arid" = "#E31A1C",
    "Semiarid" = "#F89540",
    "Humid" = "#666666"
  ))+
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.text.x = element_text(size = 20,angle=45,hjust=1,vjust=1),
    axis.text.y = element_text(size = 20)
    ) +
  labs(title = "",
       x = "",
       y = "")

# Add significance bars manually
p 

ggsave("/data/figures/Aus_HydroCinn_CPM_aridity_boxplot_annotfree.png",
       width = 4,
       height = 5,
       units = "in",
       dpi = 300)
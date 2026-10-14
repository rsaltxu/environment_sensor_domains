#!/usr/bin/env R

library(ggplot2)
library(vegan)

data  <- read.csv('/data/Australia_HK_BC_dissimilarity.csv',header = TRUE, #pre-computed matrix
  row.names = 1)
data_matrix <- as.matrix(data)

metadata <-read.csv('/data/Australia_merged_metadata_wphylum.csv',
                header = TRUE, row.names = 1)

metadata$general_env_feature[metadata$general_env_feature != "Arid"] <- 
  "Humid"

# dbRDA
keep <- complete.cases(metadata)
metadata <- metadata[keep, ]

###use continuous data
metadata$general_env_feature[metadata$Aridity.index...annual.mean > 0.65] <- 
  "Humid"
metadata$general_env_feature[
  metadata$Aridity.index...annual.mean > 0.2 &
  metadata$Aridity.index...annual.mean < 0.65
] <- "Semiarid"
metadata$general_env_feature[metadata$Aridity.index...annual.mean < 0.2] <- 
  "Arid"

data_matrix <- data_matrix[keep, ]

db <- capscale(data_matrix ~ NPP.Mean+conductivity+organic_carbon+ph, data = metadata)
# significance
anova(db)
anova(db, by = "terms")

groups <- metadata$general_env_feature
cols <- groups

sites <- as.data.frame(scores(db, display="sites", scaling=2))
bp <- as.data.frame(scores(db, display="bp", scaling=2))
sites$Environment <- metadata$general_env_feature

eig <- eigenvals(db)
var <- round(100 * eig / sum(eig), 1)

ggplot() +
  geom_point(data=sites, aes(CAP1, CAP2, color=Environment), size=5,alpha=0.55) +
  geom_segment(data=bp, aes(x=0,y=0,xend=CAP1,yend=CAP2),
               arrow=arrow(length=unit(0.2,"cm")), linewidth=0.6) +
  #geom_text(data=bp, aes(CAP1*2, CAP2*3, label=rownames(bp)), size=6) +
  geom_text(data=bp, aes(-1.2, 0, label='pH'), size=8) +
  geom_text(data=bp, aes(-0.9, 1.3, label='Conductivity'), size=8) +
  geom_text(data=bp, aes(1.4, 1, label='Organic C'), size=8) +
  geom_text(data=bp, aes(1, 0, label='NPP'), size=8) +
  coord_equal() +

  stat_ellipse(data = sites,
               aes(x = CAP1, y = CAP2, color = Environment),
               level = 0.95,
               linewidth = 0.8)+
  scale_color_manual(
    values = c(
    "Arid" = "#E31A1C",
    "Semiarid" = "#F89540",
    "Humid" = "#666666"
    )
  ) +
  theme_bw(base_size=18) +
  labs(x=paste0("dbRDA1 (",var[1],"%)"),
       y=paste0("dbRDA2 (",var[2],"%)"),
       #title="Distance-based RDA")
  )+
     theme(aspect.ratio = 1,
     axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26))
ggsave("/data/figures/australia_dbRDA_abiotic.png",
       width = 12,
       height = 8,
       units = "in",
       dpi = 300)
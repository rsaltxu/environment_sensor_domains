rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpmisc)
library(forcats)
library(stringr)
library(data.table)
library(terra)
library(tidyterra)
library(mgcv)
library(effsize)
library(colorspace)
library(ggpubr)
library(ggsignif)
library(tidyverse)
library(reshape2)

## Explore genomic features
trait.dbb = read.csv("C:/CEAB/9_Histidine_kinases/Analysis/Global_GTDB_ligand_trait_database.csv",
                            header = T)
names(trait.dbb)[names(trait.dbb) %in% c("b.type.heme.binding",
                                         "FAD",
                                         "Aminoacids",
                                         "X4.hydroxycinnamic.acid",
                                         "Biogenic.amine",
                                         "Divalent.cation")] = c("b-type heme", "FAD","Aminoacids",
                                                                 "4-hydroxycinnamic acid",
                                                                 "Biogenic amine", "Divalent cation")
trait.dbb$Lifestyle <- ifelse(
  trait.dbb$Lifestyle == "Environmental",
  "Free-living (environmental)",
  trait.dbb$Lifestyle
)
trait.dbb$Lifestyle <- ifelse(
  trait.dbb$Lifestyle == "Engineered",
  "Free-living (engineered)",
  trait.dbb$Lifestyle
)
head(trait.dbb)
nrow(trait.dbb) #31635 unique species with trait information

unique(trait.dbb$Habitat)
unique(trait.dbb$Oxygen)
unique(trait.dbb$pH)
unique(trait.dbb$Salinity)
unique(trait.dbb$Lifestyle)
unique(trait.dbb$Generalism)

# Heatmap
ligands <- c("FAD","Aminoacids",
  "4-hydroxycinnamic acid",
  "Biogenic amine", "Divalent cation"
)

habitat.hm <- reshape2::melt(
  trait.dbb,
  id.vars = "Habitat",
  measure.vars = ligands,
  variable.name = "Ligand",
  value.name = "Presence"
)

overall <- aggregate((Presence > 0) ~ Ligand, data = habitat.hm, FUN = mean)
names(overall)[2] <- "Overall"

habitat.hm <- aggregate(
  (Presence > 0) ~ Habitat + Ligand,
  data = habitat.hm,
  FUN = mean
)

names(habitat.hm)[3] <- "Prevalence"

habitat.hm <- merge(habitat.hm, overall, by = "Ligand")

habitat.hm$Habitat <- factor(
  habitat.hm$Habitat,
  levels = c(
    "Animal-associated",
    "Food/Fermentation",
    "Engineered",
    "Freshwater",
    "Marine",
    "Plant-associated",
    "Soil/Sediment",
    "Extreme",
    "Unknown"
  )
)

habitat.hm <- subset(
  habitat.hm,
  !is.na(Habitat) &
    Habitat != "Unknown"
)

habitat.hm$zscore <- ave(
  habitat.hm$Prevalence,
  habitat.hm$Ligand,
  FUN = function(x) (x - mean(x)) / sd(x)
)

ggplot(habitat.hm, aes(Habitat, Ligand, fill = zscore)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(
    low = "#BEBD7F",
    mid = "white",
    high = "#592321",
    midpoint = 0
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Stacked barplots
ligands <- c(
  "FAD",
  "Aminoacids",
  "4-hydroxycinnamic acid",
  "Biogenic amine",
  "Divalent cation"
)

plot.df <- reshape2::melt(
  trait.dbb,
  id.vars = "Habitat",
  measure.vars = ligands,
  variable.name = "Ligand",
  value.name = "Presence"
) |>
  dplyr::filter(!is.na(Habitat)) |>
  dplyr::group_by(Habitat, Ligand) |>
  dplyr::summarise(Prevalence = mean(Presence > 0), .groups = "drop")

plot.df = subset(plot.df, !Habitat == "Unknown")

ggplot(plot.df, aes(Habitat, Prevalence, fill = Ligand)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(plot.df, aes(Habitat, Prevalence, fill = Ligand)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Lifestyle
lifestyle.hm <- reshape2::melt(
  trait.dbb,
  id.vars = "Lifestyle",
  measure.vars = ligands,
  variable.name = "Ligand",
  value.name = "Presence"
)

overall <- aggregate((Presence > 0) ~ Ligand, data = lifestyle.hm, FUN = mean)
names(overall)[2] <- "Overall"

lifestyle.hm <- aggregate(
  (Presence > 0) ~ Lifestyle + Ligand,
  data = lifestyle.hm,
  FUN = mean
)

names(lifestyle.hm)[3] <- "Prevalence"

lifestyle.hm <- merge(lifestyle.hm, overall, by = "Ligand")

lifestyle.hm$Lifestyle <- factor(
  lifestyle.hm$Lifestyle,
  levels = c(
    "Host-associated",
    "Food-associated",
    "Free-living (engineered)",
    "Free-living (environmental)",
    "Unknown"
  )
)

lifestyle.hm <- subset(
  lifestyle.hm,
  !is.na(Lifestyle) &
    Lifestyle != "Unknown"
)

lifestyle.hm$zscore <- ave(
  lifestyle.hm$Prevalence,
  lifestyle.hm$Ligand,
  FUN = function(x) (x - mean(x)) / sd(x)
)

ggplot(lifestyle.hm, aes(Lifestyle, Ligand, fill = zscore)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(
    low = "#BEBD7F",
    mid = "white",
    high = "#592321",
    midpoint = 0
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot.df <- reshape2::melt(
  trait.dbb,
  id.vars = "Lifestyle",
  measure.vars = ligands,
  variable.name = "Ligand",
  value.name = "Presence"
) |>
  dplyr::filter(!is.na(Lifestyle)) |>
  dplyr::group_by(Lifestyle, Ligand) |>
  dplyr::summarise(
    Prevalence = mean(Presence > 0),
    .groups = "drop"
  )

plot.df = subset(plot.df, !Lifestyle == "Unknown")

ggplot(plot.df, aes(Lifestyle, Prevalence, fill = Ligand)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(plot.df, aes(Lifestyle, Prevalence, fill = Ligand)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Generalism (Yes/No)
generalism.hm <- reshape2::melt(
  trait.dbb,
  id.vars = "Generalism",
  measure.vars = ligands,
  variable.name = "Ligand",
  value.name = "Presence"
)

overall <- aggregate((Presence > 0) ~ Ligand, data = generalism.hm, FUN = mean)
names(overall)[2] <- "Overall"

generalism.hm <- aggregate(
  (Presence > 0) ~ Generalism + Ligand,
  data = generalism.hm,
  FUN = mean
)

names(generalism.hm)[3] <- "Prevalence"

generalism.hm <- merge(generalism.hm, overall, by = "Ligand")

generalism.hm$Generalism <- factor(
  generalism.hm$Generalism,
  levels = c("no", "yes")
)

generalism.hm <- subset(
  generalism.hm,
  !is.na(Generalism)
)

generalism.hm$zscore <- ave(
  generalism.hm$Prevalence,
  generalism.hm$Ligand,
  FUN = function(x) (x - mean(x)) / sd(x)
)
generalism.hm$Generalism <- ifelse(
  generalism.hm$Generalism == "yes",
  "Generalist",
  "Non-generalist"
)

ggplot(generalism.hm, aes(Generalism, Ligand, fill = zscore)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(
    low = "#BEBD7F",
    mid = "white",
    high = "#592321",
    midpoint = 0
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#Oxygen
oxygen.hm <- reshape2::melt(
  trait.dbb,
  id.vars = "Oxygen",
  measure.vars = ligands,
  variable.name = "Ligand",
  value.name = "Presence"
)

overall <- aggregate((Presence > 0) ~ Ligand, data = oxygen.hm, FUN = mean)
names(overall)[2] <- "Overall"

oxygen.hm <- aggregate(
  (Presence > 0) ~ Oxygen + Ligand,
  data = oxygen.hm,
  FUN = mean
)

names(oxygen.hm)[3] <- "Prevalence"

oxygen.hm <- merge(oxygen.hm, overall, by = "Ligand")

oxygen.hm$Oxygen <- factor(
  oxygen.hm$Oxygen,
  levels = c(
    "aerobic",
    "microaerophilic",
    "facultative anaerobic",
    "anaerobic",
    "aerotolerant",
    "No robust majority"
  )
)

oxygen.hm <- subset(
  oxygen.hm,
  !is.na(Oxygen) &
    Oxygen != "No robust majority" &
    Oxygen != "aerotolerant"
)

oxygen.hm$zscore <- ave(
  oxygen.hm$Prevalence,
  oxygen.hm$Ligand,
  FUN = function(x) (x - mean(x)) / sd(x)
)

ggplot(oxygen.hm, aes(Oxygen, Ligand, fill = zscore)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(
    low = "#BEBD7F",
    mid = "white",
    high = "#592321",
    midpoint = 0
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

plot.df <- reshape2::melt(
  trait.dbb,
  id.vars = "Oxygen",
  measure.vars = ligands,
  variable.name = "Ligand",
  value.name = "Presence"
) |>
  dplyr::filter(!is.na(Oxygen)) |>
  dplyr::group_by(Oxygen, Ligand) |>
  dplyr::summarise(
    Prevalence = mean(Presence > 0),
    .groups = "drop"
  )

plot.df = subset(plot.df, !Oxygen == "No robust majority")
plot.df = subset(plot.df, !Oxygen == "aerotolerant")

ggplot(plot.df, aes(Oxygen, Prevalence, fill = Ligand)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(plot.df, aes(Oxygen, Prevalence, fill = Ligand)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


#Temperature
temperature.hm <- reshape2::melt(
  trait.dbb,
  id.vars = "Temperature",
  measure.vars = ligands,
  variable.name = "Ligand",
  value.name = "Presence"
)

overall <- aggregate((Presence > 0) ~ Ligand, data = temperature.hm, FUN = mean)
names(overall)[2] <- "Overall"

temperature.hm <- aggregate(
  (Presence > 0) ~ Temperature + Ligand,
  data = temperature.hm,
  FUN = mean
)

names(temperature.hm)[3] <- "Prevalence"

temperature.hm <- merge(temperature.hm, overall, by = "Ligand")

temperature.hm$Temperature <- factor(
  temperature.hm$Temperature,
  levels = c(
    "psychrophile",
    "mesophile",
    "thermophile",
    "No robust majority"
  )
)

temperature.hm <- subset(
  temperature.hm,
  !is.na(Temperature) &
    Temperature != "No robust majority"
)

temperature.hm$zscore <- ave(
  temperature.hm$Prevalence,
  temperature.hm$Ligand,
  FUN = function(x) (x - mean(x)) / sd(x)
)

ggplot(temperature.hm, aes(Temperature, Ligand, fill = zscore)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(
    low = "#BEBD7F",
    mid = "white",
    high = "#592321",
    midpoint = 0
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))


# Salinity
salinity.hm <- reshape2::melt(
  trait.dbb,
  id.vars = "Salinity",
  measure.vars = ligands,
  variable.name = "Ligand",
  value.name = "Presence"
)

overall <- aggregate((Presence > 0) ~ Ligand, data = salinity.hm, FUN = mean)
names(overall)[2] <- "Overall"

salinity.hm <- aggregate(
  (Presence > 0) ~ Salinity + Ligand,
  data = salinity.hm,
  FUN = mean
)

names(salinity.hm)[3] <- "Prevalence"

salinity.hm <- merge(salinity.hm, overall, by = "Ligand")

salinity.hm$Salinity <- factor(
  salinity.hm$Salinity,
  levels = c(
    "non-halophile",
    "halotolerant",
    "halophile",
    "unknown"
  )
)

salinity.hm <- subset(
  salinity.hm,
  !is.na(Salinity) &
    Salinity != "unknown"
)

salinity.hm$zscore <- ave(
  salinity.hm$Prevalence,
  salinity.hm$Ligand,
  FUN = function(x) (x - mean(x)) / sd(x)
)

ggplot(salinity.hm, aes(Salinity, Ligand, fill = zscore)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(
    low = "#BEBD7F",
    mid = "white",
    high = "#592321",
    midpoint = 0
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#pH
ph.hm <- reshape2::melt(
  trait.dbb,
  id.vars = "pH",
  measure.vars = ligands,
  variable.name = "Ligand",
  value.name = "Presence"
)

overall <- aggregate((Presence > 0) ~ Ligand, data = ph.hm, FUN = mean)
names(overall)[2] <- "Overall"

ph.hm <- aggregate(
  (Presence > 0) ~ pH + Ligand,
  data = ph.hm,
  FUN = mean
)

names(ph.hm)[3] <- "Prevalence"

ph.hm <- merge(ph.hm, overall, by = "Ligand")

ph.hm$pH <- factor(
  ph.hm$pH,
  levels = c(
    "acidophile",
    "neutrophile",
    "alkaliphile"
  )
)

ph.hm <- subset(ph.hm, !is.na(pH))

ph.hm$zscore <- ave(
  ph.hm$Prevalence,
  ph.hm$Ligand,
  FUN = function(x) (x - mean(x)) / sd(x)
)

ggplot(ph.hm, aes(pH, Ligand, fill = zscore)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(
    low = "#BEBD7F",
    mid = "white",
    high = "#592321",
    midpoint = 0
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))


##GLM option
# Ligands
ligands <- c(
  "FAD",
  "Aminoacids",
  "4-hydroxycinnamic acid",
  "Biogenic amine",
  "Divalent cation"
)

# Trait-specific filters
trait.filters <- list(
  Habitat     = function(x) !is.na(x) & x != "Unknown",
  Lifestyle   = function(x) !is.na(x) & x != "Unknown",
  Oxygen      = function(x) !is.na(x) & x != "No robust majority" & x != "aerotolerant",
  Temperature = function(x) !is.na(x) & x != "No robust majority",
  Salinity    = function(x) !is.na(x) & x != "unknown",
  pH          = function(x) !is.na(x),
  Generalism  = function(x) !is.na(x)
)

# Run GLMs
glm.res <- do.call(
  rbind,
  lapply(names(trait.filters), function(trait) {
    
    do.call(
      rbind,
      lapply(ligands, function(lig) {
        
        dat <- trait.dbb[trait.filters[[trait]](trait.dbb[[trait]]), ]
        
        dat[[trait]] <- droplevels(as.factor(dat[[trait]]))
        dat$response <- as.integer(dat[[lig]] > 0)
        
        if(length(unique(dat$response)) < 2)
          return(NULL)
        
        if(nlevels(dat[[trait]]) < 2)
          return(NULL)
        
        mod <- glm(
          as.formula(paste("response ~", trait)),
          data = dat,
          family = binomial
        )
        
        out <- as.data.frame(summary(mod)$coefficients)
        
        out$Term <- rownames(out)
        out$Ligand <- lig
        out$Trait <- trait
        
        out
      })
    )
  })
)

# Rename columns
rownames(glm.res) <- NULL

names(glm.res)[names(glm.res) == "Estimate"] <- "Beta"
names(glm.res)[names(glm.res) == "Std. Error"] <- "SE"
names(glm.res)[names(glm.res) == "z value"] <- "Z"
names(glm.res)[names(glm.res) == "Pr(>|z|)"] <- "Pvalue"

# Remove intercepts
glm.res <- subset(
  glm.res,
  Term != "(Intercept)"
)

# Odds ratios
glm.res$OddsRatio <- exp(glm.res$Beta)
glm.res$OR.Lower <- exp(glm.res$Beta - 1.96 * glm.res$SE)
glm.res$OR.Upper <- exp(glm.res$Beta + 1.96 * glm.res$SE)

# FDR correction
glm.res$FDR <- p.adjust(
  glm.res$Pvalue,
  method = "BH"
)

# Significance labels
glm.res$Sig <- ""
glm.res$Sig[glm.res$FDR < 0.05] <- "*"
glm.res$Sig[glm.res$FDR < 0.01] <- "**"
glm.res$Sig[glm.res$FDR < 0.001] <- "***"

# Clean category names for plotting
glm.res$Level <- glm.res$Term
glm.res$Level <- sub("^Habitat", "", glm.res$Level)
glm.res$Level <- sub("^Lifestyle", "", glm.res$Level)
glm.res$Level <- sub("^Oxygen", "", glm.res$Level)
glm.res$Level <- sub("^Temperature", "", glm.res$Level)
glm.res$Level <- sub("^Salinity", "", glm.res$Level)
glm.res$Level <- sub("^pH", "", glm.res$Level)
glm.res$Level <- sub("^Generalism", "", glm.res$Level)

glm.res$Level <- dplyr::recode(
  glm.res$Level,
  "yes" = "Generalist",
  "no" = "Non-generalist"
)

glm.res$log2OR <- log2(glm.res$OddsRatio)

# Keep only those with significant statistical support and effect size
glm.res.sig <- subset(
  glm.res,
  FDR < 0.05 &
    abs(log2OR) > 1
)

glm.plot <- glm.res

glm.plot$Trait <- factor(
  glm.plot$Trait,
  levels = c(
    "Habitat",
    "Lifestyle",
    "Oxygen",
    "Temperature",
    "Salinity",
    "pH",
    "Generalism"
  )
)

glm.plot$log2OR.capped <- pmax(
  pmin(glm.plot$log2OR, 3),
  -3
)

glm.plot <- subset(
  glm.plot,
  !(
    Level %in% c(
      "Unknown",
      "unknown",
      "No robust majority",
      "aerotolerant"
    )
  )
)

glm.plot$Level <- dplyr::recode(
  glm.plot$Level,
  "yes" = "Generalist",
  "no" = "Non-generalist"
)

unique(glm.res$Level[glm.res$Trait == "Habitat"])
unique(glm.res$Level[glm.res$Trait == "Lifestyle"])

ggplot(glm.plot,
       aes(Level, Ligand, fill = log2OR.capped)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = Sig), size = 3) +
  scale_fill_gradient2(
    low = "#BEBD7F",
    mid = "white",
    high = "#592321",
    midpoint = 0,
    limits = c(-3,3),
    name = "log2(Odds Ratio)"
  ) +
  facet_wrap(~Trait, scales = "free_x", nrow = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


## Show all categories (even though one in each trait was chosen by the GLM as reference...)
# Start from glm results
glm.plot <- glm.res

# Clean coefficient names
glm.plot$Level <- glm.plot$Term
glm.plot$Level <- sub("^Habitat", "", glm.plot$Level)
glm.plot$Level <- sub("^Lifestyle", "", glm.plot$Level)
glm.plot$Level <- sub("^Oxygen", "", glm.plot$Level)
glm.plot$Level <- sub("^Temperature", "", glm.plot$Level)
glm.plot$Level <- sub("^Salinity", "", glm.plot$Level)
glm.plot$Level <- sub("^pH", "", glm.plot$Level)
glm.plot$Level <- sub("^Generalism", "", glm.plot$Level)

# Remove categories that should never appear
glm.plot <- subset(
  glm.plot,
  !Level %in% c(
    "Unknown",
    "unknown",
    "No robust majority",
    "aerotolerant"
  )
)

# Friendly labels
glm.plot$Level <- recode(
  glm.plot$Level,
  "yes" = "Generalist",
  "no" = "Non-generalist",
  "acidophile" = "Acidophile",
  "neutrophile" = "Neutrophile",
  "alkaliphile" = "Alkaliphile",
  "aerobic" = "Aerobic",
  "anaerobic" = "Anaerobic",
  "facultative anaerobic" = "Facultative anaerobic",
  "microaerophilic" = "Microaerophilic",
  "non-halophile" = "Non-halophile",
  "halotolerant" = "Halotolerant",
  "halophile" = "Halophile",
  "Environmental" = "Free-living (environmental)",
  "Engineered" = "Free-living (engineered)"
)

# Add the reference categories that are absent from the GLM output
ref.rows <- rbind(
  data.frame(
    Trait = "Generalism",
    Ligand = ligands,
    Level = "Non-generalist",
    log2OR = 0,
    OddsRatio = 1,
    FDR = NA,
    Sig = ""
  ),
  data.frame(
    Trait = "pH",
    Ligand = ligands,
    Level = "Neutrophile",
    log2OR = 0,
    OddsRatio = 1,
    FDR = NA,
    Sig = ""
  ),
  data.frame(
    Trait = "Oxygen",
    Ligand = ligands,
    Level = "Aerobic",
    log2OR = 0,
    OddsRatio = 1,
    FDR = NA,
    Sig = ""
  ),
  data.frame(
    Trait = "Salinity",
    Ligand = ligands,
    Level = "Non-halophile",
    log2OR = 0,
    OddsRatio = 1,
    FDR = NA,
    Sig = ""
  )
)

glm.plot <- rbind(
  glm.plot[, c(
    "Trait",
    "Ligand",
    "Level",
    "log2OR",
    "OddsRatio",
    "FDR",
    "Sig"
  )],
  ref.rows
)

# Order facets
glm.plot$Trait <- factor(
  glm.plot$Trait,
  levels = c(
    "Habitat",
    "Lifestyle",
    "Oxygen",
    "Temperature",
    "Salinity",
    "pH",
    "Generalism"
  )
)

# Biological ordering of levels
glm.plot$Level <- factor(
  glm.plot$Level,
  levels = c(
    "Animal-associated",
    "Plant-associated",
    "Food/Fermentation",
    "Free-living (engineered)",
    "Free-living (environmental)",
    "Extreme",
    
    "Generalist",
    "Non-generalist",
    
    "Aerobic",
    "Microaerophilic",
    "Facultative anaerobic",
    "Anaerobic",
    
    "Psychrophile",
    "Mesophile",
    "Thermophile",
    
    "Acidophile",
    "Neutrophile",
    "Alkaliphile",
    
    "Non-halophile",
    "Halotolerant",
    "Halophile"
  )
)

# Cap extreme ORs for cleaner visualization
glm.plot$log2OR.capped <- pmax(
  pmin(glm.plot$log2OR, 3),
  -3
)

# Plot
ggplot(glm.plot,
       aes(Level, Ligand, fill = log2OR.capped)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = Sig), size = 3) +
  scale_fill_gradient2(
    low = "#BEBD7F",
    mid = "white",
    high = "#592321",
    midpoint = 0,
    limits = c(-3, 3),
    name = "log2(OR)"
  ) +
  facet_wrap(~Trait, scales = "free_x", nrow = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



ligand_colors = c("4-hydroxycinnamic acid" : "#1B9E77",
  "Aminoacids" : "#D95F02",
  "b-type heme binding" : "#9E0202",
  "Divalent cation" :"#66A61E",
  "FAD" : "#1F78B4",
  "Biogenic amine" : "#CAB600")

## ---- Combined z-score heatmap with GLM significance overlay ----

# ============================================================
# Run GLMs
# ============================================================
glm.res <- do.call(
  rbind,
  lapply(names(trait.filters), function(trait) {
    
    do.call(
      rbind,
      lapply(ligands, function(lig) {
        
        dat <- trait.dbb[trait.filters[[trait]](trait.dbb[[trait]]), ]
        
        dat[[trait]] <- droplevels(as.factor(dat[[trait]]))
        dat$response <- as.integer(dat[[lig]] > 0)
        
        if(length(unique(dat$response)) < 2)
          return(NULL)
        
        if(nlevels(dat[[trait]]) < 2)
          return(NULL)
        
        mod <- glm(
          as.formula(paste("response ~", trait)),
          data = dat,
          family = binomial
        )
        
        out <- as.data.frame(summary(mod)$coefficients)
        
        out$Term <- rownames(out)
        out$Ligand <- lig
        out$Trait <- trait
        
        out
      })
    )
  })
)

# Rename columns
rownames(glm.res) <- NULL

names(glm.res)[names(glm.res) == "Estimate"] <- "Beta"
names(glm.res)[names(glm.res) == "Std. Error"] <- "SE"
names(glm.res)[names(glm.res) == "z value"] <- "Z"
names(glm.res)[names(glm.res) == "Pr(>|z|)"] <- "Pvalue"

# Remove intercepts
glm.res <- subset(
  glm.res,
  Term != "(Intercept)"
)

# Odds ratios
glm.res$OddsRatio <- exp(glm.res$Beta)
glm.res$OR.Lower <- exp(glm.res$Beta - 1.96 * glm.res$SE)
glm.res$OR.Upper <- exp(glm.res$Beta + 1.96 * glm.res$SE)

# FDR correction
glm.res$FDR <- p.adjust(
  glm.res$Pvalue,
  method = "BH"
)

# Significance labels
glm.res$Sig <- ""
glm.res$Sig[glm.res$FDR < 0.05] <- "*"
glm.res$Sig[glm.res$FDR < 0.01] <- "**"
glm.res$Sig[glm.res$FDR < 0.001] <- "***"

# Clean category names for plotting
glm.res$Level <- glm.res$Term
glm.res$Level <- sub("^Habitat", "", glm.res$Level)
glm.res$Level <- sub("^Lifestyle", "", glm.res$Level)
glm.res$Level <- sub("^Oxygen", "", glm.res$Level)
glm.res$Level <- sub("^Temperature", "", glm.res$Level)
glm.res$Level <- sub("^Salinity", "", glm.res$Level)
glm.res$Level <- sub("^pH", "", glm.res$Level)
glm.res$Level <- sub("^Generalism", "", glm.res$Level)

# Harmonize labels: Generalism (yes/no) + lowercase Oxygen/Temperature/Salinity/pH raw values
glm.res$Level <- dplyr::recode(
  glm.res$Level,
  "yes" = "Generalist",
  "no" = "Non-generalist",
  "aerobic" = "Aerobic",
  "anaerobic" = "Anaerobic",
  "facultative anaerobic" = "Facultative anaerobic",
  "microaerophilic" = "Microaerophilic",
  "psychrophile" = "Psychrophile",
  "mesophile" = "Mesophile",
  "thermophile" = "Thermophile",
  "acidophile" = "Acidophile",
  "neutrophile" = "Neutrophile",
  "alkaliphile" = "Alkaliphile",
  "non-halophile" = "Non-halophile",
  "halotolerant" = "Halotolerant",
  "halophile" = "Halophile"
)

glm.res$log2OR <- log2(glm.res$OddsRatio)

# Keep only those with significant statistical support and effect size
glm.res.sig <- subset(
  glm.res,
  FDR < 0.05 &
    abs(log2OR) > 1
)

glm.plot <- glm.res

glm.plot$Trait <- factor(
  glm.plot$Trait,
  levels = c(
    "Habitat",
    "Lifestyle",
    "Oxygen",
    "Temperature",
    "Salinity",
    "pH",
    "Generalism"
  )
)

glm.plot$log2OR.capped <- pmax(
  pmin(glm.plot$log2OR, 3),
  -3
)

glm.plot <- subset(
  glm.plot,
  !(
    Level %in% c(
      "Unknown",
      "unknown",
      "No robust majority",
      "aerotolerant"
    )
  )
)

# Same harmonization applied to glm.plot
glm.plot$Level <- dplyr::recode(
  glm.plot$Level,
  "yes" = "Generalist",
  "no" = "Non-generalist",
  "aerobic" = "Aerobic",
  "anaerobic" = "Anaerobic",
  "facultative anaerobic" = "Facultative anaerobic",
  "microaerophilic" = "Microaerophilic",
  "psychrophile" = "Psychrophile",
  "mesophile" = "Mesophile",
  "thermophile" = "Thermophile",
  "acidophile" = "Acidophile",
  "neutrophile" = "Neutrophile",
  "alkaliphile" = "Alkaliphile",
  "non-halophile" = "Non-halophile",
  "halotolerant" = "Halotolerant",
  "halophile" = "Halophile"
)

# ============================================================
# Build z-score heatmap data (all levels, incl. reference)
# ============================================================

build_zscore_hm <- function(trait) {
  
  df <- reshape2::melt(
    trait.dbb,
    id.vars = trait,
    measure.vars = ligands,
    variable.name = "Ligand",
    value.name = "Presence"
  )
  
  names(df)[names(df) == trait] <- "Level"
  
  df <- aggregate(
    (Presence > 0) ~ Level + Ligand,
    data = df,
    FUN = mean
  )
  names(df)[3] <- "Prevalence"
  
  df$zscore <- ave(
    df$Prevalence,
    df$Ligand,
    FUN = function(x) (x - mean(x)) / sd(x)
  )
  
  df$Trait <- trait
  df
}

traits <- names(trait.filters)

zscore.hm <- do.call(rbind, lapply(traits, build_zscore_hm))

# Drop the same "junk" categories excluded elsewhere
zscore.hm <- subset(
  zscore.hm,
  !Level %in% c("Unknown", "unknown", "No robust majority", "aerotolerant")
)

# Harmonize level labels so they match glm.res$Level exactly
zscore.hm$Level <- dplyr::recode(
  zscore.hm$Level,
  "yes" = "Generalist",
  "no" = "Non-generalist",
  "aerobic" = "Aerobic",
  "anaerobic" = "Anaerobic",
  "facultative anaerobic" = "Facultative anaerobic",
  "microaerophilic" = "Microaerophilic",
  "psychrophile" = "Psychrophile",
  "mesophile" = "Mesophile",
  "thermophile" = "Thermophile",
  "acidophile" = "Acidophile",
  "neutrophile" = "Neutrophile",
  "alkaliphile" = "Alkaliphile",
  "non-halophile" = "Non-halophile",
  "halotolerant" = "Halotolerant",
  "halophile" = "Halophile"
)

# Pull significance labels from the now-harmonized glm.res
sig.lookup <- unique(glm.res[, c("Trait", "Ligand", "Level", "Sig", "FDR")])

zscore.plot <- merge(
  zscore.hm,
  sig.lookup,
  by = c("Trait", "Ligand", "Level"),
  all.x = TRUE
)

# Levels with no GLM entry are reference categories -> blank asterisk
zscore.plot$Sig[is.na(zscore.plot$Sig)] <- ""

# Facet order
zscore.plot$Trait <- factor(
  zscore.plot$Trait,
  levels = c(
    "Habitat", "Lifestyle", "Oxygen",
    "Temperature", "Salinity", "pH", "Generalism"
  )
)

# Biological ordering of levels
zscore.plot$Level <- factor(
  zscore.plot$Level,
  levels = c(
    # Habitat
    "Animal-associated", "Plant-associated", "Food/Fermentation",
    "Freshwater", "Marine", "Soil/Sediment", "Engineered", "Extreme",
    
    # Lifestyle
    "Host-associated", "Food-associated",
    "Free-living (engineered)", "Free-living (environmental)",
    
    # Generalism
    "Generalist", "Non-generalist",
    
    # Oxygen
    "Aerobic", "Microaerophilic", "Facultative anaerobic", "Anaerobic",
    
    # Temperature
    "Psychrophile", "Mesophile", "Thermophile",
    
    # pH
    "Acidophile", "Neutrophile", "Alkaliphile",
    
    # Salinity
    "Non-halophile", "Halotolerant", "Halophile"
  )
)

# Cap extreme z-scores for cleaner color scale
zscore.plot$zscore.capped <- pmax(pmin(zscore.plot$zscore, 3), -3)

# ============================================================
# Sanity checks (run AFTER zscore.plot is built)
# ============================================================

unique(glm.res$Level[glm.res$Trait == "Oxygen"])
unique(glm.res$Level[glm.res$Trait == "Temperature"])
unique(glm.res$Level[glm.res$Trait == "Salinity"])
unique(glm.res$Level[glm.res$Trait == "pH"])

setdiff(unique(zscore.plot$Level), levels(zscore.plot$Level))  # should be empty
table(zscore.plot$Trait, is.na(zscore.plot$Sig))  # Sig should never be NA
aggregate(FDR < 0.05 ~ Trait, data = glm.res, FUN = sum)  # sanity check on significant hits per trait

# ============================================================
# Plot
# ============================================================

zscore.plot$Ligand = factor(zscore.plot$Ligand, levels = c("FAD",  "Divalent cation",
                                                           "Biogenic amine", "Aminoacids", 
                                                           "4-hydroxycinnamic acid"))

hm = ggplot(zscore.plot, aes(Level, Ligand, fill = zscore.capped)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = Sig), size = 3) +
  scale_fill_gradient2(
    low = "#BEBD7F",
    mid = "white",
    high = "#592321",
    midpoint = 0,
    limits = c(-3, 3),
    name = "Z-score"
  ) +
  facet_wrap(~Trait, scales = "free_x", nrow = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
hm

tiff("C:/CEAB/9_Histidine_kinases/Analysis/Lifestyle_ligands_heatmap.tiff", 
     units="cm", width=26, height=8, res=300)
hm
dev.off()







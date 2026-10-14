rm(list = ls())

library(dplyr)
library(ggplot2)
library(reshape2)

## ============================================================
## Load and prep data
## ============================================================
trait.dbb = read.csv(
  "C:/CEAB/9_Histidine_kinases/Analysis/Global_GTDB_ligand_trait_database.csv",
  header = T
)
head(trait.dbb)

hc <- aggregate(
  (X4.hydroxycinnamic.acid > 0) ~ Habitat,
  data = subset(trait.dbb, Habitat != "Unknown"),
  mean
)
names(hc)[2] <- "prev"
hc$z <- (hc$prev - mean(hc$prev)) / sd(hc$prev)
hc[order(-hc$z), ]

## ---- Column rename (FIXED: match() pairs old->new names correctly,
## regardless of the physical column order in trait.dbb) ----
rename_map <- c(
  "b.type.heme.binding"     = "b-type heme",
  "FAD"                     = "FAD",
  "Aminoacids"               = "Aminoacids",
  "Carboxilate"              = "Carboxylate",
  "X4.hydroxycinnamic.acid"  = "4-hydroxycinnamic acid",
  "Biogenic.amine"           = "Biogenic amine",
  "Divalent.cation"          = "Divalent cation"
)

idx <- match(names(rename_map), names(trait.dbb))
names(trait.dbb)[idx] <- rename_map

## ---- Sanity check: confirm the rename mapped correctly ----
# The raw hc$prev (from X4.hydroxycinnamic.acid, computed above) should now
# match the prevalence under the "4-hydroxycinnamic acid" label.
hc.check <- aggregate(
  `4-hydroxycinnamic acid` > 0 ~ Habitat,
  data = subset(trait.dbb, Habitat != "Unknown"),
  FUN = mean
)
names(hc.check)[2] <- "prev.after.rename"
merge(hc, hc.check, by = "Habitat")

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

ligands <- c(
  "FAD", "Aminoacids", "Carboxylate", "4-hydroxycinnamic acid",
  "Biogenic amine", "Divalent cation"
)

trait.filters <- list(
  Habitat     = function(x) !is.na(x) & x != "Unknown",
  Lifestyle   = function(x) !is.na(x) & x != "Unknown",
  Oxygen      = function(x) !is.na(x) & x != "No robust majority" & x != "aerotolerant",
  Temperature = function(x) !is.na(x) & x != "No robust majority",
  Salinity    = function(x) !is.na(x) & x != "unknown",
  pH          = function(x) !is.na(x),
  Generalism  = function(x) !is.na(x)
)

## ============================================================
## Run GLMs (needed for significance labels)
## ============================================================
glm.res <- do.call(
  rbind,
  lapply(names(trait.filters), function(trait) {
    do.call(
      rbind,
      lapply(ligands, function(lig) {
        
        dat <- trait.dbb[trait.filters[[trait]](trait.dbb[[trait]]), ]
        dat[[trait]] <- droplevels(as.factor(dat[[trait]]))
        dat$response <- as.integer(dat[[lig]] > 0)
        
        if (length(unique(dat$response)) < 2) return(NULL)
        if (nlevels(dat[[trait]]) < 2) return(NULL)
        
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

rownames(glm.res) <- NULL

names(glm.res)[names(glm.res) == "Estimate"] <- "Beta"
names(glm.res)[names(glm.res) == "Std. Error"] <- "SE"
names(glm.res)[names(glm.res) == "z value"] <- "Z"
names(glm.res)[names(glm.res) == "Pr(>|z|)"] <- "Pvalue"

glm.res <- subset(glm.res, Term != "(Intercept)")

glm.res$OddsRatio <- exp(glm.res$Beta)
glm.res$FDR <- p.adjust(glm.res$Pvalue, method = "BH")

glm.res$Sig <- ""
glm.res$Sig[glm.res$FDR < 0.05] <- "*"
glm.res$Sig[glm.res$FDR < 0.01] <- "**"
glm.res$Sig[glm.res$FDR < 0.001] <- "***"

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

## ============================================================
## Build z-score heatmap data (all levels, incl. reference)
## ============================================================
build_zscore_hm <- function(trait) {
  
  ## Apply the same filtering used elsewhere
  dat <- trait.dbb[trait.filters[[trait]](trait.dbb[[trait]]), ]
  
  out <- do.call(
    rbind,
    lapply(ligands, function(lig) {
      
      tmp <- aggregate(
        dat[[lig]] > 0,
        by = list(Level = dat[[trait]]),
        FUN = mean
      )
      
      names(tmp)[2] <- "Prevalence"
      
      tmp$Ligand <- lig
      
      tmp$zscore <- if(sd(tmp$Prevalence) > 0) {
        (tmp$Prevalence - mean(tmp$Prevalence)) /
          sd(tmp$Prevalence)
      } else {
        0
      }
      
      tmp
    })
  )
  
  out$Trait <- trait
  
  out
}
traits <- names(trait.filters)

zscore.hm <- do.call(
  rbind,
  lapply(traits, build_zscore_hm)
)

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

sig.lookup <- unique(glm.res[, c("Trait", "Ligand", "Level", "Sig", "FDR")])

zscore.plot <- merge(
  zscore.hm, sig.lookup,
  by = c("Trait", "Ligand", "Level"),
  all.x = TRUE
)
zscore.plot$Sig[is.na(zscore.plot$Sig)] <- ""

zscore.plot$Trait <- factor(
  zscore.plot$Trait,
  levels = c("Habitat", "Lifestyle", "Oxygen", "Temperature", "Salinity", "pH", "Generalism")
)

zscore.plot$Level <- factor(
  zscore.plot$Level,
  levels = c(
    "Animal-associated", "Plant-associated", "Food/Fermentation",
    "Freshwater", "Marine", "Soil/Sediment", "Engineered", "Extreme",
    "Host-associated", "Food-associated",
    "Free-living (engineered)", "Free-living (environmental)",
    "Generalist", "Non-generalist",
    "Aerobic", "Microaerophilic", "Facultative anaerobic", "Anaerobic",
    "Psychrophile", "Mesophile", "Thermophile",
    "Acidophile", "Neutrophile", "Alkaliphile",
    "Non-halophile", "Halotolerant", "Halophile"
  )
)

zscore.plot$zscore.capped <- pmax(pmin(zscore.plot$zscore, 3), -3)

zscore.plot$Ligand <- factor(
  zscore.plot$Ligand,
  levels = c("FAD", "Divalent cation", "Carboxylate", "Biogenic amine", "Aminoacids", "4-hydroxycinnamic acid")
)

## ============================================================
## Plot and save
## ============================================================
zscore.plot = subset(zscore.plot, !Ligand == "Carboxylate")
hm = ggplot(zscore.plot, aes(Level, Ligand, fill = zscore.capped)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = Sig), size = 3) +
  scale_fill_gradient2(
    low = "#BEBD7F", mid = "white", high = "#592321",
    midpoint = 0, limits = c(-3, 3), name = "Z-score"
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

tiff(
  "C:/CEAB/9_Histidine_kinases/Analysis/Lifestyle_ligands_heatmap.tiff",
  units = "cm", width = 26, height = 8, res = 300
)
hm
dev.off()

## Check z-scores for the text
subset(
  zscore.plot,
  Trait == "Oxygen" &
  Ligand == "FAD")[, c("FDR", "zscore")]

subset(
  zscore.plot,
  Trait == "Habitat" &
    Ligand == "FAD"
)[, c("Level", "zscore", "FDR")]

subset(
  zscore.plot,
  Trait == "Habitat" &
    Ligand == "Aminoacids"
)[, c("Level", "zscore", "FDR")]


subset(
  zscore.hm,
  Trait == "Habitat" &
    Ligand == "4-hydroxycinnamic acid"
)[, c("Level", "zscore")]

subset(
  zscore.hm,
  Trait == "Habitat" &
    Ligand == "Divalent cation"
)[, c("Level", "zscore")]

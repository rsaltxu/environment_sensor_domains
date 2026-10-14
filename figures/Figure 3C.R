rm(list = ls())

library(dplyr)
library(ggplot2)
library(multcompView)

## ============================================================
## Load and merge data
## ============================================================
trait.dbb = read.csv(
  "C:/CEAB/9_Histidine_kinases/Analysis/Global_GTDB_ligand_trait_database.csv",
  header = T
)

sensor.prevalence.dbb = read.csv(
  "C:/CEAB/9_Histidine_kinases/Analysis/GTDB_hk_variables.csv",
  header = T
)
names(sensor.prevalence.dbb)[names(sensor.prevalence.dbb) == "Genome"] = "accession"

complete.dbb = merge(trait.dbb, sensor.prevalence.dbb, by = "accession")

trait.filters <- list(
  Habitat     = function(x) !is.na(x) & x != "Unknown",
  Lifestyle   = function(x) !is.na(x) & x != "Unknown",
  Oxygen      = function(x) !is.na(x) & x != "No robust majority" & x != "aerotolerant",
  Temperature = function(x) !is.na(x) & x != "No robust majority",
  Salinity    = function(x) !is.na(x) & x != "unknown",
  pH          = function(x) !is.na(x),
  Generalism  = function(x) !is.na(x)
)

complete.dbb$Lifestyle <- ifelse(
  complete.dbb$Lifestyle == "Engineered",
  "Free-living (engineered)",
  complete.dbb$Lifestyle
)
complete.dbb$Lifestyle <- ifelse(
  complete.dbb$Lifestyle == "Environmental",
  "Free-living (environmental)",
  complete.dbb$Lifestyle
)

## ============================================================
## Per-trait: filter, Kruskal-Wallis, pairwise Wilcoxon, CLD
## ============================================================
build_hkmb_stats <- function(trait) {
  
  dat <- complete.dbb[trait.filters[[trait]](complete.dbb[[trait]]), ]
  dat[[trait]] <- droplevels(as.factor(dat[[trait]]))
  dat <- dat[!is.na(dat$HK.Mb), ]
  
  if (nlevels(dat[[trait]]) < 2) return(NULL)
  
  pw <- pairwise.wilcox.test(
    dat$HK.Mb,
    dat[[trait]],
    p.adjust.method = "BH"
  )
  
  pmat <- pw$p.value
  rn <- rownames(pmat)
  cn <- colnames(pmat)
  
  group.means <- aggregate(HK.Mb ~ dat[[trait]], data = dat, FUN = mean)
  names(group.means) <- c("Level", "mean")
  levs <- as.character(group.means$Level)
  
  safe.codes <- paste0("L", seq_along(levs))
  names(safe.codes) <- levs
  code.to.level <- setNames(levs, safe.codes)
  
  pvec <- c()
  for (i in seq_along(levs)) {
    for (j in seq_along(levs)) {
      if (i < j) {
        a <- levs[i]; b <- levs[j]
        p <- NA
        if (b %in% rn && a %in% cn) {
          p <- pmat[b, a]
        } else if (a %in% rn && b %in% cn) {
          p <- pmat[a, b]
        }
        safe.name <- paste(safe.codes[b], safe.codes[a], sep = "-")
        pvec[safe.name] <- p
      }
    }
  }
  pvec[is.na(pvec)] <- 1
  
  cld <- multcompLetters(pvec)$Letters
  names(cld) <- code.to.level[names(cld)]
  
  ## ---- Relabel letters so they always follow mean rank (lowest = "a") ----
  raw.letters <- unique(unlist(strsplit(cld, "")))
  
  letter.means <- sapply(raw.letters, function(L) {
    groups.with.L <- names(cld)[grepl(L, cld, fixed = TRUE)]
    mean(group.means$mean[group.means$Level %in% groups.with.L])
  })
  
  ordered.letters <- names(sort(letter.means))
  new.letter.map <- setNames(letters[seq_along(ordered.letters)], ordered.letters)
  
  remap_string <- function(s) {
    chars <- strsplit(s, "")[[1]]
    new.chars <- new.letter.map[chars]
    paste(sort(new.chars), collapse = "")
  }
  cld <- sapply(cld, remap_string)
  ## --------------------------------------------------------------------
  
  letters.df <- data.frame(
    Trait = trait,
    Level = names(cld),
    Letter = as.character(cld),
    stringsAsFactors = FALSE
  )
  
  ypos <- aggregate(HK.Mb ~ dat[[trait]], data = dat, FUN = max)
  names(ypos) <- c("Level", "ymax")
  letters.df <- merge(letters.df, ypos, by = "Level")
  
  list(
    data = data.frame(Trait = trait, Level = as.character(dat[[trait]]), HK.Mb = dat$HK.Mb),
    letters = letters.df
  )
}

traits <- names(trait.filters)
results <- lapply(traits, build_hkmb_stats)
names(results) <- traits

## ============================================================
## Combine into plotting data frames
## ============================================================

hkmb.data <- do.call(rbind, lapply(results, function(x) x$data))
hkmb.letters <- do.call(rbind, lapply(results, function(x) x$letters))

## ============================================================
## Harmonize level labels
## ============================================================

recode_levels <- function(x) {
  dplyr::recode(
    x,
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
}

hkmb.data$Level <- recode_levels(hkmb.data$Level)
hkmb.letters$Level <- recode_levels(hkmb.letters$Level)

level.order <- c(
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

hkmb.data$Level <- factor(hkmb.data$Level, levels = level.order)
hkmb.letters$Level <- factor(hkmb.letters$Level, levels = level.order)

hkmb.data$Trait <- factor(hkmb.data$Trait, levels = traits)
hkmb.letters$Trait <- factor(hkmb.letters$Trait, levels = traits)

## Sample sizes per Trait x Level
n.counts <- aggregate(HK.Mb ~ Trait + Level, data = hkmb.data, FUN = length)
names(n.counts)[3] <- "n"

hkmb.letters <- merge(hkmb.letters, n.counts, by = c("Trait", "Level"))

## Fixed y-positions for letters and sample sizes
hkmb.letters$y.letter <- rep(40, length(hkmb.letters$ymax))
hkmb.letters$y.n      <- rep(45, length(hkmb.letters$ymax))

## ============================================================
## Plot and save
## ============================================================
mb = ggplot(hkmb.data, aes(Level, HK.Mb)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width",
              fill = "grey80", colour = "grey30") +
  geom_boxplot(width = 0.12, outlier.shape = NA, colour = "grey20", fill = "white") +
  geom_text(
    data = hkmb.letters,
    aes(x = Level, y = y.letter, label = Letter),
    size = 3.2,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = hkmb.letters,
    aes(x = Level, y = y.n, label = paste0("(", n, ")")),
    size = 2.8,
    colour = "grey40",
    inherit.aes = FALSE
  ) +
  facet_grid(~Trait, scales = "free_x", space = "free_x") +
  ylim(c(0,45)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 20)
  ) +
  labs(y = "Sensor domains per Mb")
mb

tiff("C:/CEAB/9_Histidine_kinases/Analysis/Lifestyle_sensor_prevalence_violins.tiff", 
     units="cm", width=26, height=8, res=300)
mb
dev.off()


## Get means and SDs and stats
hkmb.summary <- hkmb.data %>%
  group_by(Trait, Level) %>%
  summarise(
    n = n(),
    mean = mean(HK.Mb, na.rm = TRUE),
    sd = sd(HK.Mb, na.rm = TRUE),
    .groups = "drop"
  )
hkmb.summary

hkmb.summary %>%
  mutate(text = sprintf("%.2f � %.2f", mean, sd)) %>%
  arrange(Trait, Level)

habitat.stats <- hkmb.summary %>%
  filter(Trait == "Habitat") %>%
  select(Level, mean, sd)
habitat.stats

lifestyle.stats <- hkmb.summary %>%
  filter(Trait == "Lifestyle") %>%
  select(Level, mean, sd)
lifestyle.stats

oxygen.stats <- hkmb.summary %>%
  filter(Trait == "Oxygen") %>%
  select(Level, mean, sd)
oxygen.stats

generalism.stats <- hkmb.summary %>%
  filter(Trait == "Generalism") %>%
  select(Level, mean, sd)
generalism.stats

kw.results <- lapply(traits, function(trait){
  
  dat <- complete.dbb[trait.filters[[trait]](complete.dbb[[trait]]), ]
  dat <- dat[!is.na(dat$HK.Mb), ]
  
  kw <- kruskal.test(
    as.formula(paste("HK.Mb ~", trait)),
    data = dat
  )
  
  data.frame(
    Trait = trait,
    ChiSq = unname(kw$statistic),
    df = unname(kw$parameter),
    P = kw$p.value
  )
})

kw.results <- do.call(rbind, kw.results)
kw.results

# Lifestyle
lifestyle.stats %>%
  mutate(txt = sprintf("%s (%.2f � %.2f)", Level, mean, sd))

# Habitat
habitat.stats %>%
  mutate(txt = sprintf("%s (%.2f � %.2f)", Level, mean, sd))

# Oxygen
oxygen.stats %>%
  mutate(txt = sprintf("%s (%.2f � %.2f)", Level, mean, sd))

# Generalism
generalism.stats %>%
  mutate(txt = sprintf("%s (%.2f � %.2f)", Level, mean, sd))

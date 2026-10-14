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

## Ligand annotations by GTDB genome
ligands.gtdb = fread("C:/CEAB/9_Histidine_kinases/Analysis/GTDB_metadata_ligands.csv",
                      header = T)
head(ligands.gtdb)

ligands.gtdb <- ligands.gtdb %>%
  separate(
    gtdb_taxonomy,
    into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = ";\\s*",
    fill = "right"
  ) %>%
  mutate(
    Domain = sub("^d__", "", Domain),
    Phylum = sub("^p__", "", Phylum),
    Class  = sub("^c__", "", Class),
    Order  = sub("^o__", "", Order),
    Family = sub("^f__", "", Family),
    Genus  = sub("^g__", "", Genus),
    Species= sub("^s__", "", Species)
  )
length(unique(ligands.gtdb$Species)) #51343 species
head(ligands.gtdb)

## metaTraits lifestyle classifications
meta.traits = fread("C:/CEAB/9_Histidine_kinases/Analysis/gtdb_species_summary_all.tsv.gz", header = T)
meta.traits = meta.traits[, c("Domain", "Phylum", "Class", "Order",
                "Family", "Genus", "Species") :=
              tstrsplit(taxon_lineage, "|", fixed = TRUE)]
length(unique(meta.traits$Species)) #67432 species
head(meta.traits)

# Map those taxa with lifestyle information on the ligand database
ligands.gtdb.traits = merge(ligands.gtdb, meta.traits, by = "Species")
ligands.gtdb.traits = setDT(ligands.gtdb.traits)
length(unique(ligands.gtdb.traits$Species)) #31635 species
head(ligands.gtdb.traits)

# Define traits to look for (environmental preferences and lifestyles)
oxygen = oxygen <- unique(
  ligands.gtdb.traits[
    trait_name == "oxygen preference",
    .(Species, Oxygen = consensus_value)
  ]
)
head(oxygen)

temperature <- unique(
  ligands.gtdb.traits[
    trait_name == "temperature preference",
    .(Species,
      Temperature = consensus_value)
  ]
)
head(temperature)

salinity <- unique(
  ligands.gtdb.traits[
    trait_name == "salinity preference",
    .(Species, salinity_raw = consensus_value)
  ]
)

salinity[,
  Salinity := fifelse(
    salinity_raw %in% c("halophile", "hypersaline"),
    "halophile",
    fifelse(
      salinity_raw %in% c("halotolerant", "euryhaline"),
      "halotolerant",
      fifelse(
        salinity_raw == "non-halophilic",
        "non-halophile",
        "unknown"
      )
    )
  )
]
head(salinity)

ph <- unique(
  ligands.gtdb.traits[
    trait_name == "pH preference",
    .(Species, pH_raw = consensus_value)
  ]
)

ph[
  ,
  pH := fifelse(
    pH_raw == "acidophile",
    "acidophile",
    fifelse(
      pH_raw == "alkaliphile",
      "alkaliphile",
      "neutrophile"
    )
  )
]
head(ph)

generalism <- unique(
  ligands.gtdb.traits[
    trait_name == "generalist",
    .(
      Species,
      Generalism = consensus_value
    )
  ]
)

generalism[
  ,
  Generalism := fifelse(
    Generalism == "true",
    "yes",
    fifelse(
      Generalism == "false",
      "no",
      NA_character_
    )
  )
]
head(generalism)

# Create categorical lifestyle columns from isolation source
animal_pat <- paste(
  c(
    "gut","hindgut","foregut","gastrointestinal","gastro-intestinal",
    "intestinal","intestine","jejunum","cecum","caecum","ceacal",
    "ceca","colon","feces","faeces","fecal","stool","rumen",
    "digestive tract","gut content","intestinal content",
    "alimentary canal","stomach","gastric","mucosa",
    "oral","mouth","oral cavity","dental","plaque","saliva",
    "tongue","throat","pharynx","oropharynx","nasopharynx",
    "nasopharyngeal","nasal cavity","nasal swab","throat swab",
    "tonsil","tonsil swab","trachea",
    "blood","urine","urinary tract","vagina","vaginal",
    "uterus","placenta","skin","wound","abscess",
    "sputum","lung","bronch","bal",
    "kidney","liver","eye","ocular",
    "cerebrospinal fluid","gills","hepatopancreas",
    "udder","mucus","biopsy","tissue",
    "clinical","patient","hospital",
    "human","mouse","rat","pig","cow","cattle","horse",
    "dog","cat","bird","chicken","fish","shrimp",
    "crab","sheep","equine",
    "termite","bee","ant","wasp","insect","tick",
    "beetle","malpighian","nematode","snail","larvae",
    "mussel","mussels","oncorhynchus",
    "trophosome","host","host stool","animal manure",
    "cloacal swab"
  ),
  collapse="|"
)

plant_pat <- paste(
  c(
    "rhizosphere","rhizoplane",
    "root","roots",
    "root nodule","root nodules",
    "nodule","nodules",
    "endophyte","endosphere",
    "plant","tree",
    "leaf","leaves",
    "stem","shoot",
    "bark","wood","sap",
    "phyllosphere",
    "flower","fruit",
    "seed","seedling",
    "rhizome","tuber","tubers",
    "maize field","rice field",
    "crop","crown gall"
  ),
  collapse="|"
)

food_pat <- paste(
  c(
    "food",
    "food production",
    "food processing",
    "cheese",
    "milk","raw milk",
    "honey",
    "kimchi",
    "pickle",
    "sourdough",
    "wheat sourdough",
    "silage",
    "jeotgal",
    "tea",
    "pu er tea",
    "apple juice",
    "cider",
    "fermentation",
    "fermented"
  ),
  collapse="|"
)

engineered_pat <- paste(
  c(
    "activated sludge","sludge",
    "wastewater","waste water",
    "wastewater treatment",
    "treatment plant",
    "sewage","wwtp",
    "reactor","bioreactor",
    "digester","digestate",
    "anaerobic digestion",
    "anaerobic digester",
    "anammox","ebpr",
    "biofilter","biogas","chemostat",
    "hospital shower hose",
    "household tap filter",
    "air conditioner",
    "spacecraft hardware",
    "clean room",
    "anode biofilm",
    "granular activated carbon"
  ),
  collapse="|"
)

extreme_pat <- paste(
  c(
    "hot spring",
    "acid mine","acid mine drainage",
    "acidic spring",
    "alkaline lake","alkaline salt lake",
    "hypersaline",
    "soda lake",
    "saltern","salt pond",
    "geothermal",
    "hydrothermal",
    "hydrothermal vent",
    "hydrothermal deposit",
    "hydrothermal plume",
    "sulfide chimney",
    "glacier","glacial","ice core",
    "glacier ice","snow",
    "ice-wedge polygon",
    "desert","desert environment",
    "desert regolith",
    "desert biocrust",
    "oilfield","oil field",
    "crude oil",
    "coal mine",
    "manganese mine",
    "volcan"
  ),
  collapse="|"
)

marine_pat <- paste(
  c(
    "marine","marine metagenome",
    "marine environment",
    "sea","sea water","seawater",
    "ocean","coastal",
    "estuary","estuarine",
    "intertidal",
    "tidal flat",
    "salt marsh",
    "mangrove",
    "coral",
    "sponge","marine sponge",
    "cold seep",
    "deep.?sea",
    "abyssal",
    "sea ice",
    "kelp forest",
    "phytoplankton",
    "phycosphere",
    "cyanobacterial bloom",
    "cyanobacterial mat",
    "cyanobacterial aggregates"
  ),
  collapse="|"
)

freshwater_pat <- paste(
  c(
    "lake","lake water","lake sediment",
    "freshwater","fresh water",
    "river","stream","creek",
    "pond","reservoir",
    "wetland",
    "groundwater",
    "aquifer",
    "well water",
    "spring water",
    "waterfall",
    "water column",
    "surface water",
    "water sample",
    "drinking water",
    "tap water",
    "bay water",
    "brackish water",
    "mineral water",
    "running water",
    "hyporheic zone"
  ),
  collapse="|"
)

soil_pat <- paste(
  c(
    "soil","soil sample",
    "forest soil","garden soil",
    "agricultural soil","paddy soil",
    "mineral soil","saline soil",
    "desert soil",
    "grassland",
    "forest",
    "topsoil",
    "cave soil",
    "cave",
    "permafrost",
    "peat",
    "sediment","sediments",
    "lacustrine sediment",
    "mud",
    "sand","beach sand",
    "rock",
    "stone",
    "granite",
    "clay",
    "terrestrial",
    "stromatolite",
    "lichen",
    "microbial mat",
    "microbial mat material",
    "biofilm",
    "compost"
  ),
  collapse="|"
)

ligands.gtdb.traits[, iso := tolower(trimws(ncbi_isolation_source))]

ligands.gtdb.traits[
  ,
  Habitat := fcase(
    grepl(animal_pat, iso, perl = TRUE),     "Animal-associated",
    grepl(plant_pat, iso, perl = TRUE),      "Plant-associated",
    grepl(food_pat, iso, perl = TRUE),       "Food/Fermentation",
    grepl(engineered_pat, iso, perl = TRUE), "Engineered",
    grepl(extreme_pat, iso, perl = TRUE),    "Extreme",
    grepl(marine_pat, iso, perl = TRUE),     "Marine",
    grepl(freshwater_pat, iso, perl = TRUE), "Freshwater",
    grepl(soil_pat, iso, perl = TRUE),       "Soil/Sediment",
    default = "Unknown"
  )
]

ligands.gtdb.traits[,
  Lifestyle := fcase(
    Habitat %in% c("Animal-associated", "Plant-associated"),
    "Host-associated",
    Habitat == "Food/Fermentation",
    "Food-associated",
    Habitat == "Engineered",
    "Engineered",
    Habitat %in% c(
      "Marine",
      "Freshwater",
      "Soil/Sediment",
      "Extreme"
    ),
    "Environmental",
    default = "Unknown"
  )
]

head(ligands.gtdb.traits)
length(unique(ligands.gtdb.traits$Species)) #31635 species

## Merge all trait databases
oxygen[, Oxygen := fifelse(
  Oxygen %in% c("obligate aerobic", "aerobic"), "aerobic",
  fifelse(
    Oxygen %in% c("obligate anaerobic", "anaerobic"), "anaerobic",
    Oxygen
  )
)]

temperature[, Temperature := fifelse(
  Temperature == "mesophilic", "mesophile",
  fifelse(
    Temperature %in% c("thermophilic", "hyperthermophilic"),
    "thermophile",
    fifelse(
      Temperature == "psychrophilic",
      "psychrophile",
      Temperature
    )
  )
)]

ph[pH == "No robust majority", pH := "neutrophile"]

species.traits <- Reduce(
  function(x, y) merge(x, y, by = "Species", all = TRUE),
  list(
    oxygen,
    temperature,
    salinity,
    ph,
    generalism
  )
)
head(species.traits)
nrow(species.traits) #13216 genomes have trait data

setDT(ligands.gtdb.traits)
setDT(species.traits)

global.trait.database = ligands.gtdb.traits[
  species.traits,
  `:=`(
    Oxygen = i.Oxygen,
    Temperature = i.Temperature,
    Salinity = i.Salinity,
    pH = i.pH,
    Generalism = i.Generalism
  ),
  on = "Species"
]
head(global.trait.database)
length(unique(global.trait.database$Species)) #31635 genomes, correct

global.trait.database[
  ,
  lapply(.SD, function(x) sum(!is.na(x))),
  .SDcols = c("Oxygen","Temperature","Salinity","pH","Generalism")
]

keep.cols = c("accession", "Phylum.y",
              "Class.y", "Order.y", "Family.y", "Genus.y", "Species", "genome_size", "gc_percentage",
              "checkm2_completeness", "checkm2_contamination", "ncbi_genome_category",
              "ncbi_genbank_assembly_accession", "ncbi_lat_lon", "ncbi_isolation_source",
              "protein_count", "b-type heme binding", "Purine", "FAD", "Aminoacids",
              "4-hydroxycinnamic acid", "Biogenic amine", "Carboxilate", "Divalent cation",
              "c-type heme binding", "taxon_id", "taxon_name", "Habitat", "Lifestyle",
              "Oxygen", "Temperature", "Salinity", "pH", "Generalism")
global.trait.database <- global.trait.database[, ..keep.cols]
head(global.trait.database)

global.trait.database <- unique(global.trait.database)
head(global.trait.database)

#write.csv(global.trait.database, "C:/CEAB/9_Histidine_kinases/Analysis/Global_GTDB_ligand_trait_database.csv")

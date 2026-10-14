# environment_sensor_domains

This repository contains the code and analysis workflows accompanying the study **“Environmental sensing capacity predicts bacterial ecological strategies and environmental preferences.”**

The project investigates how bacterial environmental sensing capacity relates to ecological strategies, habitat preferences, and adaptation to environmental conditions. Using thousands of bacterial genomes and soil metagenomes, the analyses characterize sensor-domain repertoires and experimentally characterized ligand-binding motifs to identify signatures associated with factors such as habitat, oxygen availability, generalist strategies, and aridity.

The repository includes the computational workflows used to process genomic and metagenomic data, characterize sensory repertoires, perform comparative analyses, and generate the figures and results presented in the manuscript.


## Repository structure

| Directory | Description |
|---|---|
| `01_metagenome_assembly_taxonomic_profiling/` | Assembly of metagenomes and taxonomic profiling |
| `02_bacterial_domain_annotation/` | Annotation and profiling of bacterial sensor domains |
| `03_ecological_phylogenetic_analysis/` | Scripts for sensor domain distribution analysis |
| `04_metatraits/` | Scripts for metaTraits |
| `05_statistical_analysis/` | Additional statistical analyses |
| `06_sequence_embedding/` | Encoding of protein sequences |
| `07_machine_learning/` | SVM-based feature importance analysis |
| `figures/` | Code to replicate figures |
| `metadata/` | Additional metadata |
| `README.md` | Project documentation |

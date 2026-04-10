# 16S rRNA analysis pipeline (DADA2)

This repository contains the R pipeline used in:

Di Nezio et al. – Multilayered human activities shape the microbial communities of groundwater-dependent ecosystems on an arid oceanic island

## Description
Pipeline for processing 16S rRNA amplicon data using DADA2, including:
- quality filtering and ASV inference
- taxonomic assignment (SILVA v138.2)
- phyloseq object construction
- diversity analyses
- functional group classification

## Requirements
- R >= 4.0
- dada2
- phyloseq
- vegan
- cutadapt

## Input files
- Paired-end FASTQ files
- SILVA reference database
- metadata table
- trait annotation table

## Outputs
- ASV table
- taxonomy table
- relative abundance tables
- figures

## Notes
Edit file paths before running the pipeline.

########################################
## 16S Analysis Pipeline
## DADA2, taxonomic assignment, phyloseq
## relative-abundance table generation
########################################

### 1. Load libraries ----

library(knitr)
library(ggplot2)
library(gridExtra)
library(dada2)
library(phangorn)
library(phyloseq)
library(DECIPHER)
library(BiocStyle)
library(Biostrings)
library(vegan)
library(usedist)
library(stringr)
library(dplyr)
library(tidyr)
library(purrr)
library(readxl)
library(writexl)

### 2. Define paths ----

# Set the working directory containing the raw FASTQ files
path <- "PATH/TO/RAW_FASTQ_FILES"

# Path to cutadapt executable
cutadapt_path <- "PATH/TO/cutadapt.exe"

# Reference databases for taxonomy assignment
silva_trainset <- "PATH/TO/silva_nr99_v138.2_toGenus_trainset.fa.gz"
silva_species  <- "PATH/TO/silva_v138.2_assignSpecies.fa.gz"

# Trait annotation table
patho_list_file <- "PATH/TO/PATHOGEN_LIST_FILE"

setwd(path)

### 3. Identify forward and reverse reads ----

fnFs <- sort(list.files(path, pattern = "_R1\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2\\.fastq\\.gz$", full.names = TRUE))

stopifnot(length(fnFs) > 0, length(fnRs) > 0, length(fnFs) == length(fnRs))

sample.names <- sub("_R1\\.fastq\\.gz$", "", basename(fnFs))
stopifnot(!any(duplicated(sample.names)))

names(fnFs) <- sample.names
names(fnRs) <- sub("_R2\\.fastq\\.gz$", "", basename(fnRs))

if (!identical(names(fnFs), names(fnRs))) {
  stop("Forward and reverse read files do not match by sample name.")
}

fnRs <- fnRs[names(fnFs)]
sample.names <- names(fnFs)

### 4. Primer trimming and read orientation with cutadapt ----

FWDn <-    # Forward primer sequence
REVn <-    # Reverse primer sequence

rc <- dada2::rc
Oreq <- as.character(floor(0.8 * min(nchar(FWDn), nchar(REVn))))

stopifnot(file.exists(cutadapt_path))
system2(cutadapt_path, "--version", stdout = TRUE, stderr = TRUE)

run_cutadapt_pair <- function(inF, inR, outF, outR, cores = 6, logFile = NULL) {
  args <- c(
    "-e", "0.2",
    "-O", Oreq,
    "--pair-filter=both",
    "--discard-untrimmed",
    "--max-n", "0",
    "-m", "50",
    "-g", FWDn,
    "-g", rc(REVn),
    "-G", REVn,
    "-G", rc(FWDn),
    "--cores", as.character(cores),
    "-o", outF,
    "-p", outR,
    inF, inR
  )
  
  out <- system2(cutadapt_path, args = args, stdout = TRUE, stderr = TRUE)
  if (!is.null(logFile)) writeLines(out, con = logFile)
  out
}

trimdir <- file.path(path, "trimmed_oriented")
dir.create(trimdir, showWarnings = FALSE)

trimFs <- file.path(trimdir, paste0(sample.names, "_F_trim.fastq.gz"))
trimRs <- file.path(trimdir, paste0(sample.names, "_R_trim.fastq.gz"))

logdir <- file.path(trimdir, "logs")
dir.create(logdir, showWarnings = FALSE)

for (i in seq_along(sample.names)) {
  run_cutadapt_pair(
    fnFs[i], fnRs[i],
    trimFs[i], trimRs[i],
    cores = 6,
    logFile = file.path(logdir, paste0(sample.names[i], "_cutadapt.log"))
  )
}

stopifnot(all(file.exists(trimFs)), all(file.exists(trimRs)))

fnFs <- trimFs
fnRs <- trimRs
names(fnFs) <- sample.names
names(fnRs) <- sample.names

### 5. Quality filtering ----

filt_path <- file.path(path, "filtered")
dir.create(filt_path, showWarnings = FALSE)

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(265, 215),
  maxN = 0,
  maxEE = c(2, 4),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

save.image("filtered.RData")

### 6. Learn error rates ----

errF <- learnErrors(filtFs, multithread = TRUE, verbose = 1)
save.image("sampleserrf.RData")

errR <- learnErrors(filtRs, multithread = TRUE, verbose = 1)
save.image("sampleserrR.RData")

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

### 7. Dereplication and ASV inference ----

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
save.image("samplesdadaf.RData")

dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
save.image("samplesdadaR.RData")

### 8. Merge paired reads ----

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
save.image("samplesmer.RData")

seqtab <- makeSequenceTable(mergers)

### 9. Remove chimeras and retain target amplicon length ----

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

save.image("samplesnc.RData")

table(nchar(getSequences(seqtab.nochim)))

seqtab.nochim2 <- seqtab.nochim[, nchar(colnames(seqtab.nochim)) %in% 401:430]

sum(seqtab.nochim2) / sum(seqtab)

### 10. Track reads through the pipeline ----

getN <- function(x) sum(getUniques(x))

track <- cbind(
  out,
  sapply(dadaFs, getN),
  sapply(dadaRs, getN),
  sapply(mergers, getN),
  rowSums(seqtab.nochim),
  rowSums(seqtab.nochim2)
)

colnames(track) <- c(
  "input",
  "filtered",
  "denoisedF",
  "denoisedR",
  "merged",
  "nonchim",
  "length_reduced"
)

rownames(track) <- sample.names

write.csv(track, "read_tracking.csv")

### 11. Save sequence table ----

write.csv(seqtab.nochim2, "sequence_tab.csv")
save.image("samples_beforeremove.RData")

rm(mergers, dadaFs, dadaRs, seqtab, errF, errR, seqtab.nochim)

### 12. Taxonomic assignment ----

taxa <- assignTaxonomy(seqtab.nochim2, silva_trainset, multithread = TRUE)
save.image("samples_taxa.RData")

taxa2 <- addSpecies(taxa, silva_species)

write.csv(taxa2, "taxonomy.csv")

taxa.print <- taxa2
rownames(taxa.print) <- NULL
head(taxa.print)

### 13. Create phyloseq object ----

seqtab.nochim2 <- as.matrix(seqtab.nochim2)

seqs <- colnames(seqtab.nochim2)
asv_ids <- paste0("ASV", seq_along(seqs))

colnames(seqtab.nochim2) <- asv_ids
tseq <- t(seqtab.nochim2)

vari <- read.csv("taxonomy.csv", check.names = FALSE)
vari <- as.data.frame(vari)
rownames(vari) <- vari[, 1]
vari <- vari[, -1, drop = FALSE]

vari <- vari[seqs, , drop = FALSE]
rownames(vari) <- asv_ids

ps <- phyloseq(
  otu_table(tseq, taxa_are_rows = TRUE),
  tax_table(as.matrix(vari))
)

dna <- DNAStringSet(seqs)
names(dna) <- asv_ids
ps <- merge_phyloseq(ps, dna)

ps <- subset_taxa(
  ps,
  !is.na(Phylum) &
    Phylum != "" &
    !Phylum %in% c("uncharacterized", "Unassigned")
)

ps <- subset_taxa(
  ps,
  !(Class %in% "Chloroplast") &
    !(Family %in% "Mitochondria")
)

### 14. Relative abundance table ----

ps_rel <- transform_sample_counts(
  ps,
  function(x) if (sum(x) == 0) x else x / sum(x)
)

otu_rel <- as(otu_table(ps_rel), "matrix")
if (!taxa_are_rows(ps_rel)) otu_rel <- t(otu_rel)

tax_mat <- as(tax_table(ps_rel), "matrix")
tax_mat[is.na(tax_mat) | tax_mat == ""] <- "Unclassified"

otu_rel <- otu_rel[rownames(tax_mat), , drop = FALSE]

relabund_table <- cbind(
  ASV = rownames(tax_mat),
  as.data.frame(tax_mat, stringsAsFactors = FALSE),
  as.data.frame(otu_rel, stringsAsFactors = FALSE)
)

### 15. Summarize relative abundances by taxon ----

is_alphanumeric <- function(x) grepl("[0-9]", x)

pick_taxon <- function(Kingdom, Phylum, Class, Order, Family, Genus, Species) {
  ranks <- c(Species, Genus, Family, Order, Class, Phylum, Kingdom)
  names(ranks) <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom")
  
  ranks <- ranks[ranks != "Unclassified"]
  
  if (length(ranks) == 0) return("Unclassified")
  
  if (names(ranks)[1] == "Species") {
    g <- if ("Genus" %in% names(ranks)) ranks["Genus"] else NA_character_
    sp <- ranks["Species"]
    if (!is.na(g) && g != "Unclassified") return(str_squish(paste(g, sp)))
    return(sp)
  }
  
  deepest <- ranks[1]
  if (is_alphanumeric(deepest) && length(ranks) >= 2) return(ranks[2])
  return(deepest)
}

taxtable_summ <- relabund_table %>%
  mutate(across(Kingdom:Species, ~ na_if(.x, ""))) %>%
  mutate(across(Kingdom:Species, ~ coalesce(.x, "Unclassified"))) %>%
  group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  mutate(
    Taxon = pmap_chr(
      list(Kingdom, Phylum, Class, Order, Family, Genus, Species),
      pick_taxon
    )
  )

writexl::write_xlsx(taxtable_summ, "taxtable_sum.xlsx")

### 16. Add pathogenicity and trait annotations ----

tax_table <- readxl::read_excel("taxtable_sum.xlsx")

Patho_List <- readxl::read_xlsx(patho_list_file, sheet = 1)

Taxa_Patho <- tax_table %>%
  left_join(
    Patho_List %>%
      select(Taxon, Pathogen, Pathogenicity_Class, Origin, Ecological_Trait),
    by = "Taxon")


### 17. Export unmatched taxa ----

Non_matching_taxa <- subset(
  Taxa_Patho,
  is.na(Pathogen) & is.na(Pathogenicity_Class) & is.na(Origin)
)

write.csv2(
  Non_matching_taxa,
  "Non_matching_taxa.csv",
  row.names = FALSE
)
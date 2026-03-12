
setwd("./Paper_Code/")

suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(readr)
})

# Parameters (change if needed)
species_dataset <- "hsapiens_gene_ensembl"  
ensembl_host <- "https://rest.ensembl.org"
min_length_to_filter <- 500
output_file <- "genes_to_filter.txt"

# Unwanted biotypes
unwanted_biotypes <- c(
  "lncRNA", "lincRNA", "antisense", "sense_intronic", "sense_overlapping", "non_coding",
  "processed_transcript", "macro_lncRNA", "ncRNA",
  "snRNA", "snoRNA", "scaRNA", "miRNA", "misc_RNA", "rRNA", "tRNA", "Mt_rRNA", "Mt_tRNA",
  "pseudogene", "processed_pseudogene", "unprocessed_pseudogene", "transcribed_unprocessed_pseudogene",
  "polymorphic_pseudogene", "unitary_pseudogene", "transcribed_processed_pseudogene", "rRNA_pseudogene",
  "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene",
  "readthrough_transcript"
)

ensembl <- NULL
try({
  ensembl <- useEnsembl(biomart = "genes", dataset = species_dataset)
}, silent = TRUE)
if (is.null(ensembl)) {
  ensembl <- useMart("ensembl", dataset = species_dataset)
}

# Attributes to retrieve
attrs <- c("ensembl_gene_id", "external_gene_name", "gene_biotype", "start_position", "end_position", "percentage_gene_gc_content")
genes <- getBM(attributes = attrs, mart = ensembl)

# Compute gene length
genes <- genes %>%
  mutate(gene_length = abs(end_position - start_position) + 1) %>%
  rename(ensembl_id = ensembl_gene_id, symbol = external_gene_name, biotype = gene_biotype)

# Some genes may lack symbol; keep them but with NA symbol
# Build filter condition
genes_to_filter <- genes %>%
  filter(
    (biotype %in% unwanted_biotypes) |
      (gene_length < min_length_to_filter) |
      (percentage_gene_gc_content > 70) |
      (grepl('^HIST', symbol, ignore.case = TRUE))
  ) %>%
  filter(
    (biotype %in% unwanted_biotypes) |
      (gene_length < min_length_to_filter) |
      (grepl('^HIST', symbol, ignore.case = TRUE))
  ) %>%
  select(symbol, ensembl_id) %>%
  distinct() %>%
  arrange(is.na(symbol), symbol)

# Write output
write_tsv(genes_to_filter, output_file, col_names = TRUE)

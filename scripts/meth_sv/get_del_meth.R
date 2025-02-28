args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Error: Please provide 'working_directory' and 'Sample id' as arguments.", call. = FALSE)
}

# Assign arguments to variables
wk_dir <- args[1]
SampleID <- args[2]


suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(rtracklayer)))
source(paste0(wk_dir, "/scripts/meth_sv/sv_extract_support.R"))

sample_folder = paste0(wk_dir, "/output/", SampleID)

sv_merged <- readRDS(paste0(sample_folder, "/", SampleID, "_sv_merged.RDS"))

read_mapping <- readRDS(paste0(sample_folder, "/", SampleID, "_read_mapping.RDS"))

modkit_extract <- fread(paste0(sample_folder, "/", SampleID, "_merged_SV_modkit_extract.tsv")) %>%
  filter(mod_code == "m") %>%
  dplyr::select(read_id, forward_read_position, ref_position,chrom, mod_strand, ref_strand, read_length, mod_qual, base_qual)

DEL <- sv_merged %>%
  filter(str_detect(ALT, "DEL"))


DEL_reads <- DEL %>%
  left_join(read_mapping%>%
              dplyr::select(ID, read_id, read_alt_pos, read_strand, read_type)) %>%
  left_join(modkit_extract) 


saveRDS(DEL_reads, file = paste0(sample_folder, "/", SampleID, "_DEL_meth.RDS"))
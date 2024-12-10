args = commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
  stop("Error: Please provide 'working_directory' and 'Sample id' as arguments.", call. = FALSE)
}

# Assign arguments to variables
wk_dir <- args[1]
SampleID <- args[2]


suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyverse)))

sample_folder = paste0(wk_dir, "/output/", SampleID)

sv_merged <- readRDS(paste0(sample_folder, "/", SampleID, "_sv_merged.RDS"))

supporting_read_mapping_nanomonsv <- fread(paste0(sample_folder, "/SV/", SampleID, ".nanomonsv.supporting_read.txt"),header = F) %>%
  `colnames<-`(c("SV_start_chr", "SV_start_pos","SV_start_strand", "SV_end_chr", "SV_end_pos", "SV_end_strand", "alt_seq", "ID", "read_id", "read_alt_pos", "read_strand")) %>%
  filter(ID %in% sv_merged$ID)

read_mapping <- fread(paste0(sample_folder, "/", SampleID, "_merged_sv_readid.tsv"), col.names = c("sv","read_id")) %>%
  separate("sv", into = c("#CHROM", "SV_start", "SV_end")) %>%
  mutate(SV_start = as.numeric(SV_start)+1,
         SV_end = as.numeric(SV_end))%>%
  left_join(sv_merged %>%
              mutate(SV_end = if_else(is.na(SV_end), SV_start, SV_end))%>%
              dplyr::select(`#CHROM`, SV_start, SV_end, ID) ) %>%
  separate_rows(read_id, sep = ",") %>%
  distinct() %>%
  left_join(supporting_read_mapping_nanomonsv %>%
              dplyr::select("ID", "read_id", "read_alt_pos", "read_strand") ) %>%
  mutate(read_type = if_else(!is.na(read_alt_pos), "alt", "ref"))


saveRDS(read_mapping, file = paste0(sample_folder, "/", SampleID, "_read_mapping.RDS"))
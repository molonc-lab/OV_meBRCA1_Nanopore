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

sv_severus_supporting_read <-  fread(paste0(sample_folder, "/severus/", "read_ids.csv") )%>%
   `colnames<-`(c("ID", "T_READ_ID")) %>%
  filter(ID %in% sv_merged$ID) %>%
  separate_rows(T_READ_ID, sep = "; ")


read_mapping <- fread(paste0(sample_folder, "/", SampleID, "_merged_sv_readid.tsv"), col.names = c("sv","read_id")) %>%
  separate("sv", into = c("#CHROM", "SV_start", "SV_end")) %>%
  mutate(SV_start = as.numeric(SV_start)+1,
         SV_end = as.numeric(SV_end))%>%
  left_join(sv_merged %>%
              mutate(SV_end = if_else(is.na(SV_end), SV_start+1, SV_end)) )%>%
  separate_rows(read_id, sep = ",") %>%
  distinct() %>%
  mutate(read_type = if_else(read_id %in% sv_severus_supporting_read$T_READ_ID, "alt", "ref"))


saveRDS(read_mapping, file = paste0(sample_folder, "/", SampleID, "_read_mapping.RDS"))
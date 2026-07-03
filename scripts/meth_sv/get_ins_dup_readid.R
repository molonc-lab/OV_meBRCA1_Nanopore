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

read_mapping <- readRDS(paste0(sample_folder, "/", SampleID, "_read_mapping.RDS"))

sv_merged %>%
  filter(SV_type == "INS") %>%
  arrange(ID) %>%
  dplyr::select(ID) %>%
  left_join(read_mapping%>%
            dplyr::select(ID, read_id, read_type)) %>%
  filter(read_type == "alt") %>%
  pull(read_id) %>%
  unique(.)%>%
  write_lines(paste0(sample_folder, "/ins_read_name.tsv"))

sv_merged %>%
  filter(SV_type == "DUP") %>%
  arrange(ID) %>%
  dplyr::select(ID) %>%
  left_join(read_mapping%>%
              dplyr::select(ID, read_id, read_type)) %>%
  filter(read_type == "alt") %>%
  pull(read_id) %>%
  unique(.)%>%
  write_lines(paste0(sample_folder, "/dup_read_name.tsv"))



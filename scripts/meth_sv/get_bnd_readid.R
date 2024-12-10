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

BND <- sv_merged %>%
  filter(str_detect(ALT, "\\]")|str_detect(ALT, "\\[")) %>%
  arrange(ID) %>%
  mutate(chr_end = str_extract(ALT, "chr[0-9XYM]+"),
         SV_end =as.numeric(str_extract(ALT, "(?<=:)[0-9]+"))) %>%
  dplyr::select(ID, ALT, chr_start = `#CHROM`, SV_start, chr_end, SV_end, BND_type)


BND %>%
  left_join(read_mapping%>%
              dplyr::select(ID, read_id, read_alt_pos, read_strand, read_type)) %>%
  filter(read_type == "alt") %>%
  pull(read_id) %>%
  unique(.)%>%
  write_lines(paste0(sample_folder,"/","bnd_read_name_saf_tag.tsv")

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

ins_primary_mapping <- fread(paste0(sample_folder,"/ins_read_name_primary_map.tsv"), col.names = c("read_id",  "flag", "rname", "pos", "cigar")) %>%
  mutate(strand = if_else(flag %in% c(16, 2064), "-", "+")) %>%
  rowwise() %>%
  mutate(poses = process_cigar_positions(strand, as.numeric(pos), cigar)) %>%
  separate(poses, into = c("ref_start", "ref_end", "read_start", "read_end")) %>%
  mutate(ref_start = as.numeric(ref_start), ref_end = as.numeric(ref_end),
         read_start = as.numeric(read_start), read_end = as.numeric(read_end)) %>%
  dplyr::select(-pos, -flag)

ins_supporting_reads <- ins_primary_mapping %>%
  left_join(modkit_extract %>%
              dplyr::select(read_id, forward_read_position, read_length, mod_qual, base_qual, ref_position) %>%
              distinct(),
            by =join_by(read_id == read_id, read_start <= forward_read_position, read_end >= forward_read_position)) %>%
  dplyr::rename(chrom = rname, ref_strand = strand) %>%
  filter(!is.na(forward_read_position) ) %>%
  mutate(ref_start = if_else(ref_strand == "-", ref_end, ref_start)) %>%
  left_join(read_mapping %>%
              dplyr::select(read_id, ID) %>%
              left_join(sv_merged %>%
              dplyr::select(ID, SV_start)) %>%
              distinct()) 

ins_shift = setDT(ins_supporting_reads %>%
                    filter(ref_position== -1))
ins_shift[, shift_position := deduce_ref_pos_ins_shift(cigar, ref_start, forward_read_position, ref_strand, SV_start)]

INS <- sv_merged %>%
  filter(str_detect(ALT, "INS"))

INS_reads <- INS %>%
  left_join(read_mapping%>%
              dplyr::select(ID, read_id, read_alt_pos, read_strand, read_type)) %>%
  left_join(modkit_extract)

merged_INS_reads <- INS_reads %>%
  left_join(ins_shift %>%
              dplyr::select(read_id, forward_read_position, ID, shift_position)) %>%
  mutate(shift_position = case_when(ref_position != -1 & ref_position <= SV_start ~ ref_position,
                                    ref_position != -1 & ref_position > SV_start ~ ref_position+ SV_len,
                                    T ~ shift_position),
         shift_position = if_else(ref_position == -1 & ref_strand == "-", shift_position +SV_len, shift_position  ))

saveRDS(merged_INS_reads, file = paste0(sample_folder, "/", SampleID, "_INS_meth.RDS"))

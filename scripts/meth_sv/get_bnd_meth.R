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

BND <- sv_merged %>%
  filter(str_detect(ALT, "\\]")|str_detect(ALT, "\\[")) %>%
  arrange(ID) %>%
  mutate(chr_end = str_extract(ALT, "chr[0-9XYM]+"),
         SV_end =as.numeric(str_extract(ALT, "(?<=:)[0-9]+"))) %>%
  dplyr::select(ID, ALT, chr_start = `#CHROM`, SV_start, chr_end, SV_end, BND_type)

bnd_primary_mapping <- fread(paste0(sample_folder, "/bnd_read_name_primary_map.tsv"), col.names = c("read_id",  "flag", "rname", "pos", "cigar")) %>%
  mutate(strand = if_else(flag %in% c(16, 2064), "-", "+")) %>%
  rowwise() %>%
  mutate(poses = process_cigar_positions(strand, as.numeric(pos), cigar)) %>%
  separate(poses, into = c("ref_start", "ref_end", "read_start", "read_end")) %>%
  mutate(ref_start = as.numeric(ref_start), ref_end = as.numeric(ref_end),
         read_start = as.numeric(read_start), read_end = as.numeric(read_end)) %>%
  dplyr::select(-pos, -flag)


bnd_saf_mapping <- fread(paste0(sample_folder, "/bnd_read_name_saf_map.tsv"), col.names = c("read_id", "saf_tag")) %>%
  mutate(saf_tag = gsub("SA:Z:", "", saf_tag)) %>%
  separate_rows(saf_tag, sep = ";") %>%
  filter(saf_tag != "", str_detect(saf_tag, ",")) %>%
  separate(saf_tag, into = c("rname", "pos", "strand", "cigar", "MAPQ", "NM"), sep = ",") %>%
  rowwise() %>%
  mutate(poses = process_cigar_positions(strand, as.numeric(pos), cigar)) %>%
  separate(poses, into = c("ref_start", "ref_end", "read_start", "read_end")) %>%
  mutate(ref_start = as.numeric(ref_start), ref_end = as.numeric(ref_end),
         read_start = as.numeric(read_start), read_end = as.numeric(read_end)) %>%
  dplyr::select(-c(pos, MAPQ, NM))%>%
  left_join(bnd_primary_mapping %>%
              dplyr::select(read_id, ref_start_primary = ref_start, ref_end_primary = ref_end, strand , rname), 
            by= join_by(rname == rname, strand == strand, read_id == read_id, ref_start <=ref_end_primary, ref_end >= ref_start_primary)) %>%
  filter(is.na(ref_end_primary)) %>%
  dplyr::select(-ref_start_primary, -ref_end_primary)


BND_supporting_reads <- rbind(bnd_saf_mapping, bnd_primary_mapping) %>%
  left_join(modkit_extract %>%
              dplyr::select(read_id, forward_read_position, read_length, mod_qual, base_qual) %>%
              distinct(),
            by =join_by(read_id == read_id, read_start <= forward_read_position, read_end >= forward_read_position)) %>%
  dplyr::rename(chrom = rname, ref_strand = strand) %>%
  filter(!is.na(forward_read_position) ) %>%
  mutate(ref_start = if_else(ref_strand == "-", ref_end, ref_start))

BND_supporting_reads = setDT(BND_supporting_reads)
BND_supporting_reads[, ref_position := deduce_reference_position(cigar, ref_start, forward_read_position, ref_strand)]


BND_supporting_reads <-BND_supporting_reads %>%
  left_join(BND %>%
              left_join(read_mapping%>%
                          dplyr::select(ID, read_id,  read_type) %>%
                          distinct(), relationship = "many-to-many")%>%
            filter(read_type == "alt") , by = join_by(read_id == read_id, chrom ==chr_start), relationship = "many-to-many") %>%
  filter(SV_start >= ref_start, ref_end >= SV_start) %>%
  dplyr::select(-c(cigar, ref_start, ref_end, read_start, read_end)) %>%
  distinct()

BND_ref_reads  <- modkit_extract %>%
  dplyr::select(chrom, read_id, ref_strand,ref_position, forward_read_position, read_length, mod_qual, base_qual) %>%
  distinct() %>%
  right_join(BND %>%
               left_join(read_mapping%>%
                           dplyr::select(ID, read_id, read_type) %>%
                           distinct()) %>%
               filter(read_type == "ref"), by = join_by(read_id == read_id, chrom ==chr_start)) 

saveRDS(BND_ref_reads, file = paste0(sample_folder, "/", SampleID, "_BND_ref_meth.RDS"))
saveRDS(BND_supporting_reads, file = paste0(sample_folder, "/", SampleID, "_BND_alt_meth.RDS"))
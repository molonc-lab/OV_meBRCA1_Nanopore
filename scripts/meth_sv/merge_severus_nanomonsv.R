#!/usr/bin/env Rscript

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

# severus
severus_folder = paste0(wk_dir,"/output/", SampleID, "/severus")

sv_severus <- fread(paste0(severus_folder, "/somatic_SVs/severus_somatic.vcf"), skip = "#CHROM") %>%
  filter(FILTER== "PASS",
         `#CHROM` %in% paste0("chr", seq(1,22)),
          if_any(matches("phased_sorted"), ~ str_starts(.x, "0/1")),
         startsWith(INFO, "PRECISE")) %>%
  mutate(SV_type = str_extract(INFO, "(?<=SVTYPE=)[^;]+"),
         SV_end = as.numeric(str_extract(INFO, "(?<=END=)[^;]+")),
         SV_len = as.numeric(str_extract(INFO, "(?<=SVLEN=)[^;]*")),
         SV_end = if_else(SV_type == "INS", POS, SV_end),
         BND_type = case_when(startsWith(ALT, "N[chr") ~1,
                              startsWith(ALT, "]chr") ~2,
                              startsWith(ALT, "N]chr") ~3,
                              startsWith(ALT, "[chr") ~4 ),
         ID = gsub("_1$|_2$", "", ID))%>%
  dplyr::rename(SV_start = POS) 

sv_severus_germ <- fread(paste0(severus_folder, "/all_SVs/severus_all.vcf"), skip = "#CHROM")%>%
  mutate(SV_type = str_extract(INFO, "(?<=SVTYPE=)[^;]+"),
         SV_start = POS,
         SV_end = as.numeric(str_extract(INFO, "(?<=END=)[^;]+")),
         SV_len = as.numeric(str_extract(INFO, "(?<=SVLEN=)[^;]*")),
         SV_end = if_else(SV_type == "INS", POS, SV_end),
         ID = gsub("_1|_2", "", ID)) %>%
  filter(!ID %in% sv_severus,
         FILTER== "PASS",
         `#CHROM` %in% paste0("chr", seq(1,22)),
         if_any(matches("phased_sorted"), ~ str_starts(.x, "0/1")),
         startsWith(INFO, "PRECISE"))


sv_severus_VNTR <- sv_severus%>%
  filter(str_detect(INFO, "INSIDE_VNTR=TRUE"))

sv_severus<-sv_severus%>%
  filter(!str_detect(INFO, "INSIDE_VNTR=TRUE"))


#nanmonSV
nanmonsv_folder = paste0(wk_dir,"/output/", SampleID, "/SV")
sv_nanomonsv <- fread(paste0(nanmonsv_folder,"/", SampleID, ".nanomonsv.result.vcf"), skip ="#CHROM")%>%
  filter(FILTER== "PASS",
         `#CHROM` %in% paste0("chr", seq(1,22))) %>%
  separate(TUMOR, into = c("TR", "VR"), sep = ":") %>%
  mutate(SV_type = str_extract(INFO, "(?<=SVTYPE=)[^;]+"),
         SV_end = as.numeric(str_extract(INFO, "(?<=END=)[^;]+")),
         SV_len = as.numeric(str_extract(INFO, "(?<=SVLEN=)[^;]*")),
         SV_ins_len = as.numeric(str_extract(INFO, "(?<=SVINSLEN=)[^;]*")),
         SV_end = if_else(SV_type == "INS", POS, SV_end),
         BND_type = case_when(endsWith(ALT, "[") ~1,
                              startsWith(ALT, "]chr") ~2,
                              endsWith(ALT, "]") ~3,
                              startsWith(ALT, "[chr") ~4 ),
         ID = gsub("_1$|_0$", "", ID),
         TR= as.numeric(TR),
         VR = as.numeric(VR),
         SV_len = if_else(is.na(SV_len), SV_ins_len, SV_len)) %>%
  dplyr::select(-SV_ins_len) %>%
  dplyr::rename(SV_start = POS) %>%
  filter(TR <100, VR/TR <0.8, VR/TR >=0.2)


sv_nanomonsv_supporting_read <- fread(paste0(nanmonsv_folder,"/", SampleID, ".nanomonsv.supporting_read.txt")) %>%
  `colnames<-`(c("start_chr", "start_pos", "start_strand", "end_chr", "end_pos", "end_strand","ins_seq", "ID", "read_id", "alt_q_pos", "read_strand")) %>%
  filter(ID %in% sv_nanomonsv$ID) %>%
  dplyr::select("ID", "read_id", "alt_q_pos", "read_strand")
  
  
  
#merge Severus and NanomonSV

sv_merged <- sv_nanomonsv %>%
  mutate(SV_start_20u = SV_start-20, 
         SV_start_20d = SV_start +20) %>%
  left_join(sv_severus_VNTR %>%
              dplyr::select(`#CHROM`, SV_start_severus =SV_start , SV_end_severus = SV_end, SV_len_severus =SV_len, SV_type_severus = SV_type),
            by = join_by(`#CHROM`==`#CHROM`, SV_start_20u <= SV_start_severus, SV_start_20d >= SV_start_severus) ) %>%
  filter(is.na(SV_start_severus)) %>%
  dplyr::select(-ends_with("_severus")) %>%
  distinct() %>%
  left_join(sv_severus_germ%>%
              dplyr::select(`#CHROM`, SV_start_severus =SV_start ,SV_end_severus = SV_end, SV_len_severus =SV_len, SV_type_severus = SV_type),
            by = join_by(`#CHROM`==`#CHROM`, SV_start_20u <= SV_start_severus, SV_start_20d >= SV_start_severus) ) %>%
  filter(!is.na(SV_start_severus)) %>%
  dplyr::select(-c(SV_start_20u, SV_start_20d, SV_start_severus, SV_end_severus,INFO)) %>%
  distinct()

# write bed file to subset bam file 
pileup_read_id_grange<- sv_merged %>%
    mutate(SV_end = if_else(is.na(SV_end), SV_start, SV_end))%>%
    dplyr::select(seqnames = `#CHROM`, start = SV_start, end = SV_end) %>%
  distinct() %>%
  makeGRangesFromDataFrame(.)

export.bed(pileup_read_id_grange, con=paste0(wk_dir, "/output/", SampleID, "/", SampleID, "_merged_sv.bed"))

saveRDS(sv_merged, file = paste0(wk_dir, "/output/", SampleID, "/", SampleID, "_sv_merged.RDS"))
           

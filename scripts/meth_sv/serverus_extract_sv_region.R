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
  dplyr::rename(TUMOUR = 10) %>%
    filter(!str_detect(INFO, "IMPRECISE")) %>%
    dplyr::select(-FORMAT) %>%
    separate(TUMOUR, into = c("GT", "VAF", "hVAF", "DR", "DV"), sep = ":") %>%
    mutate(VAF = as.numeric(VAF), DV= as.numeric(DV)) %>%
    filter(VAF > 0.1, DV >3, VAF <1, `#CHROM` %in% paste0("chr", seq(1, 22))) %>%
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


# write bed file to subset bam file 
pileup_read_id_grange<- sv_severus %>%
    mutate(SV_end = if_else(is.na(SV_end), SV_start+1, SV_end)) %>%
    dplyr::select(seqnames = `#CHROM`, start = SV_start, end = SV_end) %>%
  distinct() %>%
  makeGRangesFromDataFrame(.)

export.bed(pileup_read_id_grange, con=paste0(wk_dir, "/output/", SampleID, "/", SampleID, "_merged_sv.bed"))

saveRDS(sv_severus, file = paste0(wk_dir, "/output/", SampleID, "/", SampleID, "_sv_merged.RDS"))

 

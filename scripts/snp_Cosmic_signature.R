args = commandArgs(trailingOnly=TRUE)

working_dir = args[1]

suppressMessages(suppressWarnings(library("data.table")))
suppressMessages(suppressWarnings(library("tidyverse")))
suppressMessages(suppressWarnings(library("deconstructSigs")))
suppressMessages(suppressWarnings(library("BSgenome.Hsapiens.UCSC.hg19")))
suppressMessages(suppressWarnings(library("BSgenome.Hsapiens.UCSC.hg38")))

#function for run deconstructsigs
run_deconstructsigs <- function(data,out_path, out_name, signatures, method = "default", plot=FALSE, cutoff=0.10) {
  # Run deconstructSigs on all samples
  for (i in 1:nrow(data)){
    sample_id=row.names(data[i,])
    sigs= whichSignatures(tumor.ref = data, 
                          signatures.ref = signatures, 
                          sample.id = sample_id,
                          signature.cutoff = cutoff,
                          contexts.needed = TRUE,
                          tri.counts.method = method)
    # Plot signatures if defined
    if ( isTRUE(plot)){
      out_file <- paste0(out_path,"figures/",sample_id,"_",method,"_",out_name,".pdf")
      pdf(out_file)
      plotSignatures(sigs)
      dev.off()
    }
    # Combine sig assignments for samples to one dataframe
    if (exists("sigs_final")) { 
      sigs_final <- rbind(sigs_final,sigs$weights) 
    } else { 
      sigs_final <- sigs$weights
    }
  }
  # Subset only detected signatures
  sigs_detected <- sigs_final %>%
    rownames_to_column(var="Sample")%>%
    gather(Sig, Count, -Sample) %>%
    group_by(Sig) %>%
    filter(Count > 0) %>%
    spread(Sig, Count, fill = 0)
  
  out_file <- paste0(out_path,"tables/",method,"_",out_name,".tsv")
  # Return signature assigment dataframe
  write.table(sigs_detected, out_file ,quote = F, sep = "\t", row.names = F)
  return(sigs_detected)
}


# TCGA 
tcga_snp_meta = Sys.glob(paste0(working_dir, "/data/raw/*_TCGA_*/MAF/*.tsv"))%>%
  map_dfr(~ as.data.frame(fread(.x)))%>%
  separate(`Sample ID`, into = c("Sample.ID_1", "Sample.ID_2"), sep = ", ") %>%
  separate(`Sample Type`, into = c("Sample.Type_1", "Sample.Type_2"), sep = ", ")%>%
  mutate(name = gsub(".gz", "", `File Name`), Sample.ID = if_else(Sample.Type_1 == "Primary Tumor", Sample.ID_1, Sample.ID_2)) %>%
  dplyr::select(name, Sample.ID)

tcga_snp <- Sys.glob(paste0(working_dir, "/data/raw/*_TCGA_*/MAF/*.maf")) %>%
  map_dfr(~ as.data.frame(fread(.x, skip = "Hugo_Symbol", select = c("Variant_Type", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"), colClasses = list(character=1:20))) %>%
            mutate(name = basename(.x))) %>%
  filter(Variant_Type == "SNP") %>%
  mutate(alt = if_else(Tumor_Seq_Allele1 == Reference_Allele, Tumor_Seq_Allele2, Tumor_Seq_Allele1), Start_Position = as.numeric(Start_Position)) %>%
  left_join(tcga_snp_meta ) %>%
  dplyr::select(Sample.ID, chr=Chromosome, pos = Start_Position, ref = Reference_Allele, alt) 

# ICGC 
icgc_snp <- Sys.glob(paste0(working_dir, "/data/raw/*_ICGC_*_meth/MAFs/*.shc.maf.gz")) %>%
  map_dfr(~ as.data.frame(fread( cmd = paste0('gunzip -cq ', .x), skip = "Hugo_Symbol", select = c("Variant_Type", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"), colClasses = list(character=1:20))) %>%
            mutate(name = basename(.x))) %>%
  filter(Variant_Type == "SNP",!str_detect(string = Chromosome ,pattern = "GL|MT")) %>%
  mutate(alt = if_else(Tumor_Seq_Allele1 == Reference_Allele, Tumor_Seq_Allele2, Tumor_Seq_Allele1), Start_Position = as.numeric(Start_Position),
         chr = paste0("chr", Chromosome)) %>%
  separate(name, into = c("name1", "name2"), extra = "drop", sep = "\\.") %>%
  unite("Sample_Name", name1:name2, sep = "_") %>%
  mutate(Sample_Name = gsub("-", "_", Sample_Name)) %>%
  dplyr::select(Sample.ID = Sample_Name, chr, pos = Start_Position, ref = Reference_Allele, alt) 



#tcga
sigs_input <- mut.to.sigs.input(mut.ref = tcga_snp,
                                sample.id = "Sample.ID",
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)

#load probability data
cosmic_t <- read_delim(paste0(working_dir,"/data/supporting/signatures_probabilities.txt"),
                       delim = "\t") 
cosmic_t %<>%
  select(1:33,-`Substitution Type`,-Trinucleotide) %>%
  arrange(match(`Somatic Mutation Type`,colnames(sigs_input))) %>%
  column_to_rownames(var = "Somatic Mutation Type") %>%
  as.matrix()
cosmic_sigs <- as.data.frame(t(cosmic_t)) %>%
  rownames_to_column(var = "Signature") %>%
  mutate(Signature = str_replace(Signature,"Signature ","Sig")) %>%
  column_to_rownames(var = "Signature")

sig_cosmic_default <- run_deconstructsigs(data = sigs_input, 
                                          out_path = paste0(working_dir, "/data/processed/Signatures/"),
                                          out_name = "cosmicv2_sig_tcga",
                                          signatures = cosmic_sigs,
                                          method = "default", 
                                          plot = T)

tcga_somatic_snp_cnts <- tcga_snp %>%
  group_by(Sample.ID) %>%
  summarise(N=n())

write_tsv(tcga_somatic_snp_cnts, paste0(working_dir, "/data/processed/Signatures/tables/tcga_somatic_snp_cnt.tsv"))

#icgc
sigs_input <- mut.to.sigs.input(mut.ref = icgc_snp,
                                sample.id = "Sample.ID",
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg19)
#load probability data
cosmic_t <- read_delim(paste0(working_dir,"/data/supporting/signatures_probabilities.txt"),
                       delim = "\t") 
cosmic_t %<>%
  select(1:33,-`Substitution Type`,-Trinucleotide) %>%
  arrange(match(`Somatic Mutation Type`,colnames(sigs_input))) %>%
  column_to_rownames(var = "Somatic Mutation Type") %>%
  as.matrix()
cosmic_sigs <- as.data.frame(t(cosmic_t)) %>%
  rownames_to_column(var = "Signature") %>%
  mutate(Signature = str_replace(Signature,"Signature ","Sig")) %>%
  column_to_rownames(var = "Signature")

sig_cosmic_default <- run_deconstructsigs(data = sigs_input, 
                                          out_path = paste0(working_dir, "/data/processed/Signatures/"),
                                          out_name = "cosmicv2_sig_icgc",
                                          signatures = cosmic_sigs,
                                          method = "default", 
                                          plot = T)

icgc_somatic_snp_cnts <- icgc_snp %>%
  group_by(Sample.ID) %>%
  summarise(N=n())

write_tsv(icgc_somatic_snp_cnts, paste0(working_dir, "/data/processed/Signatures/tables/icgc_somatic_snp_cnt.tsv"))

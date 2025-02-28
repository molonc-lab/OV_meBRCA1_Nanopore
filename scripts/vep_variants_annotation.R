args = commandArgs(trailingOnly=TRUE)


library("data.table")
library("tidyverse")
working_dir = args[1]
sample = args[2]


# out_dir = paste0(working_dir,"/output/", sample, "/deepvariants")
vcf = paste0(out_dir,"/", sample, ".vcf.gz")

mutation <- fread(vcf,skip ="#CHROM") %>%
  filter(FILTER == "PASS")


#clinvar annotation

mod_vcf <- mutation %>%
  filter(`#CHROM` %in% paste0("chr", c(seq(1,22), "X", "Y"))) %>%
  select(`#CHROM`,POS,ID,REF,ALT,QUAL,FILTER,INFO) %>%
  arrange(`#CHROM`, POS) %>%
  distinct()

fwrite(mod_vcf,
       append=TRUE,
       sep="\t",
       col.names = T,
       file=paste0(out_dir, "/vep_variants_hg38.vcf"))


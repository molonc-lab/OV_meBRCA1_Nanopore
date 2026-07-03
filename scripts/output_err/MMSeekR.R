args = commandArgs(trailingOnly=TRUE)

wk_dir = args[1]
sample = args[2]

input_path = paste0(wk_dir,"output/", sample, "/Me/", sample,"_MMSeekR.tsv")
output_path = paste0(wk_dir,"output/MMSeekR/", sample, "_MMSeekR_models")
library(MMSeekR)
runMultiModel("hg38", input_path, output_path, num.cores = 16)
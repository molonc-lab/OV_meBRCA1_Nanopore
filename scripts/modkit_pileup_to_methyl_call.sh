#!/bin/bash

echo chrBase$'\t'chr$'\t'base$'\t'strand$'\t'coverage$'\t'freqC$'\t'freqT > methyl_call_format.txt
cat *_pileup_tab.bed | awk -F "\t" '{if($4=="m") printf "%s.%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$1,$2,$6,$5,$12/$5*100,$13/$5*100}' >> methyl_call_format.txt 
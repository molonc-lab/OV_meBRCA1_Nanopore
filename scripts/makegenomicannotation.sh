module load bedtools

bedtools intersect -a genomic_positions.bed -b cpg_islands.bed -wa -u | awk '{print $0 "\tCpG Island"}' > overlapping_islands.bed
bedtools intersect -a genomic_positions.bed -b cpg_shores.bed -wa -u | awk '{print $0 "\tCpG Shore"}' > overlapping_shores.bed
bedtools intersect -a genomic_positions.bed -b cpg_shelves.bed -wa -u | awk '{print $0 "\tCpG Shelf"}' > overlapping_shelves.bed

cat overlapping_islands.bed overlapping_shores.bed overlapping_shelves.bed | sort -k1,1 -k2,2n > combined_annotations.bed

awk 'BEGIN {OFS="\t"} {
    key = $1"_"$2"_"$3;
    if (!(key in seen)) {
        seen[key] = 1;
        print $0;
    }
}' combined_annotations.bed > final_annotations.bed

import pysam
import argparse

# Parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract read IDs from a BAM file for intervals in a BED file.")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file (indexed).")
    parser.add_argument("-e", "--bed", required=True, help="Input BED file.")
    parser.add_argument("-o", "--output", required=True, help="Output file to save the read IDs.")
    return parser.parse_args()

# Main function
def main():
    # Parse arguments
    args = parse_arguments()
    bam_file = args.bam
    bed_file = args.bed
    output_file = args.output

    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Open the BED file and output file
    with open(bed_file, "r") as bed, open(output_file, "w") as output:
        for line in bed:
            # Parse BED line
            fields = line.strip().split()
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            # Fetch alignments overlapping the BED interval
            reads = bam.fetch(chrom, start, end)
            read_ids = set()  # To store unique read IDs

            for read in reads:
                read_ids.add(read.query_name)

            # Write the interval and corresponding read IDs to the output file
            if read_ids:
                output.write(f"{chrom}:{start}-{end}\t{','.join(read_ids)}\n")
            print(f"Processed interval: {chrom}:{start}-{end}, Read IDs: {len(read_ids)}")

    # Close the BAM file
    bam.close()
    print(f"Output saved to {output_file}")

if __name__ == "__main__":
    main()

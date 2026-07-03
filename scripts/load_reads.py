import pysam
import pandas as pd
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--wk_dir', dest='wk_dir', type=str, help='Add root dir')
#parser.add_argument('--wk_dir', dest='wk_dir', type=str, help='Add root dir', default="/working/lab_olgak/lijunX/OV_meBRCA1_Nanopore")
parser.add_argument('--SampleID', dest='SampleID', type=str, help='Add SampleID')
#parser.add_argument('--SampleID', dest='SampleID', type=str, help='Add SampleID', default="OVCAR8")
args = parser.parse_args()

meth_modkit = args.wk_dir+"/output/"+args.SampleID+"/Me/"+args.SampleID+"_extract_reads.tsv"

methylation_reads = pd.read_table(meth_modkit, usecols = [0, 2, 10, 11, 12])
methylation_reads = methylation_reads[methylation_reads["mod_code"] == "m"]

# load vcf files to extract deletion
vcf_file = args.wk_dir+"/output/"+args.SampleID+"/deepvariants/"+args.SampleID+".filtered.vcf.gz"
deepvariants = pysam.VariantFile(vcf_file)
deepvariantsFile = deepvariants.fetch()
variants_dict = dict()

bamfile_path=args.wk_dir+"/output/"+args.SampleID+"/"+args.SampleID+"_mapped.bam"

for variants in deepvariantsFile:
    if (variants.rlen >=5) & (variants.rlen <=50)& ("PASS" in str(variants)) & ("0/1" in str(variants)):
        variants_dict[variants.chrom + "_"+ str(variants.start) + "_" + str(variants.stop)] = dict()
        del_list = list()
        match_list = list()
        other_list = list()
        bamfile =  pysam.AlignmentFile(bamfile_path, "rb")
        iter=bamfile.fetch(variants.chrom , variants.start -20 , variants.stop + 20)
        
        #sort reads depending on whether they have the exact deletion
        for reads in iter:
           
            if reads.get_overlap(variants.start-1, variants.stop+1) == 3:
                del_list.append([reads.query_name, 3])
            elif reads.get_overlap(variants.start-1, variants.stop+1) == variants.rlen + 2:
                match_list.append([reads.query_name, reads.get_overlap(variants.start-1, variants.stop+1)])
            else:
                other_list.append([reads.query_name, reads.get_overlap(variants.start-1, variants.stop+1)])
                
        variants_dict[variants.chrom + "_"+ str(variants.start) + "_" + str(variants.stop)]["Del"] = del_list
        variants_dict[variants.chrom + "_"+ str(variants.start) + "_" + str(variants.stop)]["Match"] = match_list
        variants_dict[variants.chrom + "_"+ str(variants.start) + "_" + str(variants.stop)]["Other"] = other_list

temp_file=args.wk_dir+"/output/"+args.SampleID+"/Me/"+args.SampleID+"_variants_dictionary.pkl"

with open(temp_file, "wb") as outp:
     pickle.dump(variants_dict, outp, pickle.HIGHEST_PROTOCOL)

# with open(temp_file, "rb") as input_file:
#     variants_dict = pickle.load(input_file)
    
variants_df = pd.DataFrame.from_dict(variants_dict,
                       orient='index')

del_df = pd.DataFrame(variants_df["Del"].explode()).Del.apply(pd.Series)
del_df = del_df.reset_index()
del_df.columns = ["variant_id","read_id", "overlap"]
del_df["variant_type"] = "Del"

match_df = pd.DataFrame(variants_df["Match"].explode()).Match.apply(pd.Series)
match_df = match_df.reset_index()
match_df.columns = ["variant_id","read_id", "overlap"]
match_df["variant_type"] = "Match"

other_df = pd.DataFrame(variants_df["Other"].explode()).Other.apply(pd.Series)
other_df = other_df.reset_index()
other_df.columns = ["variant_id","read_id", "overlap"]
other_df["variant_type"] = "Other"

variants_df = pd.concat([del_df, match_df, other_df])

meth_df = variants_df.set_index("read_id").join(methylation_reads.set_index("read_id"))

methylation_variants_file=args.wk_dir+"/output/"+args.SampleID+"/Me/"+args.SampleID+"_methylation_variants.csv"
meth_df.to_csv(methylation_variants_file)
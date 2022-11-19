import sys
import csv
import re

csv.field_size_limit(sys.maxsize)

# RNA based mpileup file 
RNA_mpileup_file=sys.argv[1]

# Filtered DNA VCF file
DNA_HET_VCF=sys.argv[2]

# Sample name 
sample_name=sys.argv[3]

# Path to output directory
output_dir=sys.argv[4]


my_pileup_file=csv.reader(open(RNA_mpileup_file, "r"),delimiter="\t")
DNA_HET_file=csv.reader(open(DNA_HET_VCF,"r"),delimiter="\t")

my_pileup_file_list=list(my_pileup_file)
my_het_file_list=list(DNA_HET_file)

# Output file 
my_RNA_VAF_file=csv.writer(open(output_dir+"/"+sample_name+"_RNA_BAM_VAF_FOR_DNA_HET.txt","a"),delimiter="\t")

header=["CHR","POS","REF","ALT","Ref_count","Alt_count","DP","Alt_Allele_Freq","Ref_Allele_Freq","SAMPLE"]
my_RNA_VAF_file.writerow(header)

for i in my_het_file_list:
	# Parse DNA VCF for ref alternate positions 
	chrom=str(i[0])
	pos=int(i[2])
	ref_vcf=i[4]
	alt_vcf=i[5]
	c=0
	c1=0
	for j in my_pileup_file_list:
	#Parse mpileup file for chromosome, position, and ref and depth 
		chrom_p=str(j[0])
		pos_p=int(j[1])
		ref=j[2]
		DP=float(j[3])
		# If the chromsome, position and reference alleles from the VCF and mpileup match, then
		if (chrom==chrom_p and int(pos)==int(pos_p) and str(ref_vcf)==str(ref)):
			# calculate the counts of the alternate allele
			string=j[4]
			m1=re.findall(r'[\.|\,|a|g|t|c|A|G|T|C]',string)
			# Denominator for calculation of frequency
			denom=len(m1)
			for k in string:
				if (alt_vcf.lower()==k or alt_vcf.upper()==k):
					c=c+1
				if (k=="." or k==","):
					c1=c1+1
			# if alternate allele count is not zero, calculate the counts and frequency 
			if (int(c)!=0):
				#print ("Allelic count is not zero", int(c))
				rna_comp_vaf=round((float(c)/float(denom)),2)
				my_row=[chrom+":"+str(pos),ref_vcf,alt_vcf,int(c1),int(c),int(denom),rna_comp_vaf,abs(1.0-round(float(rna_comp_vaf),2)),sample_name]
				my_RNA_VAF_file.writerow(my_row)
				
			else:		
				if ((c==0) and int(denom)==0):
					print ("Alt allele count is 0 and the depth is also 0, mpileup string has not ref or alt characters")
					my_row=[chrom+":"+str(pos),ref_vcf,alt_vcf,c1,0,int(denom),0.0,0.0,sample_name]
					my_RNA_VAF_file.writerow(my_row)
				if ((c==0) and float(denom)!=0.0):
					print ("Alt allele count is 0 but depth is not zero, the mpileup string has some ref or alt characters")
					my_row=[chrom+":"+str(pos),ref_vcf,alt_vcf,c1,0,int(denom),0.0,round((float(c1)/float(denom)),2),sample_name]
					my_RNA_VAF_file.writerow(my_row)

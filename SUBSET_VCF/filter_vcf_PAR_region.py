import sys
import csv 

befile_chrX=csv.reader(open("REF_TRACKS/XCHR_MINUS_PAR.bed","r"),delimiter="\t")
bedfile_chrX_list=list(befile_chrX)


#### Supply VQSR filtered VCF 
vcf_file=csv.reader(open(sys.argv[1],"r",encoding = 'utf-8'),delimiter="\t")
vcf_file_list=list(vcf_file)

sample=sys.argv[2]

output_dir=sys.argv[3]
my_file=csv.writer(open(output_dir+"/"+sample+"_WES_CHRX_HETS_PAR_FILTERED.vcf","a",encoding = 'utf-8'),delimiter="\t")

header = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",sample]
my_file.writerow(header)


for i in vcf_file_list:
	if not i[0].startswith("#"):
		for j in bedfile_chrX_list:
			if (int(i[1])>=int(j[1]) and  int(i[1])<=int(j[2])):
				my_file.writerow(i)

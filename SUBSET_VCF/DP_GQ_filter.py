import sys
import csv


#### Supply the PAR filtered vcf file"
vcf_file=csv.reader(open(sys.argv[1],"r", encoding='utf-8'),delimiter="\t")
vcf_file_list=list(vcf_file)

sample=sys.argv[2]

output_dir=sys.argv[3]

my_new_file=csv.writer(open(output_dir+"/"+sample+"_DP_GQ_FILTERED.vcf","a", encoding='utf-8'),delimiter="\t")

header = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",sample]
my_new_file.writerow(header)


print (sample)
for i in vcf_file_list:
	if not i[0].startswith("#"):
		#print (i[0],i[1],i[8],i[9])
		if (i[9].split(":")[2]!="."):
			DP=int(i[9].split(":")[2])
			GQ=int(i[9].split(":")[3])
			if (DP >= 10 and GQ > 20):
				my_new_file.writerow(i)
	

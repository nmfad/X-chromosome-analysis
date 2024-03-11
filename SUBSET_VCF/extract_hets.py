import sys
import csv 


my_vcf=csv.reader(open(sys.argv[1],"r", encoding = 'utf-8'),delimiter="\t")
#print (sys.argv[1])
my_vcf_list=list(my_vcf)
sample_name=sys.argv[2]
ouput_dir=sys.argv[3]

#my_new_file=csv.writer(open(sys.argv[3]+"/"+sample_name+"_WES_ALL_HETS.vcf","a", encoding='utf-8'),delimiter="\t")
my_new_file=csv.writer(open(sys.argv[3]+"/"+sample_name+"_WES_CHRX_HETS.vcf","a", encoding='utf-8'),delimiter="\t")

header = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT", sample_name]
my_new_file.writerow(header)


for i in my_vcf_list:
	if not i[0].startswith("#") and (i[0]=="chrX" or i[0]=="X") :
		my_GT=(i[9]).split(":")[0]
		if (len(str(i[3]))==1  and len(str(i[4]))==1) and str(i[3])!="-" and str(i[4])!="-": 
			if (my_GT=="1|0" or my_GT=="0|1" or my_GT=="0/1" or my_GT=="1/0"):
				if (str(i[6])=="PASS"):
					my_new_file.writerow(i)



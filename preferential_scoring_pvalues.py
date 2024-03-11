import sys
import csv


bedfile_chrX=csv.reader(open("REF_TRACKS/CHRX_genes_merged.bed","r", encoding = 'utf-8'),delimiter="\t")
bedfile_chrX_list=list(bedfile_chrX)

pvalue_file=csv.reader(open(sys.argv[1],"r", encoding = 'utf-8'),delimiter="\t")
pvalue_file_list=list(pvalue_file)
pvalue_file_list_iter=pvalue_file_list[1:]

sample=sys.argv[2]

#header=pvalue_file_list[0]

output_dir=sys.argv[3]

myfile=csv.writer(open(sys.argv[3]+"/"+sample+"_pvalues.txt","a"),delimiter="\t")
#myfile.writerow(header)

sig_list=[]

for i in pvalue_file_list_iter:
	my_pos=i[0].split(":")[1]
	for j in bedfile_chrX_list:
		if (int(my_pos)>=int(j[1]) and  int(my_pos)<=int(j[2])):
			if (str(i[3])!="NA"):
				if (float(i[3])<=0.04):
					myrow=[i[0],i[3],j[4],i[11],i[12],sample,"SIG"]
					myfile.writerow(myrow)
					sig_list.append(i[0])
			if (str(i[3])=="NA"):
				if (str(i[7])!="NA"):
					if (float(i[7])<=0.04):
						myrow=[i[0],i[7],j[4],i[11],i[12],sample,"SIG"]
						myfile.writerow(myrow)
						sig_list.append(i[0])
				if (str(i[7])=="NA"):
					if (float(i[9])<=0.04):
						myrow=[i[0],i[9],j[4],i[11],i[12],sample,"SIG"]
						myfile.writerow(myrow)
						sig_list.append(i[0])
			
for i in pvalue_file_list_iter:
	if i[0] not in sig_list:
		my_pos=i[0].split(":")[1]
		for j in bedfile_chrX_list:
			if (int(my_pos)>=int(j[1]) and  int(my_pos)<=int(j[2])):
				if (str(i[3])=="NA"):
					if (str(i[7])!="NA"):
						myrow=[i[0],i[7],j[4],i[11],i[12],sample,"INSIG"]
						myfile.writerow(myrow)
					if (str(i[7])=="NA"):
						myrow=[i[0],i[9],j[4],i[11],i[12],sample,"INSIG"]
						myfile.writerow(myrow)
				if (str(i[3])!="NA"):
					myrow=[i[0],i[3],j[4],i[11],i[12],sample,"INSIG"]
					myfile.writerow(myrow)

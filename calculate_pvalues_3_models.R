#.libPaths(c("/usr/local/biotools/rpackages/R-4.3.2-latest",.libPaths()))
#.libPaths(c("/usr/local/biotools/r/R-4.3.2",.libPaths()))
library(gamlss)
library(arsenal)
options(max.print=1000000)
options(stringsAsFactors=FALSE)

##### This portion is doing the test 
BetaBinomialOutlierAnalysis <- function(cohort_muHat,cohort_SigmaHat,
                                        test_alt_count,test_ref_count){
  
  if(cohort_muHat <= 0){
    cohort_muHat <- .Machine$double.eps
  } 
  if(cohort_muHat >= 1){
    cohort_muHat <-  1 - .Machine$double.eps
  }
  
  # find one sided p-values
 # print (cohort_SigmaHat)
  pValueIncrease <- 1 - pBB(test_alt_count-1,bd=test_alt_count+test_ref_count,
                        mu=cohort_muHat,sigma=cohort_SigmaHat)
  #print (cohort_SigmaHat)
  pValueDecrease <- pBB(test_alt_count,bd=test_alt_count+test_ref_count,
                        mu=cohort_muHat,sigma=cohort_SigmaHat)
  
  # fix numerical errors
  pValueIncrease <- min(1,max(0,pValueIncrease))
  pValueDecrease <- min(1,max(0,pValueDecrease))
  
  # Convert to two sided test
  if( pValueDecrease< pValueIncrease){
    pValue <- min(1,2*pValueDecrease)
    direction <- "decrease"
  } else{
    pValue <- min(1,2*pValueIncrease)
    direction <- "increase"
  }
  return(list(muHat=cohort_muHat,sigmaHat=cohort_SigmaHat,
              pValue=pValue,direction))
}

##### Setting command line arguments
library("optparse")
 
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="sample VAF file", metavar="character"),
	make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="sample file name [default= %default]", metavar="character"),
	make_option(c("-o","--output_directory"), type="character", default=NULL,
	      help="Path to output directory" ,metavar="character") 
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


sample_name=opt$sample

output_dir=opt$output_directory

### First open pos_specific model file, this is the database we created 
pos_specific_model=read.csv("./REF_TRACKS/BB_params_position_and_global_specific_model_chrX_exons.txt",header=TRUE, sep="\t")
str(pos_specific_model)

###### Open gene specific model file 
gene_specific_model=read.csv("./REF_TRACKS/BB_params_gene_specific_model_chrX_exons.txt",header=TRUE, sep="\t")
str(gene_specific_model)


#### choose patient sample file and create the output dataframe with the number of rows of the input file
SL2_sample_file=read.csv(opt$file,header=TRUE, sep="\t")
str(SL2_sample_file)
SL2_sample_pvalues=as.data.frame(matrix(nrow = nrow(SL2_sample_file), ncol = 14))
colnames(SL2_sample_pvalues)=c("POS","MuHat_Pos","SigmaHat_Pos","Pvalue_pos","Direction_pos","Muhat_Gene","SigmaHat_Gene","Pvalue_Gene","Direction_Gene","Global_pvalue","Direction_Global","Ref_Count","Alt_Count","Gene")
SL2_sample_pvalues$Pvalue_pos = -1.0
SL2_sample_pvalues$Pvalue_Gene= -1.0
SL2_sample_pvalues$Global_pvalue= -1.0


for (i in seq_along(SL2_sample_pvalues$Pvalue_pos)){
  my_pos=SL2_sample_file$POS[i]
  my_gene=SL2_sample_file$GENE[i]
  test_ref_count=SL2_sample_file$Ref_count[i]+1
  test_alt_count=SL2_sample_file$Alt_count[i]+1
  newdata <- pos_specific_model[pos_specific_model$pos==my_pos,c(2,3)]
  str(newdata)
  print(nrow(newdata))
  newdata_gene <- gene_specific_model[gene_specific_model$GENE==my_gene,c(2,3)]
  print(nrow(newdata_gene))
  newdata_glob <- pos_specific_model[pos_specific_model$pos=="Global",c(2,3)]
  print (nrow(newdata_glob))
  temp=newdata_glob[1,1]
  print (temp)
  if (nrow(newdata)!= 0 & newdata[1,1] != -1 & newdata[1,2] != -1){
    print ("Entered first loop")
    cohort_muHat <- newdata[1,1]
    cohort_SigmaHat <- newdata[1,2]
    res <- BetaBinomialOutlierAnalysis(cohort_muHat,cohort_SigmaHat,test_alt_count,test_ref_count)
    SL2_sample_pvalues$POS[i]=my_pos
    SL2_sample_pvalues$MuHat_Pos[i]=res[[1]]
    SL2_sample_pvalues$SigmaHat_Pos[i]=res[[2]]
    SL2_sample_pvalues$Pvalue_pos[i]=res[[3]]
    SL2_sample_pvalues$Direction_pos[i]=res[[4]]
   }
  if (nrow(newdata)==0){
        print ("Entered second loop")
	SL2_sample_pvalues$MuHat_Pos[i]="NA"
	SL2_sample_pvalues$SigmaHat_Pos[i]="NA"
	SL2_sample_pvalues$Pvalue_pos[i]="NA"
	SL2_sample_pvalues$Direction_pos[i]="NA"
   }
  if (nrow(newdata_gene)!=0 & newdata_gene[1,1] != -1 & newdata_gene [1,2] != -1){
        print ("Entered 3rd loop for gene")
     	cohort_muHat <- newdata_gene[1,1]
     	cohort_SigmaHat <- newdata_gene [1,2]
     	res <- BetaBinomialOutlierAnalysis (cohort_muHat, cohort_SigmaHat, test_alt_count, test_ref_count)
        SL2_sample_pvalues$Muhat_Gene[i]=res[[1]]
	SL2_sample_pvalues$SigmaHat_Gene[i]=res[[2]]
	SL2_sample_pvalues$Pvalue_Gene[i]=res[[3]]
	SL2_sample_pvalues$Direction_Gene[i]=res[[4]]
   }
   if (nrow(newdata_gene)==0) {
        print ("Entered 4th loop")
	SL2_sample_pvalues$Muhat_Gene[i]="NA"
	SL2_sample_pvalues$SigmaHat_Gene[i]="NA"
	SL2_sample_pvalues$Pvalue_Gene[i]="NA"
	SL2_sample_pvalues$Direction_Gene[i]="NA"
   }
   print ("We are going global")
   cohort_muHat=newdata_glob[1,1]
   print (cohort_muHat)
   cohort_SigmaHat <- newdata_glob[1,2]
   res <- BetaBinomialOutlierAnalysis(cohort_muHat,cohort_SigmaHat,test_alt_count,test_ref_count)
   SL2_sample_pvalues$Global_pvalue[i]=res[[3]]
   SL2_sample_pvalues$Direction_Global=res[[4]]
   SL2_sample_pvalues$Ref_Count[i]=test_ref_count
   SL2_sample_pvalues$Alt_Count[i]=test_alt_count
   SL2_sample_pvalues$Gene[i]=my_gene
   SL2_sample_pvalues$POS[i]=SL2_sample_file$POS[i]
}




str(SL2_sample_pvalues)
dim(SL2_sample_pvalues)
output_file=paste(output_dir,"/",sample_name,"_chrX_p_values.txt",sep='')
write.table(SL2_sample_pvalues[order(SL2_sample_pvalues$Pvalue_pos),],file=output_file,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)




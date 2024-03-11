.libPaths(c("/usr/local/biotools/rpackages/R-3.6.2-latest",.libPaths()))
library(gamlss)
library(arsenal)
options(max.print=1000000)
options(stringsAsFactors=FALSE)

estWithGAMLSS <- function(df){
  # fit parameters accross all individuals in cohort
  # notice we allow the total observations to be affected
  # by proportion of successes
  fit <- gamlss(cbind(successes, failures) ~ 1,
                data = df,family=BB)
  pred <- predictAll(fit,data = df)
  
  return(list(muHat=pred$mu[1],sigmaHat=pred$sigma[1]))
}


#### This will be the input file with the genes 
gtex_file<-read.csv("/research/bsi/projects/PI/tertiary/Klee_Eric_mrl2075/s212354.RadiaNT/ASEAnalysis/XSKEW_PROJECT_WITH_ADDN_SAMPLES_NEW_FILTERS/parameter_estimates_chrX/GTEX_BINOM_BIG_TABLE_DP_grt_3_exons.txt",header=TRUE, sep="\t")
str(gtex_file)
gene_mod_pos<-data.frame(unique(gtex_file$GENE))
colnames(gene_mod_pos) <- "GENE"
str(gene_mod_pos)

# do gene specific modeling , the lines below are 2 columns being 
# assigned a value of 0.0 
gene_mod_pos$muHat <- -1.0
gene_mod_pos$SigmaHat <- -1.0
for (i in seq_along(gene_mod_pos$muHat)){
  curr_gene=gene_mod_pos$GENE[i]
  newdata <- gtex_file[gtex_file$GENE==curr_gene,c(4,5)] + 1
  colnames(newdata) <- c("failures","successes")
  if (nrow(newdata)>=10){
    print (curr_gene)
    print (i)
    params = tryCatch({
      params<-estWithGAMLSS(newdata)
      gene_mod_pos$muHat[i] <- params$muHat
      gene_mod_pos$SigmaHat[i] <- params$sigmaHat
    },error=function(err){
      gene_mod_pos$muHat[i] <- "zero_inflated"
      gene_mod_pos$SigmaHat[i] <- "zero_inflated"
      print (paste("MY_ERROR: This position is giving problems",err))
    })
    
  }else{
    gene_mod_pos$muHat[i] <- NA
    gene_mod_pos$SigmaHat[i] <- NA
  }
}

gene_mod_pos <- gene_mod_pos[!is.na(gene_mod_pos$muHat),]
write.table(gene_mod_pos,file="/research/bsi/projects/PI/tertiary/Klee_Eric_mrl2075/s212354.RadiaNT/ASEAnalysis/XSKEW_PROJECT_WITH_ADDN_SAMPLES_NEW_FILTERS/parameter_estimates_chrX/BB_params_gene_specific_model_chrX_exons.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


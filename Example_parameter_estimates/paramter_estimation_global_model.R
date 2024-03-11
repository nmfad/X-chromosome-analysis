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


#### This will be the main big gtex file we used for position specific modelling 
gtex_file<-read.csv("/research/bsi/projects/PI/tertiary/Klee_Eric_mrl2075/s212354.RadiaNT/ASEAnalysis/XSKEW_PROJECT_WITH_ADDN_SAMPLES_NEW_FILTERS/parameter_estimates_chrX/GTEX_BINOM_BIG_TABLE_DP_grt_3_exons.txt",header=TRUE, sep="\t")
str(gtex_file)
global_pos_mod_df <- dplyr::sample_n(gtex_file, 2000)
str(global_pos_mod_df)
global_pos_mod <- data.frame(matrix(nrow=1,ncol=3))
colnames(global_pos_mod)=c("POS","muHat","SigmaHat")
#global_mod_pos<-data.frame(unique(gtex_file$CHR.POS))
#colnames(gene_mod_pos) <- "POS"
#str(gene_mod_pos)

# do global modeling , the lines below are 2 columns being 
global_pos_mod$muHat <- -1.0
global_pos_mod$SigmaHat <- -1.0
counts_df=global_pos_mod_df[,c(4,5)]+1
colnames(counts_df) <- c("failures","successes")
params = tryCatch({
      params <- estWithGAMLSS(counts_df)
      global_pos_mod$muHat[1] <- params$muHat
      global_pos_mod$SigmaHat[1] <- params$sigmaHat
    },error=function(err){
      global_pos_mod$muHat[1] <- "zero_inflated"
      global_pos_mod$SigmaHat[1] <- "zero_inflated"
      print (paste("MY_ERROR: This position is giving problems",err))
    })


global_pos_mod$POS[1]="Global"
str(global_pos_mod)
write.table(global_pos_mod,file="/research/bsi/projects/PI/tertiary/Klee_Eric_mrl2075/s212354.RadiaNT/ASEAnalysis/XSKEW_PROJECT_WITH_ADDN_SAMPLES_NEW_FILTERS/parameter_estimates_chrX/BB_params_global_position_model_new.txt",sep="\t",row.names=TRUE,col.names=c("POS","muHat","SigmaHat"))


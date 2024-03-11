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

gtex_file<-read.csv("/research/bsi/projects/PI/tertiary/Klee_Eric_mrl2075/s212354.RadiaNT/ASEAnalysis/XSKEW_PROJECT_WITH_ADDN_SAMPLES_NEW_FILTERS/parameter_estimates_chrX/GTEX_BINOM_BIG_TABLE_DP_grt_3_exons.txt",header=TRUE, sep="\t")
str(gtex_file)
gen_mod_pos<-data.frame(unique(gtex_file$CHR.POS))
colnames(gen_mod_pos) <- "POS"
str(gen_mod_pos)

# do position specific modeling , the lines below are 2 columns being 
# assigned a value of 0.0 
gen_mod_pos$muHat <- -1.0
gen_mod_pos$SigmaHat <- -1.0
for (i in seq_along(gen_mod_pos$muHat)){
  pos=gen_mod_pos$POS[i]
  newdata <- gtex_file[gtex_file$CHR.POS==pos,c(4,5)] + 1
  colnames(newdata) <- c("failures","successes")
  if (nrow(newdata)>=10){
    print (pos)
    print (i)
    params = tryCatch({
      params<-estWithGAMLSS(newdata)
      gen_mod_pos$muHat[i] <- params$muHat
      gen_mod_pos$SigmaHat[i] <- params$sigmaHat
    },error=function(err){
      gen_mod_pos$muHat[i] <- "zero_inflated"
      gen_mod_pos$SigmaHat[i] <- "zero_inflated"
      print (paste("MY_ERROR: This position is giving problems",err))
    })
    
  }else{
    gen_mod_pos$muHat[i] <- NA
    gen_mod_pos$SigmaHat[i] <- NA
  }
}
gen_mod_pos <- gen_mod_pos[!is.na(gen_mod_pos$muHat),]
write.table(gen_mod_pos,file="/research/bsi/projects/PI/tertiary/Klee_Eric_mrl2075/s212354.RadiaNT/ASEAnalysis/XSKEW_PROJECT_WITH_ADDN_SAMPLES_NEW_FILTERS/parameter_estimates_chrX/BB_params_position_specific_model_chrX_exons.txt",sep="\t",row.names=TRUE,col.names=c("pos","muHat","SigmaHat"))


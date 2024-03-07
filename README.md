# Identification of skewed X chromosome inactivation (XCI) using exome and transcriptome sequencing in patients with suspected rare genetic disease

BACKGROUND 

X-chromosome inactivation (XCI) is an epigenetic process that occurs during early development in mammalian females by randomly silencing one of two copies of the X chromosome in each cell. The preferential inactivation of either the maternal or paternal copy of the X chromosome in a majority of cells results in a skewed or non-random pattern of X inactivation and is observed in over 25% of adult females. Identifying skewed X inactivation is of clinical significance in patients with suspected rare genetic diseases due to the possibility of biased expression of disease-causing genes present on the active X chromosome. The current clinical test for the detection of skewed XCI relies on the methylation status of the methylation-sensitive restriction enzyme (Hpall) binding site present in proximity of short tandem polymorphic repeats on the androgen receptor (AR) gene. This approach using one locus results in an uninformative or inconclusive data for 10-20 % of tests. Further, recent studies have shown inconsistency between methylation of the AR locus and the state of inactivation of the X chromosome. The proposed method uses exome and transcriptome sequencing data integrated with the concept of outlier analysis for generating a chromosome wide measure of skew along the X chromosome. The method first uses a reference model for skewed XCI generated using DNA and RNA-seq data for 135 females from the GTEX reference consortium. A cohort of patient samples (or an individual sample) is then tested against the generated model to test the null hypothesis that the allelic counts mined from the patient's transcriptome using high conidence heterozygous positions from genotyping data are derived from the same distribution as the allelic counts for the reference population. A significant p-value (< 0.05) counts evaluates a position to be significantly skewed. The percentage of significantly skewed positions (>12%) classifies a sample to present with skewed or random XCI. 

REQUIREMENTS:
Python 2.7.10 or higher
R-3.6.2
Samtools
Bedtools
genome build - hg19 

STEPS:
(1) The method relies on an X-skew model generated using 135 females from the GTEX reference consortium. The model is generated using positions within the reference cohort, genes within the reference cohort and also includes a global model which is a random sample of 2000 positions present within the general population. We described these models as 
  (a) Position specific model 
  (b) gene specific model
  (c) Global model 
  These models exist as tab delimited files located within the GTEX_PARAMETER_ESTIMATES directory. These have been generated using 135 GTEX reference females (v7) on hg19. 

(2) Data Pre-processing steps
    (a) Generation of 


    (b) RNA BAM 

     

     



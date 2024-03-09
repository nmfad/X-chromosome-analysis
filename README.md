# Identification of skewed X chromosome inactivation (XCI) using exome and transcriptome sequencing in patients with suspected rare genetic disease

BACKGROUND 

X-chromosome inactivation (XCI) is an epigenetic process that occurs during early development in mammalian females by randomly silencing one of two copies of the X chromosome in each cell. The preferential inactivation of either the maternal or paternal copy of the X chromosome in a majority of cells results in a skewed or non-random pattern of X inactivation and is observed in over 25% of adult females. Identifying skewed X inactivation is of clinical significance in patients with suspected rare genetic diseases due to the possibility of biased expression of disease-causing genes present on the active X chromosome. The proposed method uses exome and transcriptome sequencing data integrated with the concept of outlier analysis for generating a chromosome wide measure of skew along the X chromosome. The method first uses a reference model for skewed XCI generated using DNA and RNA-seq data for 135 females from the GTEX reference consortium. A cohort of patient samples (or an individual sample) is then tested against the generated model to test the null hypothesis that the allelic counts mined from the patient's transcriptome using high conidence heterozygous positions from genotyping data are derived from the same distribution as the allelic counts for the reference population. A significant p-value (< 0.05) counts evaluates a position to be significantly skewed. The percentage of significantly skewed positions (>12%) classifies a sample to present with skewed or random XCI. 

REQUIREMENTS:
Python 2.7.10 or higher

R-3.6.2

Samtools

Bedtools

genome build - hg19 

STEPS:

(1) The method relies on an X-skew model generated using 135 females from the GTEX reference consortium. The model is generated using positions within the reference cohort, genes within the reference cohort and also includes a global model which is a random sample of 2000 positions present within the general population. We described these models as 

  (a) Position specific model - GTEX_PARAMETER_ESTIMATES/BB_params_position_and_global_specific_model_chrX_exons.txt
  
  (b) gene specific model - GTEX_PARAMETER_ESTIMATES/BB_params_gene_specific_model_chrX_exons.txt
  
  (c) Global model - Same file as position specific model ; GTEX_PARAMETER_ESTIMATES/BB_params_position_and_global_specific_model_chrX_exons.txt
  
  These models have been pre-generated for users and exist as tab delimited files located within the GTEX_PARAMETER_ESTIMATES directory. They need not be regenerated. The models were generated using 135 GTEX reference females (v7) using human genome reference build hg19. If users would like to generate their own reference models using their choice of reference cohort ; please see user guide on wiki page for details on generation a reference model. The online wiki page provides step by step examples for creation of models using the GTEx reference cohorts but can be applied to any cohorts. 

(2) Data Pre-processing steps for samples undergoing evaluation for X-skew testing. 

    (a) Generation of variant calls from DNA - Users may use any DNA aligner and variant caller of choice. We have used an updated version of the TREAT worfklow [22088845]. 
         (1) Next, Subset the VCF files to include calls that lie on the X chromosome, genotyped as HET (variants genotyped as 1/0 or 0/1 within the GT field) and are "PASS" in the "QUAL" field. Users can do this either using GATK SelectVariants or a utility script for this has been provided within the SUBSET_VCF directory - extract_hets.py
         Usage : python extract_hets.py <input file> <sample name> <output directory> 

         (2) The VCF file should exclude PAR regions. These regions are difficult to genotype owing to high homology and thus are excluded from X-skew analysis. A bedfile for PAR regions corresponding to human genome reference build hg19 is provided within the REF_TRACKS directory (XCHR_MINUS_PAR.bed) . Filter the VCF file to exclude these regions. Users can do this using a tool of thier choice such as Bedtools or a utility script is provided within the SUBSET_VCF directory - filter_vcf_PAR_region.py
         Usage : python filter_vcf_PAR_region.py <vcf file from step (1)> <sample name> <output directory>

         (3) The VCF file should be filtered to a depth of (DP>=10) and genotype quality (GQ > 20). Users can use any filtering tool of their choice such as GATK SelectVariants to achieve this or a utility script is provided for the same within the SUBSET_VCF directory - DP_GQ_filter.py
          Usage: python DP_GQ_filter.py <input vcf> <sample name> <output directory> 

         (4) 


    (b) RNA BAM 
filter_vcf_PAR_region.py
     

     



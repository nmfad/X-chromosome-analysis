# Identification of skewed X chromosome inactivation (XCI) using exome and transcriptome sequencing in patients with suspected rare genetic disease

BACKGROUND 

X-chromosome inactivation (XCI) is an epigenetic process that occurs during early development in mammalian females by randomly silencing one of two copies of the X chromosome in each cell. The preferential inactivation of either the maternal or paternal copy of the X chromosome in a majority of cells results in a skewed or non-random pattern of X inactivation and is observed in over 25% of adult females. Identifying skewed X inactivation is of clinical significance in patients with suspected rare genetic diseases due to the possibility of biased expression of disease-causing genes present on the active X chromosome. The proposed method uses exome and transcriptome sequencing data integrated with the concept of outlier analysis for generating a chromosome wide measure of skew along the X chromosome. The method first uses a reference model for skewed XCI generated using DNA and RNA-seq data for 135 females from the GTEX reference consortium. A cohort of patient samples (or an individual sample) is then tested against the generated model to test the null hypothesis that the allelic counts mined from the patient's transcriptome using high conidence heterozygous positions from genotyping data are derived from the same distribution as the allelic counts for the reference population. A significant p-value (< 0.05) counts evaluates a position to be significantly skewed. The percentage of significantly skewed positions (>12%) classifies a sample to present with skewed or random XCI. 

REQUIREMENTS:

1. Python 2.7.10 or higher

2. R-3.6.2 or higher 

3. Samtools

4. Bedtools

5. genome build - hg19 

STEPS:

(1) The method relies on an X-skew model generated using 135 females from the GTEX reference consortium. The model is generated using positions within the reference cohort, genes within the reference cohort and also includes a global model which is a random sample of 2000 positions present within the general population. These models have been pre-generated for users and exist as tab delimited files located within the GTEX_PARAMETER_ESTIMATES directory. They need not be regenerated. The models were generated using 135 GTEX reference females (v7) using human genome reference build hg19. We described these models as follows and are located within the GTEX_PARAMETER_ESTIMATES directory

  (a) Position specific model - GTEX_PARAMETER_ESTIMATES/BB_params_position_and_global_specific_model_chrX_exons.txt
  
  (b) gene specific model - GTEX_PARAMETER_ESTIMATES/BB_params_gene_specific_model_chrX_exons.txt
  
  (c) Global model -  GTEX_PARAMETER_ESTIMATES/BB_params_position_and_global_specific_model_chrX_exons.txt (Same file as position specific model file) 
  
   If users would like to generate their own reference models using their choice of reference cohort ; please see user guide on wiki page for details on generation a reference model. The online wiki page provides step by step examples for creation of models using the GTEx reference cohorts but can be applied to any cohorts. 

(2) Data Pre-processing steps for samples undergoing evaluation for X-skew testing. 

    (a) Generation of variant calls from DNA - Users may use any DNA aligner and variant caller of choice. We have used an updated version of the TREAT worfklow [22088845]. 
         (1) Next, Subset the VCF files to include SNV calls that lie on the X chromosome, genotyped as HET (variants genotyped as 1/0 or 0/1 within the GT field) and are "PASS" in the "QUAL" field. Users can do this either using GATK SelectVariants or a utility script for this has been provided within the SUBSET_VCF directory - extract_hets.py
         Usage : python extract_hets.py <input file> <sample name> <output directory> 

         (2) The VCF file should exclude PAR regions. These regions are difficult to genotype owing to high homology and thus are excluded from X-skew analysis. A bedfile for PAR regions corresponding to human genome reference build hg19 is provided within the REF_TRACKS directory (XCHR_MINUS_PAR.bed) . Filter the VCF file to exclude these regions. Users can do this using a tool of thier choice such as Bedtools or a utility script is provided within the SUBSET_VCF directory - filter_vcf_PAR_region.py
         Usage : python filter_vcf_PAR_region.py <vcf file from step (1)> <sample name> <output directory>

         (3) The VCF file should be filtered to a depth of (DP>=10) and genotype quality (GQ > 20). Users can use any filtering tool of their choice such as GATK SelectVariants to achieve this or a utility script is provided for the same within the SUBSET_VCF directory - DP_GQ_filter.py
          Usage: python DP_GQ_filter.py <input vcf> <sample name> <output directory> 

         (4) Subset the filtered file from (1) , (2) and (3) to only exon/coding regions on the X chromosome. A bedfile for all exon only regions on the X chromosome is provided within the REF_TRACKS directory REF_TRACKS/XCHR_MINUS_PAR.bed

         (5) Using the VCF file generated in step 4 , prepare a bedfile that will be used for generating mpileups using RNA bams. Ensure the bedfile in a valid format from the VCF file such as
            chrX  pos-1  pos


    (b) RNA BAM 
        (1) Users may use any splice-aware aligner of choice. This study makes use of the MAPR-seq pipeline [24972667] that uses the Tophat aligner (https://ccb.jhu.edu/software/tophat/index.shtml) and generates sorted aligned files compatible with the BAM (Binary Alignment Map) format.  (Note: Duplicate reads are not marked or removed). 

        (2) Using the bedfile generated in steps(a)(5) above using high confidence heterozygous positions from the DNA seqeuncing data for a sample, and the alignment file from RNA-seq data, mpileups are generated to enable generating counts of alt and ref alleles within the RNA-seq data. Samtools v1.9 or higher was used for generation of mpileups. Example command given below 
        samtools mpileup -A -d 8000 -q 20 -Q 30 -f <ref.fa> -l <sample_DNA_filtered_het_calls.bed> <same_RNA.bam> -o <sample_mpileup.txt>

        (3) Generate computed variant allele counts (CVACs) using mpileup data
            Usage: python CVAC_calculation.py <sample_mpileup.txt> <filtered DNA HET VCF from step (a)(4)> <sample name> <output directory>
            Output file generated by script : <sample_name>_RNA_BAM_VAF_FOR_DNA_HET.txt
            Important: To the output generated above, add a final column for gene name called "GENE" by overlapping the output with REF_TRACKS/chrX_genes_exons_transcripts_sorted.bed. Name the final file such as <samplename>_RNA_BAM_VAF_FOR_DNA_HET_sorted_gene.txt
            The format of the final file with the "GENE" column included should be tabdelimited and as follows:
            POS	REF	ALT	Ref_count	Alt_count	DP	Alt_Allele_Freq	Ref_Allele_Freq	SAMPLE	GENE		
	    chrX:100  T	A	0	30	30	1.0	0.0	<sample>	TRMT2B
             
            Column definitions:
              POS - chrX:100  - high confidence HET position derived from the DNA seq data of the indivdiual 
              REF - "T" - reference allele of the high confidence HET position derived from the DNA seq data of the individual
              ALT - "A" - alternate allele of the high confidence HET position derived from the DNA seq data of the individual
	      Ref_count - Count of the reference allele using the RNA sequencing reads derived from the RNA seq data of the individual
              Alt_count - Count of the alternate allele using the RNA sequencing reads derived from the RNA seq data of the individual
	      DP - Positional depth from RNA sequencing reads derived from the RNA seq data of the individual
              Alt_Allele_Freq - Fraction of Alt_count divided by the DP 
	      Ref_Allele_Freq - Fraction of the Ref_count divided by the DP
       	      Sample - SampleID or name of the sample 
	      GENE - Gene in which the position lies. If the position lies within overlapping genes, you can include both the genes separated by a comma. 

(3) Outlier based analysis - The CVAC file generated using DNA and RNA sequencing data in step (b)(3) is then used to query the reference model or GTEX model using the script calculate_pvalues_3_models.R 
	Usage: Rscript calculate_pvalues_3_models.R -f <samplename_RNA_BAM_VAF_FOR_DNA_HET_sorted_gene.txt> -s <samplID> -o <output directory>
  	Output: <sampleID>_chrX_p_values.txt
        Column definition
	POS - position queried from the input file
        MuHat_Pos - parameter estimate 1 for the position from the GTEX reference cohort
        SigmaHat_Pos - parameter estimate 2 for the position from the GTEX reference cohort 
        Pvalue_pos - Probability (P-value) based on parameter estimates (1 and 2) from the GTEX reference cohort 
        Direction_pos - Direction of the allelic bias observed (increase - allelic bias observed is in the direction of the alternate allele)
        Muhat_Gene - parameter estimate 1 for the gene being queried from the GTEX reference cohort 
        SigmaHat_Gene - paramerter estimate 2 for the gene being queried from the GTEX reference cohort
        Pvalue_Gene - Probability (P-value) based on parameter estimates (1 and 2) from the GTEX reference cohort for the gene being queried 
        Direction_Gene - Direction of the allelic bias observed (increase - allelic bias observed is in the direction of the alternate allele)
        Global_pvalue - Probability (P-value) based on parameter estimates (1 and 2) from the GTEX reference cohort  
        Direction_Global
        Ref_Count
        Alt_Count
        Gene

(4) Calculation of X-skew 
     Using the file generated in step (3) above <sampleID>_chrX_pvalues.txt
     (a) Count the total number of positions at DP >= 10 
     (b) Count the total number of positions at DP >=10 and those that are annotated as "SIG" meaning they had significant p-values 
     (c) Calculate the percentage of significantly skewed positions by doing ((b)/(a))*100.0 
     (d) If the percentage of significantly skewed positions is greater than 12%, then the sample is classified as a skewed sample or else it is classified as a random sample. 
       




    

     



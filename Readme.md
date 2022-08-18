# Calling mutations 

This code runs basic QC for tumor & normal bams, calls somatic mutations with strelks, filters and annotates with VEP, adds custom filters and gets segments for tumor samples and annotates mutations according to CNA status. 

# Requirements

Requires R, Python3.8, bash and set of various tools. 

## Python3.8 

See requirements.txt 

## R 

### Facets

Facets to get CNA, purity, ploidy. 
Facets library and all of its dependencies as in https://github.com/mskcc/facets/

## Tools 

### Fastqc 

Fastqc tool to check bam quality. 
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
Install as in instruction

### Mosdepth 

Mosdepth to check coverage. 
Mosdepth binary can be downloaded for latest realeases. https://github.com/brentp/mosdepth/releases 

### Snp-pileup

Snp-pileup code required for creating pileups for QC and for facets tool. 
Snp-pilep script from https://github.com/mskcc/facets/tree/master/inst/extcode (clone repo and compile snp-pileup as instructed)

### Strelka

Strelka somatic caller. 
Download from releases or compiled for code. https://github.com/Illumina/strelka

### Manta 

Indel and structural variation tool for strelka best performance. 
Download from releases or compiled for code.  https://github.com/Illumina/manta

### VEP

Ensembl VEP annotator for vcfs. 
Download and install command line tool as inctructed. Installing cashe for hg38 (not merged) is required. http://www.ensembl.org/info/docs/tools/vep/index.html 

# Input 

Requires tumor and normal sorted bams, reference sequence, reference dbsnp vcf file, regions of sequencing as input (see how to input below in Running)

# Output 

Only endfiles are listed here, intermediaries are not listed.   
  
SAMPLE_NAME_fastqc.html - fastqc reports  
normal-mosdepth.mosdepth.summary.txt - normal coverage  
tumor-mosdepth.mosdepth.summary.txt tumor - tumor coverage  
tumor-normal-correlation.txt - correlation between tumor&normal  
segments.txt - facets CNA segments  
facets_plot.pdf - segments plotted  
purity_ploidy.txt - sample purity and ploidy by facets  
strelka_somatic_merged_unfiltered.vcf.gz - unfiltered and not annotated full vcf  
strelka_somatic_merged_filtered_postprocessed.vcf - filtered (Strelka FILTER == PASS), VEP-annotated and custom annotated vcf.  

# Preparation 

To prepare for running, enter tools locations in tool_config_template.sh and save as tool_config.sh 


# Running 

```sh
full_process.sh PATIENT_DIRECTORY REFERENCE.FASTA REFERENCE.VCF REGIONS.BED TUMOR.BAM NORMAL.BAM PREPARATION
```
PATIENT_DIRECTORY - link to patient directory, where all files will be created  
REFERENCE.FASTA - link to reference sequence  
REFERENCE.VCF - link to reference vcf  
REGIONS.BED - link to sequencing regions bed file  
TUMOR.BAM - link to tumor bam  
NORMAL.BAM - link to normal bam  
PREPARATION - preparation type, FFPE for parafine, FF - for fresh/fresh frozen  

Or can be run step by step   

```sh
process_patient.sh PATIENT_DIRECTORY REFERENCE.FASTA REFERENCE.VCF REGIONS.BED TUMOR.BAM NORMAL.BAM  
prepare_n_annotate.sh PATIENT_DIRECTORY REFERENCE.FASTA
postprocess.sh PATIENT_DIRECTORY PREPARATION
```
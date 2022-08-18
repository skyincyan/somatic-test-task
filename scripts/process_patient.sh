#!/bin/bash
# Runs QC and calls somatic mutations from tumor 

set -o errexit
set -o pipefail
set -o nounset

PATIENT_DIR=$1
REFERENCE=$2
REFERENCE_VCF=$3
REGIONS=$4
TUMOR=$5
NORMAL=$6
GIT_DIR=$7

source "$GIT_DIR"/tool_config.sh

# quick qc 

"$FASTQC_DIR"/fastqc "$TUMOR" "$NORMAL" -o "$PATIENT_DIR"

touch "$TUMOR".bai
touch "$NORMAL".bai

"$MOSDEPTH_DIR"/mosdepth --by "$REGIONS" "$PATIENT_DIR"/tumor-mosdepth "$TUMOR" 
"$MOSDEPTH_DIR"/mosdepth --by "$REGIONS" "$PATIENT_DIR"/normal-mosdepth "$NORMAL"

"$PILEUP_DIR"/snp-pileup -g -q15 -Q20 -P100 -r25,0 "$REFERENCE_VCF" "$PATIENT_DIR"/normal-pileup "$NORMAL" 
"$PILEUP_DIR"/snp-pileup -g -q15 -Q20 -P100 -r25,0 "$REFERENCE_VCF" "$PATIENT_DIR"/tumor-pileup "$TUMOR"

# calling variants

# manta for better indels

"$MANTA_DIR"/bin/configManta.py --exome \
--normalBam "$NORMAL" \
--tumorBam "$TUMOR" \
--referenceFasta "$REFERENCE" \
--callRegions "$REGIONS" \
--runDir "$PATIENT_DIR"/manta_dir \

"$PATIENT_DIR"/manta_dir/runWorkflow.py -j 2

# strelka calling

"$STRELKA_DIR"/bin/configureStrelkaSomaticWorkflow.py --exome \
--normalBam "$NORMAL" \
--tumorBam "$TUMOR" \
--referenceFasta "$REFERENCE" \
--callRegions "$REGIONS" \
--runDir "$PATIENT_DIR"/strelka_dir \
--indelCandidates "$PATIENT_DIR"/manta_dir/results/variants/candidateSmallIndels.vcf.gz

"$PATIENT_DIR"/strelka_dir/runWorkflow.py -m local -j 3 -g 4

# merge vcfs 

"$BCFTOOLS_DIR"/bcftools merge --merge all "$PATIENT_DIR"/strelka_dir/results/variants/somatic.indels.vcf.gz "$PATIENT_DIR"/strelka_dir/results/variants/somatic.snvs.vcf.gz --force-samples -O z -o "$PATIENT_DIR"/merged_vcf.vcf.gz

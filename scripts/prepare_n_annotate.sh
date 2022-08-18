# !/bin/bash
# Prepares for VEP and runs VEP 

set -o errexit
set -o pipefail
set -o nounset

PATIENT_DIR=$1
REFERENCE=$2

source tool_config.sh 

python3.8 $UTILITY_DIR/clear_merge_n_filter.py $PATIENT_DIR

$VEP_DIR/vep -i $PATIENT_DIR/strelka_somatic_merged_filtered.vcf.gz \
-o $PATIENT_DIR/strelka_somatic_merged_filtered_annotated.vcf.gz \
--format vcf --vcf --species homo_sapiens -a GRCh38 --no_stats --fork 2 --cache --fasta $REFERENCE \
--sift b --polyphen b --ccds --hgvs --symbol --regulatory --canonical --protein --biotype  --af --af_1kg --af_gnomad \
--max_af --uniprot --tsl --appris --gene_phenotype --pubmed --domains --numbers --no_escape \
--domains --numbers --no_escape --compress_output gzip
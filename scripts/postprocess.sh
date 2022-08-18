# !/bin/bash
# Processes pileups and calls facets (used in annotation), then adds custom filter and annotations to vcf 

set -o errexit
set -o pipefail
set -o nounset

PATIENT_DIR=$1
PREPARATION=$2

source tool_config.sh 

# process pileups, get snp correlation and prepare pileups for facets
python3.8 $UTILITY_DIR/pileups_process.py $PATIENT_DIR

# run facets for segments and purity-ploidy
Rscript $UTILITY_DIR/run_facets.R $PATIENT_DIR

# custom filters for vcf and total/minor annotation for mutations 
python3.8 $UTILITY_DIR/postprocess_vcf.py $PATIENT_DIR $PREPARATION
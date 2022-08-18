#!/bin/bash
# All patient processing. See Readme.md for instruction

set -o errexit
set -o pipefail
set -o nounset

PATIENT_DIR=$1
REFERENCE=$2
REFERENCE_VCF=$3
REGIONS=$4
TUMOR=$5
NORMAL=$6
PREPARATION=$7

# run everything 

./process_patient.sh "$PATIENT_DIR" "$REFERENCE" "$REFERENCE_VCF" "$REGIONS" "$TUMOR" "$NORMAL"  
./prepare_n_annotate.sh "$PATIENT_DIR" "$REFERENCE"
./postprocess.sh "$PATIENT_DIR" "$PREPARATION"
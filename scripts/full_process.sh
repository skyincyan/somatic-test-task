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
GIT_DIR=$8

# run everything 

"$GIT_DIR"/process_patient.sh "$PATIENT_DIR" "$REFERENCE" "$REFERENCE_VCF" "$REGIONS" "$TUMOR" "$NORMAL" "$GIT_DIR"
"$GIT_DIR"/prepare_n_annotate.sh "$PATIENT_DIR" "$REFERENCE" "$GIT_DIR"
"$GIT_DIR"/postprocess.sh "$PATIENT_DIR" "$PREPARATION" "$GIT_DIR"
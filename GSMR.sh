#!/bin/bash

# ==============================================================================
# Generalized Summary-data-based Mendelian Randomization Analysis
#
# Description:
#   This script runs GCTA-GSMR for all combinations of exposure and outcome
#   summary-statistics files.
#
# ==============================================================================

set -euo pipefail

# ============================================
# Run GSMR for exposure-outcome pairs
# ============================================

GCTA="<SOFT_DIR>/gcta64"
BFILE="<LDREF_DIR>/ld_reference"
EXPOSURE_DIR="<PROJECT_DIR>/gsmr_input/exposure"
OUTCOME_DIR="<PROJECT_DIR>/gsmr_input/outcome"
RESULT_DIR="<PROJECT_DIR>/gsmr_results"
NCPU=20

mkdir -p "${RESULT_DIR}"

#--------------------------------------------------
# Run GSMR for one exposure-outcome pair
#--------------------------------------------------
run_gsmr() {
    exp_file=$1
    out_file=$2

    exp_name=$(basename "$exp_file" _gwas.cojo.txt)
    out_name=$(basename "$out_file" .cojo.txt)

    exp_dir=$RESULT_DIR/$exp_name
    tmp_dir=$exp_dir/tmp
    mkdir -p "$tmp_dir"

    exp_list=$tmp_dir/${exp_name}.exposure.txt
    out_list=$tmp_dir/${out_name}.outcome.txt

    # Prepare exposure and outcome file lists required by GCTA-GSMR
    echo -e "${exp_name}\t${exp_file}" > "$exp_list"
    echo -e "${out_name}\t${out_file}" > "$out_list"

    out_prefix=$exp_dir/${exp_name}_to_${out_name}

    $GCTA \
      --bfile "$BFILE" \
      --gsmr-file "$exp_list" "$out_list" \
      --gsmr-direction 0 \
      --gsmr-snp-min 1 \
      --gwas-thresh 5e-8 \
      --out "$out_prefix"
}

export -f run_gsmr
export GCTA BFILE RESULT_DIR

#--------------------------------------------------
# Run all exposure-outcome pairs
#--------------------------------------------------
parallel -j $NCPU run_gsmr ::: \
    $EXPOSURE_DIR/*_gwas.cojo.txt ::: \
    $OUTCOME_DIR/*.cojo.txt
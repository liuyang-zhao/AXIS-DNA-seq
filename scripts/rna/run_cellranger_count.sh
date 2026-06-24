#!/bin/bash
# AXIS spatial RNA Cell Ranger pipeline.
#
# Converts raw RNA FASTQ (R1=cDNA, R2=spatial barcode) into Cell Ranger layout,
# then runs `cellranger count`.
#
# Usage:
#   bash scripts/rna/run_cellranger_count.sh <sample_id> <species> [start_step]
#
# Steps:
#   1  Prepare Cell Ranger FASTQs (Split_BC from R2 -> R1; copy cDNA R1 -> R2)
#   2  cellranger count
#
# Input (default: data/Spatial_RNA/<sample_id>/):
#   <sample_id>_R1.{fq,fastq}.gz   cDNA read  -> Cell Ranger R2
#   <sample_id>_R2.{fq,fastq}.gz   barcode read source -> Cell Ranger R1
#
# Output:
#   data/Spatial_RNA/<sample_id>/cellranger_fastq/              prepared FASTQs
#   data/Spatial_RNA/<sample_id>/cellranger/<sample_id>/        Cell Ranger run dir
#   data/Spatial_RNA/<sample_id>/raw_feature_bc_matrix/         raw matrix (3 files for downstream R)
#
# Environment variables:
#   AXIS_CELLRANGER_HUMAN_REF / AXIS_CELLRANGER_MOUSE_REF  Cell Ranger transcriptome
#   AXIS_RNA_DIR           RNA data root (default: <repo>/data/Spatial_RNA)
#   AXIS_NUM_CORES         localcores for cellranger (default: 30)
#   AXIS_CELLRANGER_MEM    localmem GB for cellranger (default: 64)
#   AXIS_CELLRANGER_BIN    cellranger executable (default: cellranger)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PREPROCESS_DIR="${REPO_ROOT}/scripts/preprocess"

sample_id="${1:?sample_id required}"
species="${2:?species required (human|mouse)}"
step="${3:-1}"

if [ "${species}" = "human" ]; then
    transcriptome="${AXIS_CELLRANGER_HUMAN_REF:-}"
elif [ "${species}" = "mouse" ]; then
    transcriptome="${AXIS_CELLRANGER_MOUSE_REF:-}"
else
    echo "Error: species must be 'human' or 'mouse'"
    exit 1
fi

if [ "${step}" -ge 2 ] && [ -z "${transcriptome}" ]; then
    echo "Error: set AXIS_CELLRANGER_HUMAN_REF or AXIS_CELLRANGER_MOUSE_REF"
    exit 1
fi

rna_root="${AXIS_RNA_DIR:-${REPO_ROOT}/data/Spatial_RNA}"
sample_dir="${rna_root}/${sample_id}"
fastq_dir="${sample_dir}/cellranger_fastq"
output_dir="${sample_dir}/cellranger"
cellranger_bin="${AXIS_CELLRANGER_BIN:-cellranger}"
num_cores="${AXIS_NUM_CORES:-30}"
local_mem="${AXIS_CELLRANGER_MEM:-64}"

mkdir -p "${fastq_dir}" "${output_dir}" logs

resolve_input() {
    local suffix="$1"
    for ext in fq fastq; do
        local candidate="${sample_dir}/${sample_id}_${suffix}.${ext}.gz"
        if [ -f "${candidate}" ]; then
            echo "${candidate}"
            return 0
        fi
    done
    return 1
}

r1_path="$(resolve_input R1)" || {
    echo "Cannot find R1 FASTQ for ${sample_id} in ${sample_dir}"
    echo "Expected: ${sample_id}_R1.fq.gz or ${sample_id}_R1.fastq.gz"
    exit 1
}
r2_path="$(resolve_input R2)" || {
    echo "Cannot find R2 FASTQ for ${sample_id} in ${sample_dir}"
    echo "Expected: ${sample_id}_R2.fq.gz or ${sample_id}_R2.fastq.gz"
    exit 1
}

cr_r1="${fastq_dir}/${sample_id}_S1_L001_R1_001.fastq.gz"
cr_r2="${fastq_dir}/${sample_id}_S1_L001_R2_001.fastq.gz"
cr_r1_plain="${fastq_dir}/${sample_id}_S1_L001_R1_001.fastq"

compress_fastq() {
    local plain="$1"
    local gz="${plain}.gz"
    if command -v pigz >/dev/null 2>&1; then
        pigz -p "${num_cores}" -f "${plain}"
    else
        gzip -f "${plain}"
    fi
    [ -f "${gz}" ] || { echo "Failed to compress ${plain}"; exit 1; }
}

publish_raw_matrix() {
    local raw_src="${output_dir}/${sample_id}/raw_feature_bc_matrix"
    local raw_dst="${sample_dir}/raw_feature_bc_matrix"
    local matrix_file

    [ -d "${raw_src}" ] || {
        echo "Missing Cell Ranger raw matrix directory: ${raw_src}"
        exit 1
    }

    mkdir -p "${raw_dst}"
    for matrix_file in barcodes.tsv.gz features.tsv.gz matrix.mtx.gz; do
        [ -f "${raw_src}/${matrix_file}" ] || {
            echo "Missing ${raw_src}/${matrix_file}"
            exit 1
        }
        cp -f "${raw_src}/${matrix_file}" "${raw_dst}/${matrix_file}"
    done

    echo "Published raw matrix for downstream R:"
    echo "  ${raw_dst}/"
    echo "    barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz"
}

if [ "${step}" -le 1 ]; then
    echo "[step 1] Prepare Cell Ranger FASTQs for ${sample_id}"
    python3 "${PREPROCESS_DIR}/split_rna_barcode.py" \
        -i "${r2_path}" \
        -o "${cr_r1_plain}" \
        2>"logs/${sample_id}.cellranger.split_bc.log"
    compress_fastq "${cr_r1_plain}"
    cp -f "${r1_path}" "${cr_r2}"
    echo "Prepared:"
    echo "  ${cr_r1}  (barcode read from R2)"
    echo "  ${cr_r2}  (cDNA read from R1)"
fi

if [ "${step}" -le 2 ]; then
    for req in "${cr_r1}" "${cr_r2}"; do
        [ -f "${req}" ] || { echo "Missing prepared FASTQ: ${req}. Run step 1 first."; exit 1; }
    done

    if ! command -v "${cellranger_bin}" >/dev/null 2>&1; then
        echo "Error: '${cellranger_bin}' not found in PATH"
        echo "Install Cell Ranger or set AXIS_CELLRANGER_BIN"
        exit 1
    fi

    echo "[step 2] cellranger count for ${sample_id}"
    "${cellranger_bin}" count \
        --id="${sample_id}" \
        --output-dir="${output_dir}" \
        --transcriptome="${transcriptome}" \
        --fastqs="${fastq_dir}" \
        --sample="${sample_id}" \
        --create-bam=true \
        --localcores="${num_cores}" \
        --localmem="${local_mem}" \
        2>&1 | tee "${REPO_ROOT}/logs/${sample_id}.cellranger_count.log"

    mri="${output_dir}/${sample_id}/${sample_id}.mri.tgz"
    if [ -f "${mri}" ]; then
        publish_raw_matrix
        echo "Done. Cell Ranger output: ${output_dir}/${sample_id}/"
    else
        echo "cellranger finished but ${mri} not found; check logs/${sample_id}.cellranger_count.log"
        exit 1
    fi
fi

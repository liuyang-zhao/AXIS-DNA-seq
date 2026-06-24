#!/bin/bash
# Generate genomic bin BED/TSV files from chromosome size files.
#
# Usage:
#   bash scripts/reference/generate_genomic_bins.sh
#
# Requires: bedtools

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

REFERENCE_DIR="${REPO_ROOT}/reference/chrom_sizes"
OUTPUT_DIR="${REPO_ROOT}/reference/genomic_bins"
PYTHON_SCRIPT="${SCRIPT_DIR}/bins_from_chrom_sizes.py"
version_mouse="GRCm39"
version_human="GRCh38"

mkdir -p "${OUTPUT_DIR}"

# Optional: tab-delimited bins with indices via bins_from_chrom_sizes.py
# python "${PYTHON_SCRIPT}" "${REFERENCE_DIR}/${version_mouse}.chrom.sizes" 1000000 > "${OUTPUT_DIR}/${version_mouse}_1Mb_bins.txt"

MOUSE_CHROM="${REFERENCE_DIR}/${version_mouse}.chrom.sizes"
HUMAN_CHROM="${REFERENCE_DIR}/${version_human}.chrom.sizes"

bedtools makewindows -g "${MOUSE_CHROM}" -w 10000000 > "${OUTPUT_DIR}/${version_mouse}_10Mb_bins.bed"
bedtools makewindows -g "${MOUSE_CHROM}" -w 5000000 > "${OUTPUT_DIR}/${version_mouse}_5Mb_bins.bed"
bedtools makewindows -g "${MOUSE_CHROM}" -w 2500000 > "${OUTPUT_DIR}/${version_mouse}_2.5Mb_bins.bed"
bedtools makewindows -g "${MOUSE_CHROM}" -w 1000000 > "${OUTPUT_DIR}/${version_mouse}_1Mb_bins.bed"
bedtools makewindows -g "${MOUSE_CHROM}" -w 500000 > "${OUTPUT_DIR}/${version_mouse}_500kb_bins.bed"
bedtools makewindows -g "${MOUSE_CHROM}" -w 250000 > "${OUTPUT_DIR}/${version_mouse}_250kb_bins.bed"
bedtools makewindows -g "${MOUSE_CHROM}" -w 200000 > "${OUTPUT_DIR}/${version_mouse}_200kb_bins.bed"
bedtools makewindows -g "${MOUSE_CHROM}" -w 100000 > "${OUTPUT_DIR}/${version_mouse}_100kb_bins.bed"
bedtools makewindows -g "${MOUSE_CHROM}" -w 50000 > "${OUTPUT_DIR}/${version_mouse}_50kb_bins.bed"

bedtools makewindows -g "${HUMAN_CHROM}" -w 10000000 > "${OUTPUT_DIR}/${version_human}_10Mb_bins.bed"
bedtools makewindows -g "${HUMAN_CHROM}" -w 1000000 > "${OUTPUT_DIR}/${version_human}_1Mb_bins.bed"
bedtools makewindows -g "${HUMAN_CHROM}" -w 250000 > "${OUTPUT_DIR}/${version_human}_250kb_bins.bed"
bedtools makewindows -g "${HUMAN_CHROM}" -w 100000 > "${OUTPUT_DIR}/${version_human}_100kb_bins.bed"
bedtools makewindows -g "${HUMAN_CHROM}" -w 50000 > "${OUTPUT_DIR}/${version_human}_50kb_bins.bed"
bedtools makewindows -g "${HUMAN_CHROM}" -w 25000 > "${OUTPUT_DIR}/${version_human}_25kb_bins.bed"
bedtools makewindows -g "${HUMAN_CHROM}" -w 10000 > "${OUTPUT_DIR}/${version_human}_10kb_bins.bed"

# Convert 1Mb BED to tabular bins with indices for downstream R/Python tools
for species in "${version_human}" "${version_mouse}"; do
    python "${PYTHON_SCRIPT}" "${REFERENCE_DIR}/${species}.chrom.sizes" 1000000 > "${OUTPUT_DIR}/${species}_1Mb_bins.txt"
done

echo "Genomic bins written to ${OUTPUT_DIR}"

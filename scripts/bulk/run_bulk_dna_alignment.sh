#!/bin/bash
# AXIS-DNA-seq bulk DNA alignment for normal-control generation.
#
# Usage:
#   bash scripts/bulk/run_bulk_dna_alignment.sh <run_id> <species> [start_step]
#
# Expects bulk FASTQ files in data/bulk/:
#   <run_id>_R1.{fq,fastq}.gz and <run_id>_R2.{fq,fastq}.gz
#
# Environment variables:
#   AXIS_GENOME_HUMAN_INDEX / AXIS_GENOME_MOUSE_INDEX  Bowtie2 index prefix
#   AXIS_BULK_DIR  Bulk data directory (default: <repo>/data/bulk)
#   AXIS_NUM_CORES Number of parallel jobs (default: 60)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

run="${1:?run_id required}"
species="${2:?species required (human|mouse)}"
step="${3:-1}"

if [ "${species}" = "human" ]; then
    genome="${AXIS_GENOME_HUMAN_INDEX:-}"
    bins="${REPO_ROOT}/reference/genomic_bins/GRCh38_1Mb_bins.bed"
elif [ "${species}" = "mouse" ]; then
    genome="${AXIS_GENOME_MOUSE_INDEX:-}"
    bins="${REPO_ROOT}/reference/genomic_bins/GRCm39_1Mb_bins.bed"
else
    echo "Error: Unknown species '${species}'"
    exit 1
fi

if [ -z "${genome}" ]; then
    echo "Error: Bowtie2 genome index not set."
    echo "Export AXIS_GENOME_HUMAN_INDEX or AXIS_GENOME_MOUSE_INDEX before running."
    exit 1
fi

home_dir="${AXIS_BULK_DIR:-${REPO_ROOT}/data/bulk}"
output="${home_dir}/${run}"
tmp="${home_dir}/tmp/"
num_cores="${AXIS_NUM_CORES:-60}"
split_size=4000000

mkdir -p "${output}" "${tmp}" logs

if [ -f "${home_dir}/${run}_R1.fastq.gz" ]; then
    r1_path="${home_dir}/${run}_R1.fastq.gz"
    r2_path="${home_dir}/${run}_R2.fastq.gz"
elif [ -f "${home_dir}/${run}_R1.fq.gz" ]; then
    r1_path="${home_dir}/${run}_R1.fq.gz"
    r2_path="${home_dir}/${run}_R2.fq.gz"
else
    echo "Cannot find R1 for run ${run} in ${home_dir}"
    exit 1
fi

r1="${run}_R1"
r2="${run}_R2"

echo "Run: ${run}"
echo "Genome path: ${genome}"
echo "Bins: ${bins}"
echo "Start step: ${step}"

## 1. Split FASTQ
mkdir -p "${tmp}/split/"

if [ "${step}" -le 1 ]; then
    echo "${run} - Splitting FASTQs (1)"
    fastqc -t 4 -o "${output}" "${r1_path}" "${r2_path}"
    fastp -i "${r1_path}" -o "${tmp}/split/${r1}.fastq.gz" -S "${split_size}" --thread 1 -d 4 -A -G -L -Q 2>logs/split_R1.log &
    fastp -i "${r2_path}" -o "${tmp}/split/${r2}.fastq.gz" -S "${split_size}" --thread 1 -d 4 -A -G -L -Q 2>logs/split_R2.log &
    wait
fi

ls "${tmp}/split/" | grep "${run}" | grep "R1" | grep -P -o "^[0-9]{4}" > "${output}/${run}.split_list.txt"

## 2. Trim adapters
mkdir -p "${tmp}/trim/"

if [ "${step}" -le 2 ]; then
    echo "${run} - Trimming adapters (2)"
    parallel --will-cite --jobs "${num_cores}" --colsep '\t' \
        trim_galore -j 12 --phred33 --length 10 -e 0.1 --stringency 4 --paired -o "${tmp}/trim/" \
        "${tmp}/split/{1}.${r1}.fastq.gz" "${tmp}/split/{1}.${r2}.fastq.gz" \
        1>logs/trim.log :::: "${output}/${run}.split_list.txt"
fi

## 3. Alignment
mkdir -p "${tmp}/aln/"

if [ "${step}" -le 3 ]; then
    echo "${run} - Aligning paired-end reads (3)"
    parallel --will-cite --jobs "${num_cores}" --colsep '\t' \
        bowtie2 -X2000 -p 1 --rg-id "${run}" \
        -x "${genome}" \
        -1 "${tmp}/trim/{1}.${r1}_val_1.fq.gz" \
        -2 "${tmp}/trim/{1}.${r2}_val_2.fq.gz" '|' \
        samtools view -bS - -o "${tmp}/aln/{1}.${run}.aln.bam" \
        2>logs/align.log :::: "${output}/${run}.split_list.txt"
fi

## 4. Sort
mkdir -p "${tmp}/sort/"

if [ "${step}" -le 4 ]; then
    echo "${run} - Sorting (4)"
    parallel --will-cite --jobs "${num_cores}" --colsep '\t' \
        samtools sort "${tmp}/aln/{1}.${run}.aln.bam" -o "${tmp}/sort/{1}.${run}.sort.bam" \
        2>logs/sort.log :::: "${output}/${run}.split_list.txt"
    parallel --will-cite --jobs "${num_cores}" --colsep '\t' \
        samtools index "${tmp}/sort/{1}.${run}.sort.bam" \
        2>logs/index_sort.log :::: "${output}/${run}.split_list.txt"
fi

## 5. Filter
mkdir -p "${tmp}/flt/"

if [ "${step}" -le 5 ]; then
    echo "Getting chromosome list..."
    chrs=$(samtools view -H "${tmp}/sort/0001.${run}.sort.bam" | grep chr | cut -f2 | sed 's/SN://g' | awk '{if(length($0)<6)print}')
    echo "Chromosomes to keep: ${chrs}"

    echo "${run} - Filtering by quality, chromosome, pairing (5)"
    parallel --will-cite --jobs "${num_cores}" --colsep '\t' \
        samtools view -@ "${num_cores}" -b -q 30 -f 0x2 "${tmp}/sort/{1}.${run}.sort.bam" -o "${tmp}/flt/{1}.${run}.flt.bam" ${chrs} \
        2>logs/filter.log :::: "${output}/${run}.split_list.txt"
    parallel --will-cite --jobs "${num_cores}" --colsep '\t' \
        samtools index "${tmp}/flt/{1}.${run}.flt.bam" \
        2>logs/index_flt.log :::: "${output}/${run}.split_list.txt"
fi

## 6. Merge
if [ "${step}" -le 6 ]; then
    echo "${run} - Merging all BAM files (6)"
    ls "${tmp}/flt/"*"${run}.flt.bam" > "${output}/${run}.merge_list.txt"
    samtools merge -@ "${num_cores}" -f -b "${output}/${run}.merge_list.txt" --threads "${num_cores}" "${output}/${run}.merge.bam" 2>logs/merge.log
    samtools index "${output}/${run}.merge.bam"
fi

## 7. 1Mb window coverage
if [ "${step}" -le 7 ]; then
    echo "${run} - Calculate 1Mb window coverage (7)"
    mosdepth --by 1000000 --threads "${num_cores}" "${output}/${run}" "${output}/${run}.merge.bam"
    zcat "${output}/${run}.regions.bed.gz" | awk '{print $1"\t"$2"\t"$3"\t"$4}' > "${output}/${run}.1Mb_coverage.txt"
    bedtools coverage -a "${bins}" -b "${output}/${run}.merge.bam" -counts > "${output}/${run}_1Mb_total_reads.txt"
fi

echo "==========================================="
echo "Bulk pipeline completed!"
echo "Results are in ${output}/"
echo "  - coverage: $(head -n 1 "${output}/${run}.1Mb_coverage.txt")"
echo "  - total reads: $(head -n 1 "${output}/${run}_1Mb_total_reads.txt")"
echo "==========================================="

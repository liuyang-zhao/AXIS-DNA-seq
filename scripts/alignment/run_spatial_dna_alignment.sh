#!/bin/bash
# AXIS-DNA-seq spatial DNA alignment pipeline.
#
# Usage:
#   bash scripts/alignment/run_spatial_dna_alignment.sh <project_dir> <run_id> <species> [start_step]
#
# Steps:
#   1  Split FASTQ          7  Sort BAM
#   2  Extract barcodes      8  Filter BAM (MAPQ>=30, properly paired)
#   3  Attach barcode        9  Write UM tag
#   4  Trim adapters        10  Merge BAM chunks
#   5  Concat clean FASTQ   11  UMI deduplication
#   6  Bowtie2 alignment    12  Sparse count matrix

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BIN_DIR="${REPO_ROOT}/scripts/preprocess"

home_dir="${1:?project_dir required}"
run="${2:?run_id required}"
species="${3:?species required (human|mouse)}"
step="${4:-1}"

if [ "${species}" = "human" ]; then
    genome="${AXIS_GENOME_HUMAN_INDEX:-}"
    bins="${REPO_ROOT}/reference/genomic_bins/GRCh38_1Mb_bins.txt"
elif [ "${species}" = "mouse" ]; then
    genome="${AXIS_GENOME_MOUSE_INDEX:-}"
    bins="${REPO_ROOT}/reference/genomic_bins/GRCm39_1Mb_bins.txt"
else
    echo "Error: species must be 'human' or 'mouse'"
    exit 1
fi

if [ -z "${genome}" ]; then
    echo "Error: set AXIS_GENOME_HUMAN_INDEX or AXIS_GENOME_MOUSE_INDEX"
    exit 1
fi

processed_dir="${home_dir}/processed"
output_dir="${home_dir}/results"
tmp="${home_dir}/tmp"
num_cores="${AXIS_NUM_CORES:-30}"
split_size=1000000
bin_size=1000000
bin_size_name="1Mb"
barcodes="${processed_dir}/${run}.spatial_barcodes_location.csv"

mkdir -p "${processed_dir}" "${output_dir}" "${tmp}" logs
cd "${home_dir}"

r1_path="${home_dir}/rawdata/${run}_R1.fq.gz"
r2_path="${home_dir}/rawdata/${run}_R2.fq.gz"
r1="${run}_R1"
r2="${run}_R2"
r3="${run}_R3"

for req in "${barcodes}" "${r1_path}" "${r2_path}"; do
    [ -f "${req}" ] || { echo "Missing required file: ${req}"; exit 1; }
done

if [ "${step}" -le 1 ]; then
    mkdir -p "${tmp}/split"
    seqkit stat "${r1_path}" "${r2_path}"
    fastp -i "${r1_path}" -o "${tmp}/split/${r1}.fastq.gz" -S "${split_size}" --thread 1 -d 4 -A -G -L -Q \
        2>"logs/${run}.1.split_R1.log" &
    fastp -i "${r2_path}" -o "${tmp}/split/${r2}.fastq.gz" -S "${split_size}" --thread 1 -d 4 -A -G -L -Q \
        2>"logs/${run}.1.split_R2.log" &
    wait
fi

ls "${tmp}/split/" | grep "${run}" | grep "R1" | grep -P -o "^[0-9]{4}" > "${processed_dir}/${run}_split_list.txt"

if [ "${step}" -le 2 ]; then
    mkdir -p "${tmp}/ext"
    parallel --will-cite --jobs "${num_cores}" \
        python "${BIN_DIR}/extract_spatial_barcodes.py" \
        --r1 "${tmp}/split/{1}.${r1}.fastq.gz" \
        --r2 "${tmp}/split/{1}.${r2}.fastq.gz" \
        --barcode "${barcodes}" \
        --out "${tmp}/ext/{1}.${run}" \
        --threads 3 \
        2>"logs/${run}.2.extract_barcode.log" \
        :::: "${processed_dir}/${run}_split_list.txt"
fi

if [ "${step}" -le 3 ]; then
    mkdir -p "${tmp}/umi"
    parallel --will-cite --jobs "${num_cores}" \
        python "${BIN_DIR}/attach_barcode_to_headers.py" \
        "${tmp}/ext/{1}.${r1}.ext.fastq.gz" "${tmp}/ext/{1}.${r2}.ext.fastq.gz" "${tmp}/ext/{1}.${r3}.ext.fastq.gz" \
        "${tmp}/umi/{1}.${r1}.umi.fastq.gz" "${tmp}/umi/{1}.${r2}.umi.fastq.gz" \
        2>"logs/${run}.3.attach_barcode.log" \
        :::: "${processed_dir}/${run}_split_list.txt"
fi

if [ "${step}" -le 4 ]; then
    mkdir -p "${tmp}/trim"
    parallel --will-cite --jobs "${num_cores}" \
        trim_galore -a CTGTCTCTTATACACATCT -a2 CTGTCTCTTATACACATCT -j 6 \
        --phred33 --length 10 -e 0.1 --stringency 4 --paired -o "${tmp}/trim/" \
        "${tmp}/umi/{1}.${r1}.umi.fastq.gz" "${tmp}/umi/{1}.${r2}.umi.fastq.gz" \
        2>"logs/${run}.4.trim.log" \
        :::: "${processed_dir}/${run}_split_list.txt"
    wait
fi

if [ "${step}" -le 5 ]; then
    cat "${tmp}/trim/"*."${r1}.umi_val_1.fq.gz" > "${processed_dir}/${r1}.clean.fastq.gz"
    cat "${tmp}/trim/"*."${r2}.umi_val_2.fq.gz" > "${processed_dir}/${r2}.clean.fastq.gz"
    cat "${tmp}/ext/"*."${r3}.ext.fastq.gz" > "${processed_dir}/${r3}.clean.fastq.gz"
fi

if [ "${step}" -le 6 ]; then
    mkdir -p "${tmp}/aln"
    parallel --will-cite --jobs "${num_cores}" \
        bowtie2 -X2000 -p1 --rg-id "${run}" -x "${genome}" \
        -1 "${tmp}/trim/{1}.${r1}.umi_val_1.fq.gz" \
        -2 "${tmp}/trim/{1}.${r2}.umi_val_2.fq.gz" '|' \
        samtools view -bS - -o "${tmp}/aln/{1}.${run}.aln.bam" \
        2>"logs/${run}.6.align.log" \
        :::: "${processed_dir}/${run}_split_list.txt"
fi

if [ "${step}" -le 7 ]; then
    mkdir -p "${tmp}/sort"
    parallel --will-cite --jobs "${num_cores}" \
        samtools sort "${tmp}/aln/{1}.${run}.aln.bam" -o "${tmp}/sort/{1}.${run}.sort.bam" \
        2>"logs/${run}.7.sort.log" \
        :::: "${processed_dir}/${run}_split_list.txt"
    parallel --will-cite --jobs "${num_cores}" \
        samtools index "${tmp}/sort/{1}.${run}.sort.bam" \
        2>"logs/${run}.7.index.log" \
        :::: "${processed_dir}/${run}_split_list.txt"
fi

if [ "${step}" -le 8 ]; then
    mkdir -p "${tmp}/flt"
    chrs=$(samtools view -H "${tmp}/sort/0001.${run}.sort.bam" | awk '/^@SQ/ {gsub("SN:", "", $2); if (length($2) < 6) print $2}')
    parallel --will-cite --jobs "${num_cores}" \
        samtools view -b -q 30 -f 0x2 "${tmp}/sort/{1}.${run}.sort.bam" -o "${tmp}/flt/{1}.${run}.flt.bam" ${chrs} \
        2>"logs/${run}.8.filter.log" \
        :::: "${processed_dir}/${run}_split_list.txt"
    parallel --will-cite --jobs "${num_cores}" \
        samtools index "${tmp}/flt/{1}.${run}.flt.bam" \
        2>"logs/${run}.8.index.log" \
        :::: "${processed_dir}/${run}_split_list.txt"
fi

if [ "${step}" -le 9 ]; then
    mkdir -p "${tmp}/tag"
    parallel --will-cite --jobs "${num_cores}" \
        python "${BIN_DIR}/write_barcode_umi_tag.py" \
        "${tmp}/flt/{1}.${run}.flt.bam" "${tmp}/tag/{1}.${run}.tag.bam" \
        2>"logs/${run}.9.umitag.log" \
        :::: "${processed_dir}/${run}_split_list.txt"
    parallel --will-cite --jobs "${num_cores}" \
        samtools index "${tmp}/tag/{1}.${run}.tag.bam" \
        2>"logs/${run}.9.index.log" \
        :::: "${processed_dir}/${run}_split_list.txt"
fi

if [ "${step}" -le 10 ]; then
    ls "${tmp}/tag/"*."${run}.tag.bam" > "${processed_dir}/${run}_merge_list.txt"
    samtools merge -f -b "${processed_dir}/${run}_merge_list.txt" --threads "${num_cores}" \
        "${processed_dir}/${run}.merge.bam"
    samtools index "${processed_dir}/${run}.merge.bam"
fi

if [ "${step}" -le 11 ]; then
    umi_tools group \
        -I "${processed_dir}/${run}.merge.bam" \
        --extract-umi-method=tag --umi-tag=UM:Z \
        --method=cluster --edit-distance-threshold=2 \
        --group-out="${processed_dir}/${run}.group.tsv" \
        --log="logs/${run}.11.group.log" \
        --paired --output-bam \
        -S "${processed_dir}/${run}.group.bam"
    grep -v "contig" "${processed_dir}/${run}.group.tsv" | cut -f2,3,7,8,9 | sort -u | sort -k5,5 -n | grep -Pv "\tN\t" \
        > "${output_dir}/${run}.reads.txt"
    cut -f3 "${output_dir}/${run}.reads.txt" | sort | uniq -c | sort -n -r > "${output_dir}/${run}.barcodes.txt"
fi

if [ "${step}" -le 12 ]; then
    python "${BIN_DIR}/generate_sparse_counts.py" \
        "${barcodes}" "${bins}" "${bin_size}" "${output_dir}/${run}.reads.txt" \
        > "${output_dir}/${run}.sparse_counts_${bin_size_name}.txt"
fi

echo "Pipeline completed: ${run}"

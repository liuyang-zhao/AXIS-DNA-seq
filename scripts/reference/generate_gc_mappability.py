#!/usr/bin/env python3
"""Calculate GC content and mappability annotations for genomic bins."""

import argparse
import os
import subprocess
import sys

import pysam


def calculate_gc_content(fasta_file, bins_file, output_file):
    print("Calculating GC content...")
    if not os.path.exists(fasta_file + ".fai"):
        print(f"Indexing {fasta_file} ...")
        pysam.faidx(fasta_file)

    fasta = pysam.FastaFile(fasta_file)

    with open(bins_file, "r", encoding="utf-8") as f_in, open(output_file, "w", encoding="utf-8") as f_out:
        f_in.readline()
        for line in f_in:
            fields = line.strip().split("\t")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            seq = fasta.fetch(chrom, start - 1, end)
            g_count = seq.upper().count("G")
            c_count = seq.upper().count("C")
            gc_total = g_count + c_count
            gc_fraction = 0.0 if len(seq) == 0 else gc_total / len(seq)
            f_out.write("\t".join(fields + [str(gc_total), str(gc_fraction)]) + "\n")

    fasta.close()
    print(f"GC content saved to {output_file}")


def calculate_signal_from_bigwig(feature_name, bigwig_file, bins_file, output_file):
    print(f"Calculating {feature_name} from {bigwig_file} ...")
    if not os.path.exists(bigwig_file):
        raise FileNotFoundError(f"bigWig file not found: {bigwig_file}")

    bins_data = {}
    original_header = []
    with open(bins_file, "r", encoding="utf-8") as f_in:
        original_header = f_in.readline().strip().split("\t")
        for line in f_in:
            fields = line.strip().split("\t")
            bin_key = f"{fields[0]}:{fields[1]}-{fields[2]}"
            bins_data[bin_key] = fields

    temp_bed_file = output_file + ".tmp.bed"
    with open(temp_bed_file, "w", encoding="utf-8") as f_out:
        for key, fields in bins_data.items():
            chrom, start, end = fields[0], int(fields[1]) - 1, int(fields[2])
            f_out.write(f"{chrom}\t{start}\t{end}\t{key}\n")

    command = f"bigWigAverageOverBed {bigwig_file} {temp_bed_file} stdout"
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        with open(output_file, "w", encoding="utf-8") as f_out:
            new_header = original_header + [f"{feature_name}_score"]
            f_out.write("\t".join(new_header) + "\n")
            for line in result.stdout.strip().split("\n"):
                res_fields = line.split("\t")
                bin_key, mean_signal = res_fields[0], res_fields[4]
                if bin_key in bins_data:
                    f_out.write("\t".join(bins_data[bin_key] + [mean_signal]) + "\n")
    finally:
        if os.path.exists(temp_bed_file):
            os.remove(temp_bed_file)

    print(f"{feature_name} saved to {output_file}")


def parse_args():
    parser = argparse.ArgumentParser(description="Generate GC and mappability reference files.")
    parser.add_argument("--species", choices=["human", "mouse"], required=True)
    parser.add_argument("--mode", choices=["gc", "map", "all"], default="all")
    parser.add_argument("--repo-root", default=None, help="AXIS-DNA-seq repository root")
    parser.add_argument("--fasta", default=None, help="Reference genome FASTA")
    parser.add_argument("--mappability-bw", default=None, help="Mappability bigWig file")
    return parser.parse_args()


def main():
    args = parse_args()
    repo_root = args.repo_root or os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

    if args.species == "human":
        genome_tag = "GRCh38"
        map_name = "hg38_1Mb_map.txt"
        default_fasta = os.environ.get("AXIS_HUMAN_FASTA", "")
        default_bw = os.environ.get("AXIS_HUMAN_MAP_BW", "")
    else:
        genome_tag = "GRCm39"
        map_name = "GRCm39_1Mb_map.txt"
        default_fasta = os.environ.get("AXIS_MOUSE_FASTA", "")
        default_bw = os.environ.get("AXIS_MOUSE_MAP_BW", "")

    bins_file = os.path.join(repo_root, "reference", "genomic_bins", f"{genome_tag}_1Mb_bins.txt")
    gc_file = os.path.join(repo_root, "reference", "gc_content", f"{genome_tag}_1Mb_gc.txt")
    map_file = os.path.join(repo_root, "reference", "mappability", map_name)

    fasta_file = args.fasta or default_fasta
    map_bw = args.mappability_bw or default_bw

    if args.mode in ("gc", "all"):
        if not fasta_file:
            sys.exit("FASTA path required. Use --fasta or export AXIS_HUMAN_FASTA / AXIS_MOUSE_FASTA.")
        calculate_gc_content(fasta_file, bins_file, gc_file)

    if args.mode in ("map", "all"):
        if not map_bw:
            sys.exit("Mappability bigWig required. Use --mappability-bw or export AXIS_HUMAN_MAP_BW / AXIS_MOUSE_MAP_BW.")
        calculate_signal_from_bigwig("mappability", map_bw, bins_file, map_file)


if __name__ == "__main__":
    main()

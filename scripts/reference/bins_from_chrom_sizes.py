#!/usr/bin/env python3
"""Convert a chrom.sizes file into fixed-width bins with indices."""

import math
import sys

import natsort


def main(chrom_sizes_path, bin_size):
    print("\t".join(["chr", "bin_start", "bin_end", "bin_len", "chr_ind", "bin_ind"]))

    chrom_sizes = {}
    with open(chrom_sizes_path, encoding="utf-8") as handle:
        for line in handle:
            chrom, size = line.rstrip().split()[:2]
            if "." not in chrom:
                chrom_sizes[chrom] = int(size)

    chrom_ind = 1
    bin_ind = 1
    ordered = [c for c in natsort.natsorted(chrom_sizes) if c != "chrM"]
    if "chrM" in chrom_sizes:
        ordered.append("chrM")

    for chrom in ordered:
        chrom_size = chrom_sizes[chrom]
        num_bins = int(math.ceil(chrom_size / bin_size))
        for i in range(num_bins):
            bin_start = max(1, i * bin_size)
            bin_end = min(chrom_size, (i + 1) * bin_size)
            print("\t".join([
                chrom,
                str(bin_start),
                str(bin_end),
                str(bin_end - bin_start),
                str(chrom_ind),
                str(bin_ind),
            ]))
            bin_ind += 1
        chrom_ind += 1


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <chrom.sizes> <bin_size>")
    main(sys.argv[1], int(sys.argv[2]))

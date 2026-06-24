#!/usr/bin/env python3
"""
Build a sparse spot x bin count matrix from deduplicated read assignments.

Output (stdout): row<TAB>col<TAB>count
  row = spot index (1-based, order in barcode file)
  col = genomic bin index (bin_ind from reference)
"""

import sys


def load_barcodes(path):
    barcodes = {}
    with open(path, encoding="utf-8") as handle:
        for idx, line in enumerate(handle, start=1):
            if not line.strip() or "barcode" in line.lower():
                continue
            barcode = line.strip().split(",")[0]
            barcodes[barcode] = idx
    return barcodes


def load_bins(path):
    bins = []
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().strip().split()
        try:
            bin_ind_col = header.index("bin_ind")
            chr_col = header.index("chr")
            start_col = header.index("bin_start")
            end_col = header.index("bin_end")
        except ValueError as exc:
            raise SystemExit(
                "Bins file must contain chr, bin_start, bin_end, bin_ind columns"
            ) from exc

        for line in handle:
            fields = line.rstrip().split()
            if len(fields) <= bin_ind_col:
                continue
            bins.append((
                fields[chr_col],
                int(fields[start_col]),
                int(fields[end_col]),
                int(fields[bin_ind_col]),
            ))
    return bins


def build_bin_lookup(bins):
    by_chr = {}
    for chrom, start, end, bin_ind in bins:
        by_chr.setdefault(chrom, []).append((start, end, bin_ind))
    for chrom in by_chr:
        by_chr[chrom].sort(key=lambda x: x[0])
    return by_chr


def locate_bin(by_chr, chrom, pos):
    intervals = by_chr.get(chrom)
    if not intervals:
        return None
    for start, end, bin_ind in intervals:
        if start <= pos <= end:
            return bin_ind
    return None


def main():
    if len(sys.argv) != 5:
        sys.exit(
            f"Usage: {sys.argv[0]} <barcodes.csv> <bins.txt> <bin_size> <reads.txt>"
        )

    barcodes = load_barcodes(sys.argv[1])
    bins = load_bins(sys.argv[2])
    by_chr = build_bin_lookup(bins)

    counts = {}
    with open(sys.argv[4], encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip().split()
            if len(fields) < 3:
                continue
            chrom, pos, barcode = fields[0], int(fields[1]), fields[2]
            spot_ind = barcodes.get(barcode)
            bin_ind = locate_bin(by_chr, chrom, pos)
            if spot_ind is None or bin_ind is None:
                continue
            key = (spot_ind, bin_ind)
            counts[key] = counts.get(key, 0) + 1

    for spot_ind, bin_ind in sorted(counts):
        print(f"{spot_ind}\t{bin_ind}\t{counts[(spot_ind, bin_ind)]}")


if __name__ == "__main__":
    main()

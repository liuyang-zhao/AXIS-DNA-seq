#!/usr/bin/env python3
"""
Extract spatial barcode + UMI from AXIS spatial RNA R2 reads for Cell Ranger.

Read structure on R2 (0-based slices):
  BC2  bp 1-8   (0:8)
  BC1  bp 39-46 (38:46)
  UMI  bp 77-86 (76:86)

Output is written as plain FASTQ (Cell Ranger R1: barcode read).
Adapted from /nas23/lab/yinshimei/s1/spatial_rna/CRC/src/python/Split_BC.py
"""

import argparse
import gzip
import sys

BC2 = slice(0, 8)
BC1 = slice(38, 46)
UMI = slice(76, 86)


def open_read(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(
        description="Extract BC2+BC1+UMI from AXIS spatial RNA R2 FASTQ."
    )
    parser.add_argument("-i", "--input", required=True, help="Input R2 FASTQ (.fq.gz)")
    parser.add_argument("-o", "--output", required=True, help="Output barcode FASTQ (uncompressed)")
    args = parser.parse_args()

    n_reads = 0
    with open_read(args.input) as in_handle, open(args.output, "w", encoding="utf-8") as out_handle:
        while True:
            header = in_handle.readline()
            if not header:
                break
            seq = in_handle.readline().rstrip("\n")
            plus = in_handle.readline()
            qual = in_handle.readline().rstrip("\n")
            if not qual:
                sys.exit(f"Malformed FASTQ near read {n_reads + 1} in {args.input}")

            barcode = seq[BC2] + seq[BC1] + seq[UMI]
            barcode_qual = qual[BC2] + qual[BC1] + qual[UMI]
            out_handle.write(f"{header}{barcode}\n+\n{barcode_qual}\n")
            n_reads += 1

    print(f"Wrote {n_reads} reads to {args.output}")


if __name__ == "__main__":
    main()

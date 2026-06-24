#!/usr/bin/env python3
"""Write spatial barcode from read name into the UM BAM tag for UMI-tools."""

import sys

import pysam


def main():
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <input.bam> <output.bam>")

    in_bam, out_bam = sys.argv[1], sys.argv[2]
    with pysam.AlignmentFile(in_bam, "rb") as src, pysam.AlignmentFile(out_bam, "wb", template=src) as dst:
        for read in src.fetch(until_eof=True):
            parts = read.query_name.split("+", 1)
            if len(parts) == 2:
                read.set_tag("UM", parts[1], value_type="Z")
            dst.write(read)


if __name__ == "__main__":
    main()

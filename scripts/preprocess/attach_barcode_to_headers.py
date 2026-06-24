#!/usr/bin/env python3
"""Append extracted spatial barcode sequences to R1/R2 read headers."""

import gzip
import sys
from itertools import zip_longest


def process_fastqs(r1_path, r2_path, r3_path, o1_path, o2_path):
    with gzip.open(r1_path, "rt") as r1, \
         gzip.open(r2_path, "rt") as r2, \
         gzip.open(r3_path, "rt") as r3, \
         gzip.open(o1_path, "wt") as o1, \
         gzip.open(o2_path, "wt") as o2:

        for bundle in zip_longest(*[r1] * 4, *[r2] * 4, *[r3] * 4):
            if None in bundle:
                continue

            (r1_header, r1_seq, _, r1_qual,
             r2_header, r2_seq, _, r2_qual,
             _, r3_seq, _, _) = bundle

            barcode = r3_seq.strip()
            r1_parts = r1_header.strip().split(maxsplit=1)
            r2_parts = r2_header.strip().split(maxsplit=1)

            new_r1_id = f"{r1_parts[0]}+{barcode}"
            new_r2_id = f"{r2_parts[0]}+{barcode}"
            new_r1_header = f"{new_r1_id} {r1_parts[1]}" if len(r1_parts) > 1 else new_r1_id
            new_r2_header = f"{new_r2_id} {r2_parts[1]}" if len(r2_parts) > 1 else new_r2_id
            new_r1_plus = f"+{new_r1_header[1:]}\n"
            new_r2_plus = f"+{new_r2_header[1:]}\n"

            o1.write(f"{new_r1_header}\n{r1_seq}{new_r1_plus}{r1_qual}")
            o2.write(f"{new_r2_header}\n{r2_seq}{new_r2_plus}{r2_qual}")


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(
            f"Usage: {sys.argv[0]} <in_r1.gz> <in_r2.gz> <in_r3.gz> <out_r1.gz> <out_r2.gz>"
        )
    process_fastqs(*sys.argv[1:6])

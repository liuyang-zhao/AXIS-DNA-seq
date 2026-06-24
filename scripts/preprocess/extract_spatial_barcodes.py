#!/usr/bin/env python3
"""
Extract spatial barcodes from R2 reads (1 mismatch allowed) and trim adapter sequence.

AXIS-DNA-seq read structure on R2:
  - bp 1-8 and 39-46: 16 bp spatial barcode
  - bp 96+: biological insert retained after trimming
"""

import argparse
import gzip
from functools import partial
from itertools import islice
from multiprocessing import Pool, cpu_count

PROGRESS_INTERVAL = 100_000


def load_valid_barcodes(barcode_file):
    """Read canonical barcode sequences from the first column of a CSV file."""
    barcodes = []
    with open(barcode_file, encoding="utf-8") as handle:
        for line_num, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            barcode = line.split(",")[0].strip()
            if line_num == 1 and "barcode" in barcode.lower():
                continue
            if barcode:
                barcodes.append(barcode)
    return barcodes


def build_barcode_index(valid_barcodes, max_mismatch=1):
    """Map all 1-mismatch neighbors to their canonical barcode."""
    bases = "ACGT"
    index = {}
    for bc in valid_barcodes:
        index.setdefault(bc, bc)
        if max_mismatch >= 1:
            for i, base in enumerate(bc):
                for alt in bases:
                    if alt == base:
                        continue
                    neighbor = bc[:i] + alt + bc[i + 1:]
                    index.setdefault(neighbor, bc)
    return index


def match_barcode(seq, barcode_index):
    """Extract 16 bp barcode from R2 and match against the index."""
    barcode = seq[:8] + seq[38:46]
    return barcode_index.get(barcode)


def process_read_block(read_block, barcode_index):
    r1_lines, r2_lines = read_block
    new_r1_lines, new_r2_lines, r3_lines = [], [], []

    for r1, r2 in zip(r1_lines, r2_lines):
        seq = r2[1].strip()
        barcode = match_barcode(seq, barcode_index)
        if not barcode:
            continue

        barcode_qual = r2[3].strip()[:8] + r2[3].strip()[38:46]
        header = r2[0].rstrip()
        r3_lines.extend([
            f"{header} {barcode}\n",
            f"{barcode}\n",
            "+\n",
            f"{barcode_qual}\n",
        ])

        new_seq = seq[95:]
        new_qual = r2[3].strip()[95:]
        new_r2_lines.extend([r2[0], f"{new_seq}\n", r2[2], f"{new_qual}\n"])
        new_r1_lines.extend(r1)

    return new_r1_lines, new_r2_lines, r3_lines


def read_fastq(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        while True:
            record = list(islice(handle, 4))
            if not record:
                break
            yield record


def chunked(iterator, size=10000):
    chunk = []
    for item in iterator:
        chunk.append(item)
        if len(chunk) == size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def main(r1_file, r2_file, barcode_file, out_prefix, threads, chunk_size=10000):
    valid_barcodes = load_valid_barcodes(barcode_file)
    if not valid_barcodes:
        raise SystemExit(f"No barcodes found in {barcode_file}")

    barcode_index = build_barcode_index(valid_barcodes, max_mismatch=1)
    out_r1 = gzip.open(out_prefix + "_R1.ext.fastq.gz", "wt")
    out_r2 = gzip.open(out_prefix + "_R2.ext.fastq.gz", "wt")
    out_r3 = gzip.open(out_prefix + "_R3.ext.fastq.gz", "wt")

    total_reads = 0
    next_report = PROGRESS_INTERVAL

    def chunk_generator():
        nonlocal total_reads, next_report
        for r1_chunk, r2_chunk in zip(
            chunked(read_fastq(r1_file), chunk_size),
            chunked(read_fastq(r2_file), chunk_size),
        ):
            total_reads += len(r1_chunk)
            while total_reads >= next_report:
                print(f"Processed {next_report} input reads")
                next_report += PROGRESS_INTERVAL
            yield (r1_chunk, r2_chunk)

    pool = Pool(threads)
    worker = partial(process_read_block, barcode_index=barcode_index)

    try:
        for new_r1_lines, new_r2_lines, r3_lines in pool.imap(worker, chunk_generator(), chunksize=10):
            out_r1.writelines(new_r1_lines)
            out_r2.writelines(new_r2_lines)
            out_r3.writelines(r3_lines)
    finally:
        pool.close()
        pool.join()
        out_r1.close()
        out_r2.close()
        out_r3.close()

    print(f"Finished. Total input reads: {total_reads}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract spatial barcodes from split FASTQ chunks.")
    parser.add_argument("--r1", required=True, help="Input R1 FASTQ (.gz supported)")
    parser.add_argument("--r2", required=True, help="Input R2 FASTQ (.gz supported)")
    parser.add_argument("--barcode", required=True, help="Sample spatial_barcodes_location.csv (first column = barcode)")
    parser.add_argument("--out", required=True, help="Output prefix, e.g. tmp/ext/0001.sample01")
    parser.add_argument("--threads", type=int, default=cpu_count(), help="Worker processes per chunk")
    parser.add_argument("--chunk-size", type=int, default=10000, help="Reads per in-memory block")
    args = parser.parse_args()
    main(args.r1, args.r2, args.barcode, args.out, args.threads, args.chunk_size)

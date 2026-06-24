# AXIS-DNA-seq

Spatial DNA sequencing (AXIS-DNA-seq) analysis pipeline for human and mouse samples. This repository provides alignment, bulk reference generation, CNV-based spatial clustering, multi-sample evolution analysis (MEDICC2), and RNA-DNA multi-omics alignment.

## Repository layout

```
AXIS-DNA-seq/
├── scripts/
│   ├── alignment/          # Spatial DNA read alignment
│   ├── bulk/               # Bulk DNA alignment (normal control)
│   ├── preprocess/         # Barcode extraction, UMI tagging, sparse counts
│   ├── reference/          # Genomic bin & annotation generation
│   ├── rna/                # Spatial RNA Cell Ranger pipeline
│   └── R/                  # CNV clustering, evolution, RNA-DNA registration
├── reference/              # Bundled 1 Mb bins, GC, mappability
├── data/
│   ├── samples/            # Sample manifest (tracked)
│   ├── bulk/               # Bulk FASTQ / counts (not tracked)
│   └── Spatial_DNA/        # Barcode CSVs (tracked); sparse counts per sample (not tracked)
└── config/                 # Example configuration
```

## Script index

| Script | Purpose |
|--------|---------|
| `scripts/alignment/run_spatial_dna_alignment.sh` | Main spatial DNA alignment → sparse counts |
| `scripts/bulk/run_bulk_dna_alignment.sh` | Bulk DNA alignment → 1 Mb bin counts |
| `scripts/reference/generate_genomic_bins.sh` | Generate genomic bin BED/TSV files |
| `scripts/reference/generate_gc_mappability.py` | Compute GC content and mappability |
| `scripts/preprocess/extract_spatial_barcodes.py` | Extract spatial barcodes from R2 |
| `scripts/preprocess/attach_barcode_to_headers.py` | Attach barcode to read headers |
| `scripts/preprocess/write_barcode_umi_tag.py` | Write barcode into UM BAM tag |
| `scripts/preprocess/generate_sparse_counts.py` | Build sparse spot × bin matrix |
| `scripts/preprocess/split_rna_barcode.py` | Extract BC2+BC1+UMI from spatial RNA R2 for Cell Ranger |
| `scripts/rna/run_cellranger_count.sh` | Spatial RNA FASTQ prep + `cellranger count` |
| `scripts/R/human_spatial_dna_analysis.R` | Human CNV clustering & multi-omics |
| `scripts/R/mouse_spatial_dna_analysis.R` | Mouse CNV clustering |
| `scripts/R/medicc2_evolution.R` | MEDICC2 export for merged samples |
| `scripts/R/spatial_rna_dna_registration.R` | RNA-DNA Procrustes registration QC |

## Prerequisites

### Command-line tools

- `fastp`, `seqkit`, `parallel`, `trim_galore`, `bowtie2`, `samtools`, `umi_tools`
- `bedtools`, `mosdepth` (bulk pipeline)
- `bigWigAverageOverBed` (UCSC tools, for mappability)
- `medicc2` (optional, for evolution tree)
- `cellranger` (optional, for spatial RNA gene expression; 10x Genomics)
- `pigz` (optional, faster FASTQ compression during RNA prep)

### Python

- Python 3.8+
- `pysam`, `natsort`

### R (>= 4.2 recommended)

Core: `data.table`, `Matrix`, `ggplot2`, `dplyr`, `FNN`, `patchwork`, `cluster`, `Rtsne`, `viridis`, `tidyr`, `fields`, `hexbin`, `ape`

Extended: `ComplexHeatmap`, `circlize`, `Morpho`, `monocle`, `ggtree`, `treeio`

RNA workflow (optional): `Seurat`, `SpotClean`, `infercnv`, `MASS`, `arrow`

## Reference genome index (not included)

Bowtie2 indexes are too large for GitHub. Build or download them locally, then export:

```bash
export AXIS_GENOME_HUMAN_INDEX=/path/to/GRCh38_index
export AXIS_GENOME_MOUSE_INDEX=/path/to/GRCm39_index
export AXIS_HUMAN_FASTA=/path/to/GRCh38.p14.genome.fa
export AXIS_MOUSE_FASTA=/path/to/GRCm39.genome.fa
export AXIS_HUMAN_MAP_BW=/path/to/hg38.k100.umap.bw
export AXIS_MOUSE_MAP_BW=/path/to/GRCm39_umap.bigWig
export AXIS_CELLRANGER_HUMAN_REF=/path/to/refdata-gex-GRCh38-2024-A
export AXIS_CELLRANGER_MOUSE_REF=/path/to/refdata-gex-GRCm39-2024-A
```

### Build Bowtie2 index example

```bash
bowtie2-build GRCh38.genome.fa GRCh38_index
bowtie2-build GRCm39.genome.fa GRCm39_index
```

## Data availability

Raw sequencing reads are deposited in the GSA database:

| Species | Accession | Browse |
|---------|-----------|--------|
| Human | HRA019095 | [GSA-Human HRA019095](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA019095) |
| Mouse | CRA044730 | [GSA CRA044730](https://ngdc.cncb.ac.cn/gsa/browse/CRA044730) |

Per-sample spatial barcode coordinate files (`spatial_barcodes_location.csv`) for DNA-seq samples are included under `data/Spatial_DNA/`. The full sample list is in [`data/samples/sample_manifest.tsv`](data/samples/sample_manifest.tsv).

| Sample name | Sample ID | Species | Assay | Sample title | Barcode file |
|-------------|-----------|---------|-------|--------------|--------------|
| H01 | 16_1 | human | DNA-seq | LSIL | `data/Spatial_DNA/16_1.spatial_barcodes_location.csv` |
| H02 | 18_1 | human | DNA-seq | HCC | `data/Spatial_DNA/18_1.spatial_barcodes_location.csv` |
| H03 | 11_2 | human | DNA-seq | CRC1 | `data/Spatial_DNA/11_2.spatial_barcodes_location.csv` |
| H04 | 11_1 | human | DNA-seq | CRC2 | `data/Spatial_DNA/11_1.spatial_barcodes_location.csv` |
| H05 | 17_1 | human | DNA-seq | CRC_biopsy_D | `data/Spatial_DNA/17_1.spatial_barcodes_location.csv` |
| H06 | 20_3 | human | DNA-seq | CRC_biopsy_C | `data/Spatial_DNA/20_3.spatial_barcodes_location.csv` |
| H07 | 20_1 | human | DNA-seq | CRC_biopsy_P | `data/Spatial_DNA/20_1.spatial_barcodes_location.csv` |
| H08 | 20_2 | human | DNA-seq | CRC_biopsy_L | `data/Spatial_DNA/20_2.spatial_barcodes_location.csv` |
| H09 | 19_3 | human | DNA-seq | CRC_surgical_C | `data/Spatial_DNA/19_3.spatial_barcodes_location.csv` |
| H10 | 19_1 | human | DNA-seq | CRC_surgical_P | `data/Spatial_DNA/19_1.spatial_barcodes_location.csv` |
| H11 | 19_2 | human | DNA-seq | CRC_surgical_L | `data/Spatial_DNA/19_2.spatial_barcodes_location.csv` |
| H12 | total_test10 | human | DNA-seq | CRC | `data/Spatial_DNA/total_test10.spatial_barcodes_location.csv` |
| H13 | 1_3 | human | RNA-seq | CRC | `data/Spatial_RNA/1_3/position_1_3.txt` |
| M01 | 11_4 | mouse | DNA-seq | normal cerebellum | `data/Spatial_DNA/11_4.spatial_barcodes_location.csv` |
| M02 | 13_1 | mouse | DNA-seq | lung tumor | `data/Spatial_DNA/13_1.spatial_barcodes_location.csv` |
| M03 | 12_4 | mouse | DNA-seq | subcutaneous tumor (rep1) | `data/Spatial_DNA/12_4.spatial_barcodes_location.csv` |
| M04 | 12_2 | mouse | DNA-seq | subcutaneous tumor (rep2) | `data/Spatial_DNA/12_2.spatial_barcodes_location.csv` |

Sample names starting with **H** are human (HRA019095); **M** are mouse (CRA044730). H13 (1_3) is RNA-seq only: place `<sample_id>_R1.fq.gz` and `<sample_id>_R2.fq.gz` under `data/Spatial_RNA/<sample_id>/`, then run `scripts/rna/run_cellranger_count.sh`. Also provide raw `position_<sample>.txt` (comma-separated `rowxcol` tissue outline); Section 0 of `spatial_rna_dna_registration.R` converts it to tab-delimited slide coordinates using the bundled `reference/barcode_index/spatial_barcodes_index.txt`.

To use a bundled barcode file for alignment, copy or symlink it into your project `processed/` directory:

```bash
cp data/Spatial_DNA/17_1.spatial_barcodes_location.csv /path/to/project/processed/17_1.spatial_barcodes_location.csv
```

## Quick start

### 1. Generate reference annotations (optional)

Bundled 1 Mb bins/GC/mappability are included. To regenerate:

```bash
bash scripts/reference/generate_genomic_bins.sh

python scripts/reference/generate_gc_mappability.py --species human --mode all \
  --fasta $AXIS_HUMAN_FASTA --mappability-bw $AXIS_HUMAN_MAP_BW

python scripts/reference/generate_gc_mappability.py --species mouse --mode all \
  --fasta $AXIS_MOUSE_FASTA --mappability-bw $AXIS_MOUSE_MAP_BW
```

### 2. Bulk DNA reference (normal control)

Place paired FASTQ in `data/bulk/`:

```
data/bulk/<sample>_R1.fq.gz
data/bulk/<sample>_R2.fq.gz
```

```bash
bash scripts/bulk/run_bulk_dna_alignment.sh <sample> human 1
# Output: data/bulk/<sample>/<sample>_1Mb_total_reads.txt
```

### 3. Spatial DNA alignment

Prepare per-sample inputs:

```
<project>/rawdata/<run>_R1.fq.gz
<project>/rawdata/<run>_R2.fq.gz
<project>/processed/<run>.spatial_barcodes_location.csv   # copy from data/Spatial_DNA/ or provide your own
```

```bash
export AXIS_GENOME_HUMAN_INDEX=/path/to/GRCh38_index

bash scripts/alignment/run_spatial_dna_alignment.sh \
  /path/to/project <run_id> human 1
```

Output: `results/<run>.sparse_counts_1Mb.txt`

For R analysis, place sparse counts under `data/Spatial_DNA/<expt_id>/`. Barcode coordinates are read from `data/Spatial_DNA/<expt_id>.spatial_barcodes_location.csv` (flat layout in repo root).

### 4. Spatial RNA Cell Ranger

Place paired RNA FASTQ under `data/Spatial_RNA/<sample_id>/`:

```
data/Spatial_RNA/<sample_id>/<sample_id>_R1.fq.gz   # cDNA
data/Spatial_RNA/<sample_id>/<sample_id>_R2.fq.gz   # spatial barcode read
data/Spatial_RNA/<sample_id>/position_<sample_id>.txt   # tissue outline (optional, for registration)
```

Example for H13 (`1_3`):

```bash
export AXIS_CELLRANGER_HUMAN_REF=/path/to/refdata-gex-GRCh38-2024-A

bash scripts/rna/run_cellranger_count.sh 1_3 human 1
# Step 1 only: prepare Cell Ranger FASTQs
# bash scripts/rna/run_cellranger_count.sh 1_3 human 2   # resume from cellranger count

# Outputs:
#   data/Spatial_RNA/1_3/cellranger_fastq/          prepared FASTQs
#   data/Spatial_RNA/1_3/cellranger/1_3/            Cell Ranger run directory
#   data/Spatial_RNA/1_3/raw_feature_bc_matrix/     raw matrix for downstream R
#     barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz
```

**Prep logic (from `src_rna/data_prepare.smk`):** R2 is split with `split_rna_barcode.py` (BC2+BC1+UMI) to Cell Ranger R1; original R1 is copied to Cell Ranger R2.

Downstream R scripts (`spatial_rna_dna_registration.R`) read **`raw_feature_bc_matrix`** (not filtered). After `cellranger count`, the shell script copies the three raw matrix files to `data/Spatial_RNA/<sample_id>/raw_feature_bc_matrix/`.

### 5. Human spatial CNV analysis

```r
source("scripts/R/human_spatial_dna_analysis.R")

base_dir <- "/path/to/AXIS-DNA-seq"
res <- run_unified_spatial_dna(
  expt_id = "sample_01",
  base_dir = base_dir,
  resolution = 80,
  bulk_control_path = file.path(base_dir, "data/bulk/normal_1Mb_total_reads.txt"),
  bins_path = file.path(base_dir, "reference/genomic_bins/GRCh38_1Mb_bins.txt"),
  gc_path = file.path(base_dir, "reference/gc_content/GRCh38_1Mb_gc.txt"),
  map_path = file.path(base_dir, "reference/mappability/hg38_1Mb_map.txt"),
  max_k_silhouette = 6,
  rna_loc_path = file.path(base_dir, "data/Spatial_RNA/sample_01/position_sample_01.txt")  # optional
)
```

**Outputs:** `<expt_id>_Spatial_Clusters.pdf`, `<expt_id>.Genome_Tracks_Overlap.pdf`, `<expt_id>_Silhouette_Scores.pdf`

**Analysis workflow:** KNN spatial smoothing → depth normalization → bulk log2 CNV → GC/mappability correction → variable-bin PCA → t-SNE → k-means (silhouette auto-k with `max_k_silhouette`, or manual `k_clusters`).

### 6. Mouse spatial CNV analysis

```r
source("scripts/R/mouse_spatial_dna_analysis.R")

run_slide_dna_mouse(
  expt_id = "13_1",
  base_dir = "/path/to/AXIS-DNA-seq",
  normal_bulk_path = file.path(base_dir, "data/bulk/mouse_blood_1Mb_total_reads.txt"),
  gc_path = file.path(base_dir, "reference/gc_content/GRCm39_1Mb_gc.txt"),
  map_path = file.path(base_dir, "reference/mappability/GRCm39_1Mb_map.txt")
)
```

### 7. MEDICC2 evolution analysis

```r
source("scripts/R/human_spatial_dna_analysis.R")

sample_ids <- c("Pre_D_1", "Pre_D_2", "Post_D_1", "Post_D_2")
res_list <- lapply(sample_ids, load_spatial_dna_res, base_dir = base_dir)
export_for_medicc2(res_list, file.path(base_dir, "results/MEDICC2_Input"), "Chemo_Evolution")
```

```bash
medicc2 ./results/MEDICC2_Input/Chemo_Evolution_medicc2_input.tsv \
  ./results/MEDICC2_Output/ \
  -j 8 --plot both --total-copy-numbers --normal-name Diploid_Root
```

For merged multi-sample clustering, see `scripts/R/medicc2_evolution.R`.

### 8. RNA-DNA multi-omics alignment

1. **Integrated** — pass `rna_loc_path` to `run_unified_spatial_dna()`.
2. **Standalone** — `scripts/R/spatial_rna_dna_registration.R` section 4 (Procrustes registration QC).

## Input file formats

| File | Description |
|------|-------------|
| `*.spatial_barcodes_location.csv` | Per-sample spot table; column 1 = barcode sequence, plus coordinates |
| `*.sparse_counts_1Mb.txt` | Sparse count matrix: `row<TAB>col<TAB>count` |
| `*_1Mb_total_reads.txt` | Bulk 1 Mb bin counts (chr, start, end, count) |
| `<sample_id>_R1.fq.gz` | Spatial RNA cDNA read (Cell Ranger R2 input) |
| `<sample_id>_R2.fq.gz` | Spatial RNA barcode read (source for Cell Ranger R1) |
| `position_<sample_id>.txt` | Comma-separated `rowxcol` tissue outline for RNA spots |
| `raw_feature_bc_matrix/` | Cell Ranger raw count matrix: `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` |

## Environment variables

| Variable | Purpose |
|----------|---------|
| `AXIS_GENOME_HUMAN_INDEX` | Bowtie2 index prefix (human) |
| `AXIS_GENOME_MOUSE_INDEX` | Bowtie2 index prefix (mouse) |
| `AXIS_REPO_ROOT` | Repository root for R scripts |
| `AXIS_BULK_DIR` | Bulk FASTQ directory (default: `data/bulk`) |
| `AXIS_RNA_DIR` | Spatial RNA FASTQ root (default: `data/Spatial_RNA`) |
| `AXIS_CELLRANGER_HUMAN_REF` | Cell Ranger transcriptome (human) |
| `AXIS_CELLRANGER_MOUSE_REF` | Cell Ranger transcriptome (mouse) |
| `AXIS_CELLRANGER_BIN` | `cellranger` executable path (default: `cellranger`) |
| `AXIS_CELLRANGER_MEM` | Cell Ranger `--localmem` in GB (default: 64) |
| `AXIS_NUM_CORES` | Parallel job count (default: 30/60) |

## Notes

- Barcode coordinates for published samples are in `data/Spatial_DNA/` (see [Data availability](#data-availability)). Before alignment, place `<run>.spatial_barcodes_location.csv` under your project `processed/` directory. Barcode extraction uses **column 1** of this file (one barcode sequence per spot).
- Raw FASTQ, sparse count matrices, and genome indexes are excluded via `.gitignore`.
- Replication timing (`rep_timing`) reference generation is excluded.

## License

MIT — see [LICENSE](LICENSE).

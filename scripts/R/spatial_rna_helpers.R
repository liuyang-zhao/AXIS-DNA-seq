# AXIS-DNA-seq: Spatial RNA preprocessing helpers.

DEFAULT_SPATIAL_RESOLUTION <- 80L
GRID_CELL_MM <- 1.5
POINT_OVERLAP_FACTOR <- 1.6
PLOT_LEGEND_MARGIN_INCH <- 1.5

spatial_plot_dim_inch <- function(resolution) {
  res_num <- as.numeric(resolution)
  (res_num * GRID_CELL_MM) / 25.4 + PLOT_LEGEND_MARGIN_INCH
}

spatial_point_size <- function(resolution, panel_mm = 89) {
  res_num <- as.numeric(resolution)
  (panel_mm / res_num) * POINT_OVERLAP_FACTOR
}

spatial_grid_point_size <- function(resolution) {
  plot_dim_inch <- spatial_plot_dim_inch(resolution)
  panel_mm <- plot_dim_inch * 25.4 - PLOT_LEGEND_MARGIN_INCH * 25.4
  spatial_point_size(resolution, panel_mm = panel_mm)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

# Load per-sample barcode grid from spatial_barcodes_location.csv.
# Uses column 1 as barcode; prefers explicit row/col, otherwise xcoord/ycoord.
load_spatial_barcode_index <- function(path, resolution = DEFAULT_SPATIAL_RESOLUTION) {
  if (!file.exists(path)) {
    stop("Barcode location file not found: ", path)
  }

  df <- read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(df) == 0 || nrow(df) == 0) {
    stop("Empty barcode location file: ", path)
  }

  first_col_name <- as.character(names(df)[1])
  if (grepl("^[ACGT]+$", first_col_name, ignore.case = TRUE)) {
    df <- read.csv(path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
    default_names <- c("barcode", "xcoord", "ycoord", "row", "col")
    names(df) <- default_names[seq_len(ncol(df))]
  }

  names(df) <- tolower(names(df))
  barcodes <- df[[1]]

  if ("row" %in% names(df) && "col" %in% names(df)) {
    return(data.frame(
      barcode = barcodes,
      row = df$row,
      col = df$col,
      stringsAsFactors = FALSE
    ))
  }

  if ("xcoord" %in% names(df) && "ycoord" %in% names(df)) {
    return(data.frame(
      barcode = barcodes,
      row = df$xcoord,
      col = df$ycoord,
      stringsAsFactors = FALSE
    ))
  }

  grid_res <- as.integer(resolution)
  if (is.na(grid_res) || grid_res <= 0) {
    grid_res <- DEFAULT_SPATIAL_RESOLUTION
  }
  idx <- seq_along(barcodes)
  data.frame(
    barcode = barcodes,
    row = ((idx - 1) %% grid_res) + 1,
    col = ((idx - 1) %/% grid_res) + 1,
    stringsAsFactors = FALSE
  )
}

resolve_spatial_barcode_location <- function(sample_id, repo_root = ".", species = NULL) {
  candidates <- c(
    file.path(repo_root, "data/Spatial_DNA", paste0(sample_id, ".spatial_barcodes_location.csv")),
    file.path(repo_root, "data/Spatial_RNA", sample_id, paste0(sample_id, ".spatial_barcodes_location.csv")),
    file.path(repo_root, "data/Spatial_DNA", sample_id, paste0(sample_id, ".spatial_barcodes_location.csv")),
    file.path(repo_root, "processed", paste0(sample_id, ".spatial_barcodes_location.csv"))
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    stop(
      "Could not find spatial_barcodes_location.csv for sample ", sample_id,
      ". Expected one of:\n  ", paste(candidates, collapse = "\n  ")
    )
  }
  found[[1]]
}

# Parse raw tissue-outline file: one comma-separated line of rowxcol tokens (e.g. 9x1,10x1,...).
parse_raw_tissue_position <- function(path) {
  if (!file.exists(path)) {
    stop("Tissue position file not found: ", path)
  }
  raw <- read.table(path, sep = ",", header = FALSE, stringsAsFactors = FALSE)
  tokens <- as.character(unlist(raw[1, ], use.names = FALSE))
  tokens <- tokens[!is.na(tokens) & nzchar(tokens)]
  if (length(tokens) == 0) {
    stop("No tissue coordinates found in: ", path)
  }
  parts <- strsplit(tokens, "x", fixed = TRUE)
  invalid <- vapply(parts, length, integer(1)) != 2L
  if (any(invalid)) {
    stop("Invalid coordinate token(s) in ", path, ": ", paste(tokens[which(invalid)[1]], collapse = ", "))
  }
  data.frame(
    xcoord = as.integer(vapply(parts, `[`, 1, FUN.VALUE = character(1))),
    ycoord = as.integer(vapply(parts, `[`, 2, FUN.VALUE = character(1))),
    stringsAsFactors = FALSE
  )
}

resolve_spatial_rna_barcode_index <- function(sample_id, repo_root = ".", resolution = DEFAULT_SPATIAL_RESOLUTION) {
  sample_path <- tryCatch(
    resolve_spatial_barcode_location(sample_id, repo_root),
    error = function(e) NULL
  )
  if (!is.null(sample_path)) {
    return(load_spatial_barcode_index(sample_path, resolution = resolution))
  }
  bundled_index <- file.path(repo_root, "reference/barcode_index/spatial_barcodes_index.txt")
  if (file.exists(bundled_index)) {
    idx <- read.table(bundled_index, sep = "\t", header = FALSE, col.names = c("barcode", "row", "col"))
    return(idx[, c("barcode", "row", "col"), drop = FALSE])
  }
  stop(
    "Could not resolve barcode index for sample ", sample_id,
    ". Provide a per-sample spatial_barcodes_location.csv or bundled reference/barcode_index/spatial_barcodes_index.txt."
  )
}

# Map tissue pixel coordinates to spatial barcode grid coordinates.
# barcode_index: data.frame with columns barcode, row, col
# location: pixel coordinates from tissue image (x, y)
# tissue_flag: 1 = in tissue
# Resolve Cell Ranger raw_feature_bc_matrix directory (barcodes/features/matrix .gz trio).
resolve_cellranger_raw_matrix_dir <- function(sample_id, repo_root = ".") {
  required <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
  candidates <- c(
    file.path(repo_root, "data/Spatial_RNA", sample_id, "raw_feature_bc_matrix"),
    file.path(repo_root, "data/Spatial_RNA", sample_id, "cellranger", sample_id, "raw_feature_bc_matrix")
  )

  for (path in candidates) {
    if (!dir.exists(path)) {
      next
    }
    missing <- required[!file.exists(file.path(path, required))]
    if (length(missing) == 0L) {
      return(normalizePath(path, winslash = "/"))
    }
  }

  stop(
    "Could not find raw_feature_bc_matrix for sample ", sample_id,
    ". Expected barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz under:\n  ",
    paste(candidates, collapse = "\n  "),
    "\nRun scripts/rna/run_cellranger_count.sh first."
  )
}

fetch_side_infor <- function(barcode_index, location, tissue_flag = 1) {
  colnames(barcode_index) <- c("barcode", "row", "col")
  location <- as.data.frame(location)
  if (!all(c("xcoord", "ycoord") %in% names(location))) {
    names(location)[1:2] <- c("xcoord", "ycoord")
  }

  location$tissue <- tissue_flag
  location$barcode <- barcode_index$barcode[match(
    paste(location$xcoord, location$ycoord),
    paste(barcode_index$row, barcode_index$col)
  )]

  matched <- !is.na(location$barcode)
  if (!any(matched)) {
    stop("No tissue coordinates matched the barcode grid. Check resolution and barcode index.")
  }

  location$row <- location$xcoord
  location$col <- location$ycoord
  location <- location[matched, c("barcode", "row", "col", "xcoord", "ycoord", "tissue")]
  rownames(location) <- NULL
  location
}

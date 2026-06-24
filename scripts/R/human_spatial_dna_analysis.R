# AXIS-DNA-seq: Human spatial DNA CNV analysis
# Pipeline: KNN smoothing -> bulk-normalized log2 CNV -> PCA/t-SNE -> k-means clustering
# Optional: RNA-DNA Procrustes alignment, MEDICC2 export

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(FNN)
  library(patchwork)
  library(viridis)
  library(cluster)
  library(Rtsne)
  library(tidyr)
  library(fields)
})

optional_r_packages <- c("Morpho", "monocle")
for (pkg in optional_r_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

require_optional_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for this step. Please install it first.", pkg), call. = FALSE)
  }
}

select <- dplyr::select

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

# ===================== 1. Helper Functions ===============================
TRACK_CLUSTER_PALETTE <- c("1" = "#4A90E2", "2" = "#7ED321", "3" = "#FF9500", "4" = "#E41A1C",
                           "5" = "#377EB8", "6" = "#4DAF4A", "7" = "#F3E5E2", "Unmapped" = "#E0E0E0")

smooth_matrix_knn <- function(data_mat, coords, k=50) {
  cat(sprintf("   -> Running KNN smoothing (k=%d)...\n", k))
  knn_res <- FNN::get.knn(coords, k = k)
  indices <- knn_res$nn.index
  
  smoothed <- matrix(0, nrow = nrow(data_mat), ncol = ncol(data_mat))
  pb <- txtProgressBar(min = 0, max = nrow(data_mat), style = 3)
  for(i in 1:nrow(data_mat)) {
    smoothed[i, ] <- colMeans(data_mat[indices[i, ], , drop=FALSE])
    if(i %% 500 == 0) setTxtProgressBar(pb, i)
  }
  close(pb)
  return(smoothed)
}

plot_spatial_cns <- function(beads_df, resolution, color_col, title_str, palette_type = "viridis", discrete = FALSE) {
  res_num <- as.numeric(resolution)
  plot_dim_inch <- spatial_plot_dim_inch(res_num)
  panel_mm <- plot_dim_inch * 25.4 - PLOT_LEGEND_MARGIN_INCH * 25.4
  point_size <- spatial_point_size(res_num, panel_mm = panel_mm)

  p <- ggplot(beads_df, aes(x = xcoord, y = ycoord, color = .data[[color_col]])) +
    geom_point(size = point_size, shape = 16) +
    coord_fixed(ratio = 1, xlim = c(0.5, res_num + 0.5), ylim = c(res_num + 0.5, 0.5), expand = FALSE) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 5, 5, 5),
      axis.line = element_line(color = "black", linewidth = 0.1), # axis line
      axis.ticks = element_line(color = "black", linewidth = 0.5), # axis ticks
      axis.text = element_text(color = "black", size = 10),        # axis labels
    ) +
    labs(title = title_str, color = "")
  
  if (discrete) {
    p <- p + scale_color_manual(values = c("#4A90E2", "#7ED321", "#FF9500", "#E41A1C", "#377EB8", "#4DAF4A", "#F3E5E2"))
  } else {
    if (palette_type == "viridis") p <- p + scale_color_viridis_c(option = "magma")
    else if (palette_type == "redblue") p <- p + scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0)
  }
  return(p)
}

plot_alignment_cns <- function(df, color_col, title_str, discrete = TRUE) {
  p <- ggplot(df, aes(x = x_norm, y = y_norm, color = .data[[color_col]])) +
    geom_point(size = 3.0, shape = 16, alpha = 0.8) +
    scale_x_continuous(limits = c(-0.05, 1.05), expand = c(0, 0)) + 
    scale_y_continuous(limits = c(-0.05, 1.05), expand = c(0, 0)) + 
    coord_fixed() + theme_void() + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "bottom",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5) 
    ) + labs(title = title_str, color = "Cluster")
  if (discrete) p <- p + scale_color_manual(values = TRACK_CLUSTER_PALETTE, na.value = "#E0E0E0")
  return(p)
}

# LOWESS normalization: fit on aggregated selected beads, then divide each bin by fit_all.
lowess_norm_suit_style <- function(counts_mat, track_vec, sel_fit, sel_beads, f = 0.1) {
  counts_mat <- as.matrix(counts_mat)
  track_vec <- as.numeric(track_vec)
  sel_fit <- as.logical(sel_fit)
  sel_beads <- as.logical(sel_beads)
  
  if (ncol(counts_mat) != length(track_vec)) {
    stop("lowess_norm_suit_style: counts_mat columns must match track length.")
  }
  if (nrow(counts_mat) != length(sel_beads)) {
    stop("lowess_norm_suit_style: sel_beads length must match counts_mat rows.")
  }
  
  agg_signal <- colSums(counts_mat[sel_beads, , drop = FALSE], na.rm = TRUE)
  x_fit <- track_vec[sel_fit]
  y_fit <- agg_signal[sel_fit]
  keep <- is.finite(x_fit) & is.finite(y_fit)
  
  if (sum(keep) < 20) {
    fit_all <- rep(1, length(track_vec))
    return(list(counts_norm = counts_mat, fit_all = fit_all))
  }
  
  lo <- stats::lowess(x_fit[keep], y_fit[keep], f = f, iter = 3)
  fit_all <- stats::approx(
    x = lo$x, y = lo$y, xout = track_vec, method = "linear", ties = mean, rule = 2
  )$y
  
  fallback <- median(fit_all[is.finite(fit_all) & fit_all > 0], na.rm = TRUE)
  if (!is.finite(fallback) || fallback <= 0) fallback <- 1
  fit_all[!is.finite(fit_all) | fit_all <= 0] <- fallback
  
  counts_norm <- sweep(counts_mat, 2, fit_all, "/")
  list(counts_norm = counts_norm, fit_all = fit_all)
}




# ===================== 2. Main Analysis Pipeline ==============================

run_unified_spatial_dna <- function(expt_id, base_dir, resolution = DEFAULT_SPATIAL_RESOLUTION, bulk_control_path, bins_path, gc_path, map_path,
                                    k_clusters = NULL, max_k_silhouette = 6L, rna_loc_path = NULL) {
  
  # Setup parameter
  gc_thresh <- 0.35
  map_thresh <- 0.7
  
  # --- Setup Directories ---
  spatial_dna_root <- file.path(base_dir, "data/Spatial_DNA")
  data_dir <- file.path(spatial_dna_root, expt_id)
  result_dir <- file.path(base_dir, "results/Spatial_DNA", expt_id)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat(sprintf("\n=== Starting Pipeline for %s ===\n", expt_id))
  
  # ---- 1. Load Data ----
  cat("1. Loading Matrix and Coordinates...\n")
  sparse_path <- file.path(data_dir, paste0(expt_id, ".sparse_counts_1Mb.txt"))
  beads_path <- file.path(spatial_dna_root, paste0(expt_id, ".spatial_barcodes_location.csv"))
  
  sparse_dt <- fread(sparse_path, header = FALSE, col.names = c("row", "col", "val"))
  beads <- fread(beads_path)
  bins <- fread(bins_path, header = TRUE)
  
  # Ensure bins are sorted by chromosome and position
  chr_levels <- c(as.character(1:22), "X", "Y")
  bins$chr_clean <- gsub("chr", "", bins$chr)
  bins$chr_factor <- factor(bins$chr_clean, levels = chr_levels)
  bins <- bins[order(bins$chr_factor, bins$bin_start), ]
  
  # load and align GC / mappability annotations
  cat("   -> Loading GC and Mappability annotations...\n")
  # GC file: no header; column 8 is gc_ratio
  gc_dt <- fread(gc_path, header = FALSE)
  gc_dt <- gc_dt[, c(1, 2, 3, 8), with = FALSE]
  colnames(gc_dt) <- c("chr", "bin_start", "bin_end", "gc_ratio")
  
  # mappability file has header
  map_dt <- fread(map_path, header = TRUE)
  map_dt <- map_dt[, c("chr", "bin_start", "bin_end", "mappability_score"), with = FALSE]
  
  # join to main bins table
  bins <- left_join(bins, gc_dt, by = c("chr", "bin_start", "bin_end"))
  bins <- left_join(bins, map_dt, by = c("chr", "bin_start", "bin_end"))
  
  # impute missing GC with mean; missing map -> 0 (untrusted)
  bins$gc_ratio[is.na(bins$gc_ratio)] <- mean(bins$gc_ratio, na.rm=TRUE)
  bins$mappability_score[is.na(bins$mappability_score)] <- 0
  
  # bin selection: gc > 0.35 & map > 0.7
  sel_bins <- bins$gc_ratio > gc_thresh & bins$mappability_score > map_thresh
  sel_bins[is.na(sel_bins)] <- FALSE
  if (sum(sel_bins) == 0) stop("No bins pass gc/map thresholds. Check gc_path/map_path or thresholds.")
  cat(sprintf("   -> Selected %d / %d bins by thresholds (gc > %.2f, map > %.2f)\n", sum(sel_bins), length(sel_bins), gc_thresh, map_thresh))

  # Calculate Cumulative Position for Plotting (Manhattan Plot style)
  bins_grouped <- aggregate(bin_end ~ chr_factor, data = bins, FUN = max)
  colnames(bins_grouped)[colnames(bins_grouped) == "bin_end"] <- "max_pos"
  bins_grouped <- bins_grouped[order(bins_grouped$chr_factor), , drop = FALSE]
  bins_grouped$chr_shift <- c(0, head(cumsum(as.numeric(bins_grouped$max_pos)), -1))

  bins <- merge(bins, bins_grouped[,c("chr_factor","chr_shift")], by="chr_factor", all.x=TRUE)
  bins$cum_start <- bins$bin_start + bins$chr_shift
  bins <- bins[order(bins$chr_factor, bins$bin_start), ] # Ensure order persists
  
  counts_mat <- sparseMatrix(i = sparse_dt$row, j = sparse_dt$col, x = sparse_dt$val)
  
  
  # === Load Normal Bulk ===
  cat("   -> Loading and aligning normal bulk data...\n")
  bulk_raw_df <- fread(bulk_control_path, header = FALSE)
  colnames(bulk_raw_df)[1:4] <- c("bulk_chr", "bulk_start", "bulk_end", "bulk_count")
  # align bulk to bins by chromosome and bin_end
  aligned_bulk <- dplyr::left_join(bins, bulk_raw_df, by = c("chr" = "bulk_chr", "bin_end" = "bulk_end"))
  norm_bulk_vec <- aligned_bulk$bulk_count
  norm_bulk_vec[is.na(norm_bulk_vec)] <- 0
  
  cat(sprintf("   -> Aligned bulk data to %d genomic bins.\n", length(norm_bulk_vec)))
  
  # normalize bulk baseline
  bulk_mean <- mean(norm_bulk_vec[norm_bulk_vec > 0], na.rm = TRUE)
  if(is.na(bulk_mean) || bulk_mean == 0) bulk_mean <- 1 
  norm_bulk_baseline <- norm_bulk_vec / bulk_mean
  norm_bulk_baseline[norm_bulk_baseline <= 0] <- 1e-6
  
  
  # ---- 2. Filter & Preprocessing ----
  cat("2. Preprocessing...\n")
  #beads$ycoord <- -beads$ycoord
  
  # Keep original matrix for V3-style clustering branch
  counts_mat_all <- counts_mat
  bins_all <- bins
  
  # Coverage filtering is shared by both branches
  keep_cells <- rowSums(counts_mat_all) > 100
  beads_filt <- beads[keep_cells, ]
  counts_filt_all <- counts_mat_all[keep_cells, ]
  coords <- as.matrix(beads_filt[, .(xcoord, ycoord)])
  
  # Shared XY smoothing
  counts_smo_xy_all <- smooth_matrix_knn(as.matrix(counts_filt_all), coords, k = 50)

  
  # ---- 3. Clustering ----
  # ---- Clustering branch A: spatially smoothed data ----
  cat("A. Clustering ...\n")
  cell_means_smo <- rowMeans(counts_smo_xy_all)
  counts_norm_smo <- counts_smo_xy_all / cell_means_smo
  
  bin_means <- colMeans(counts_norm_smo)
  bin_sds <- apply(counts_norm_smo, 2, sd)
  bin_cv <- bin_sds / bin_means
  pc_bins <- which(bin_means > 0 & bin_means < 5 & bin_cv < 1)
  cat(sprintf("   -> Selected %d genomic bins for PCA\n", length(pc_bins)))
  
  pca_res <- prcomp(counts_norm_smo[, pc_bins], center = TRUE, scale. = FALSE, rank. = 50)
  
  cell_means_raw <- rowMeans(counts_filt_all)
  counts_norm_raw <- counts_filt_all / cell_means_raw
  raw_data_subset <- as.matrix(counts_norm_raw[, pc_bins, drop = FALSE])
  pc_scores <- scale(raw_data_subset, center = pca_res$center, scale = FALSE) %*% pca_res$rotation
  
  pc_dist_input <- pc_scores[, 1:10, drop = FALSE]
  pc_smo_pc <- smooth_matrix_knn(pc_scores, pc_dist_input, k = 50)
  pc_smo_both <- smooth_matrix_knn(pc_smo_pc, coords, k = 10)
  
  set.seed(123)
  tsne_input <- pc_smo_both[, 1:10, drop = FALSE]
  
  # adaptive t-SNE perplexity
  n_samples <- nrow(tsne_input)
  max_perplexity <- floor((n_samples - 1) / 3)
  
  # default 80; reduce for small n; minimum 2
  current_perplexity <- min(80, max_perplexity)
  if (current_perplexity < 5) {
    current_perplexity <- max(2, current_perplexity)
    warning(sprintf("   -> Very few spots (%d); perplexity set to %d", n_samples, current_perplexity))
  } else if (current_perplexity < 80) {
    cat(sprintf("   -> Sample size %d; reduced perplexity from 80 to %d\n", n_samples, current_perplexity))
  }
  
  set.seed(123)
  tsne_res <- Rtsne(tsne_input, dims = 2, perplexity = current_perplexity, verbose = TRUE, check_duplicates = FALSE)
  
  tsne_coords <- tsne_res$Y
  beads_filt$tSNE_1 <- tsne_coords[, 1]
  beads_filt$tSNE_2 <- tsne_coords[, 2]
  
  
  # --- Clustering ---
  max_k_silhouette <- as.integer(max_k_silhouette)
  if (is.na(max_k_silhouette) || max_k_silhouette < 2) {
    max_k_silhouette <- 6L
  }
  max_k_feasible <- min(max_k_silhouette, n_samples - 1L)
  if (max_k_feasible < 2) {
    stop(sprintf("Too few spots (%d) for silhouette clustering (need at least 3).", n_samples))
  }
  ks <- seq.int(2L, max_k_feasible)

  sil_scores <- vapply(ks, function(k) {
    km <- kmeans(tsne_coords, centers = k, nstart = 25)
    mean(silhouette(km$cluster, dist(tsne_coords))[, 3])
  }, numeric(1))

  sil_df <- data.frame(K = ks, Silhouette = sil_scores)
  p_sil <- ggplot(sil_df, aes(x = K, y = Silhouette)) +
    geom_line(color = "#377EB8", linewidth = 1) +
    geom_point(color = "#E41A1C", size = 1.2) +
    theme_minimal() +
    labs(
      title = "Optimal K Selection via Silhouette Score",
      x = "Number of Clusters (K)",
      y = "Average Silhouette Score"
    ) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(
    file.path(result_dir, paste0(expt_id, "_Silhouette_Scores.pdf")),
    p_sil, width = 6, height = 4
  )

  if (!is.null(k_clusters)) {
    best_k <- as.integer(k_clusters)
    if (is.na(best_k) || best_k < 2) {
      stop("k_clusters must be an integer >= 2.")
    }
    cat(sprintf("   -> User specified Clusters: %d\n", best_k))
  } else {
    best_k <- ks[which.max(sil_scores)]
    cat(sprintf("   -> Optimal Clusters (Auto-detected): %d\n", best_k))
  }

  # final K-means clustering
  set.seed(123) # fixed seed for reproducible clustering
  final_km <- kmeans(tsne_coords, centers = best_k, nstart = 25)
  cluster_raw <- as.character(final_km$cluster)
  beads_filt$cluster <- as.factor(cluster_raw)

  # --- Spatial cluster plots (point map with touching spots) ---
  plot_dim <- spatial_plot_dim_inch(resolution)

  p1 <- plot_spatial_cns(beads_filt, resolution, "cluster", "Spatial Clusters", discrete = TRUE)
  ggsave(
    file.path(result_dir, paste0(expt_id, "_Spatial_Clusters.pdf")),
    p1, width = plot_dim, height = plot_dim
  )
  # --- CNV branch B: unsmoothed counts for local CNV ---
  cat("B. CNV normalization ...\n")
  bins <- bins_all[sel_bins, , drop = FALSE]
  norm_bulk_baseline <- norm_bulk_baseline[sel_bins]

  # Use unsmoothed filtered counts for CNV normalization (align with reference).
  counts_cnv_input <- as.matrix(counts_filt_all[, sel_bins, drop = FALSE])
  counts_norm <- counts_cnv_input
  
  cat("LOWESS normalization: GC -> Mappability...\n")
  chr_num <- suppressWarnings(as.numeric(gsub("^chr", "", bins$chr)))
  sel_bins_auto <- !is.na(chr_num) & chr_num <= 22
  gc_track_sel <- bins$gc_ratio
  map_track_sel <- bins$mappability_score
  
  # fit on selected beads + selected autosomal bins, then divide by fitted track trend.
  sel_beads_lowess <- rowSums(counts_cnv_input, na.rm = TRUE) > 0
  step_gc <- lowess_norm_suit_style(counts_cnv_input, gc_track_sel, sel_bins_auto, sel_beads_lowess, f = 0.1)
  step_map <- lowess_norm_suit_style(step_gc$counts_norm, map_track_sel, sel_bins_auto, sel_beads_lowess, f = 0.1)
  counts_norm <- step_map$counts_norm
  
  
  # Build bulk baseline normalization.
  bulk_sel <- as.numeric(norm_bulk_vec[sel_bins])
  bulk_mat <- matrix(bulk_sel, nrow = 1)
  bulk_norm <- bulk_sel

  bulk_gc <- lowess_norm_suit_style(bulk_mat, gc_track_sel, sel_bins_auto, TRUE, f = 0.1)
  bulk_map <- lowess_norm_suit_style(bulk_gc$counts_norm, map_track_sel, sel_bins_auto, TRUE, f = 0.1)
  bulk_norm <- as.numeric(bulk_map$counts_norm)
  
  bulk_norm[!is.finite(bulk_norm) | bulk_norm <= 0] <- median(bulk_norm[is.finite(bulk_norm) & bulk_norm > 0], na.rm = TRUE)

  # Convert back to bead-comparable scale for track plotting.
  # Direct lowess_norm output can be very small per bead, so we depth-normalize per bead
  # before comparing to the bulk baseline.
  counts_norm_depth <- counts_norm / pmax(rowMeans(counts_norm, na.rm = TRUE), 1e-6)
  bulk_scale <- median(bulk_norm[is.finite(bulk_norm) & bulk_norm > 0], na.rm = TRUE)
  if (!is.finite(bulk_scale) || bulk_scale <= 0) bulk_scale <- 1
  bulk_ref <- bulk_norm / bulk_scale
  bulk_ref[!is.finite(bulk_ref) | bulk_ref <= 0] <- 1
  
  # Apply sample/bulk ratio after both are in the same normalization space.
  counts_cnv <- sweep(counts_norm_depth, 2, bulk_ref, "/")
  
  # visualize_coverage(..., 2, "mode") effectively recenters each profile
  # around a diploid baseline. Without this step, ratios can collapse to near-zero.
  chr_num_plot <- suppressWarnings(as.numeric(gsub("^chr", "", bins$chr)))
  sel_bins_auto_plot <- !is.na(chr_num_plot) & chr_num_plot <= 22
  if (sum(sel_bins_auto_plot) > 10) {
    bead_baseline <- apply(counts_cnv[, sel_bins_auto_plot, drop = FALSE], 1, median, na.rm = TRUE)
    bead_baseline[!is.finite(bead_baseline) | bead_baseline <= 0] <- 1
    counts_cnv <- counts_cnv / bead_baseline
  }
  
  counts_log2 <- log2(counts_cnv + 1e-6)
  cat(sprintf("   -> CNV ranges: counts_norm_depth [%.4f, %.4f], bulk_ref [%.4f, %.4f], log2 [%.4f, %.4f]\n",
              min(counts_norm_depth, na.rm = TRUE), max(counts_norm_depth, na.rm = TRUE),
              min(bulk_ref, na.rm = TRUE), max(bulk_ref, na.rm = TRUE),
              min(counts_log2, na.rm = TRUE), max(counts_log2, na.rm = TRUE)))
  
  counts_log2[counts_log2 > 3] <- 3
  counts_log2[counts_log2 < -3] <- -3

  # ---- 4. Genome Tracks & Intra-Sample Correlation ----
  cat("4. Generating Diffused CNV Tracks & Pearson Correlation Labels...\n")
  # exclude X/Y; analyze autosomes only
  valid_idx <- which(!bins$chr_clean %in% c("X", "Y"))
  bins_plot <- bins[valid_idx, ]
  counts_log2_plot <- counts_log2[, valid_idx]
  
  # recalculate genome-axis offsets without X/Y
  bins_grouped_p <- aggregate(bin_end ~ chr_factor, data = bins_plot, FUN = max)
  colnames(bins_grouped_p)[colnames(bins_grouped_p) == "bin_end"] <- "max_pos"
  bins_grouped_p <- bins_grouped_p[order(bins_grouped_p$chr_factor), , drop = FALSE]
  bins_grouped_p$chr_shift <- c(0, head(cumsum(as.numeric(bins_grouped_p$max_pos)), -1))

  if ("chr_shift" %in% colnames(bins_plot)) bins_plot$chr_shift <- NULL
  bins_plot <- merge(bins_plot, bins_grouped_p[, c("chr_factor", "chr_shift"), drop = FALSE], by = "chr_factor", all.x = TRUE, sort = FALSE)
  bins_plot$cum_start <- bins_plot$bin_start + bins_plot$chr_shift
  bins_plot <- bins_plot[order(bins_plot$chr_factor, bins_plot$bin_start), ]
  
  x_axis_breaks <- (bins_grouped_p$chr_shift + (bins_grouped_p$max_pos / 2))
  x_axis_labels <- bins_grouped_p$chr_factor
  
  # cluster colors
  # optional external cluster color mapping
  cluster_colors <- c( "#4A90E2", "#7ED321", "#FF9500", "#E41A1C", "#377EB8", "#4DAF4A", "#F3E5E2")
  unique_clusters <- sort(unique(beads_filt$cluster))
  plot_data_list <- list()
  
  for (i in seq_along(unique_clusters)) {
    cid <- as.character(unique_clusters[i])
    cells_idx <- which(beads_filt$cluster == cid)
    
    # per-cluster spatial mean
    spatial_mean <- colMeans(counts_log2_plot[cells_idx, , drop = FALSE])

    # point plot of cluster means (avoids dense line artifacts)
    df_cluster_mean <- data.frame(
      x_pos = bins_plot$cum_start,
      log2_val = spatial_mean,
      cluster_id = cid
    )
    df_cluster_mean <- df_cluster_mean[df_cluster_mean$log2_val >= -1.2, ]
    df_cluster_mean$copy_number <- pmin(pmax(2 * (2 ^ df_cluster_mean$log2_val), 0), 6)
    
    # ggplot track: cluster mean points
    title_text <- paste0("Cluster ", cid, " (n=", length(cells_idx), ")")
    p_track <- ggplot() +
      # cluster track as points (not diffused lines)
      geom_point(data = df_cluster_mean, aes(x = x_pos, y = log2_val), color = cluster_colors[i],
        alpha = 0.9, size = 1.5, shape = 16, stroke = 0) +
      scale_x_continuous(breaks = x_axis_breaks, labels = x_axis_labels, expand = c(0,0)) +
      scale_y_continuous(limits = c(-2, 2), oob = scales::squish) +
      labs(title = title_text, y = "Log2 Ratio", x = NULL) +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 10, face = "bold")
      )
    
    p_track_cn <- ggplot() +
      geom_point(data = df_cluster_mean, aes(x = x_pos, y = copy_number), color = cluster_colors[i],
        alpha = 0.9, size = 1.5, shape = 16, stroke = 0) +
      scale_x_continuous(breaks = x_axis_breaks, labels = x_axis_labels, expand = c(0,0)) +
      scale_y_continuous(limits = c(0, 6), oob = scales::squish) +
      labs(title = paste0(title_text, " | Copy Number"), y = "Copy Number", x = NULL) +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 10, face = "bold")
      )
    # store for patchwork assembly
    plot_data_list[[paste0("C", cid)]] <- list(
      plot = p_track,
      plot_cn = p_track_cn,
      mean = spatial_mean,
      mean_points = df_cluster_mean
    )
  }
  
  # 1B. Add one merged track across all clusters (keep cluster colors)
  all_cells_idx <- seq_len(nrow(counts_log2_plot))
  merged_diffused_long <- dplyr::bind_rows(lapply(plot_data_list, function(x) x$mean_points))
  # Randomized display (fixed seed): reduce deterministic over-plotting in merged view.
  merged_show_seed <- 2026
  merged_show_frac <- 0.8
  if (nrow(merged_diffused_long) > 0) {
    set.seed(merged_show_seed)
    keep_n <- max(1, floor(nrow(merged_diffused_long) * merged_show_frac))
    keep_idx <- sample.int(nrow(merged_diffused_long), size = keep_n, replace = FALSE)
    merged_diffused_long <- merged_diffused_long[keep_idx, , drop = FALSE]
    merged_diffused_long <- merged_diffused_long[sample.int(nrow(merged_diffused_long)), , drop = FALSE]
  }
  # Add a small per-cluster x-offset to avoid complete point overlap in merged view.
  cluster_offsets <- setNames(
    seq_along(unique_clusters) - (length(unique_clusters) + 1) / 2,
    as.character(unique_clusters)
  )
  merged_diffused_long$x_pos_shifted <- merged_diffused_long$x_pos +
    as.numeric(cluster_offsets[as.character(merged_diffused_long$cluster_id)]) * 120000
  merged_diffused_long$copy_number <- pmin(pmax(2 * (2 ^ merged_diffused_long$log2_val), 0), 6)
  cluster_color_map <- setNames(cluster_colors[seq_along(unique_clusters)], as.character(unique_clusters))
  
  merged_track_plot <- ggplot() +
    geom_point(data = merged_diffused_long, aes(x = x_pos_shifted, y = log2_val, color = cluster_id),
      alpha = 1, size = 1.5, shape = 16, stroke = 0) +
    scale_color_manual(values = cluster_color_map, guide = "none") +
    scale_x_continuous(breaks = x_axis_breaks, labels = x_axis_labels, expand = c(0, 0)) +
    scale_y_continuous(limits = c(-2, 2), oob = scales::squish) +
    labs(title = paste0("Merged Clusters (n=", length(all_cells_idx), ")"), y = "Log2 Ratio", x = NULL) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = 10, face = "bold")
    )
  
  # show x-axis labels on the last track
  if (length(plot_data_list) > 0) {
    last_idx <- length(plot_data_list)
    plot_data_list[[last_idx]]$plot <- plot_data_list[[last_idx]]$plot + 
      theme(axis.text.x = element_text(size = 8, angle = 0))
  }
  

  # 3. Assemble with Patchwork
  cat("   -> Assembling Final PDF...\n")
  all_grobs <- c(list(merged_track_plot), lapply(plot_data_list, function(x) x$plot))
  
  final_design <- wrap_plots(all_grobs, ncol = 1)
  ggsave(file.path(result_dir, paste0(expt_id, ".Genome_Tracks_Overlap.pdf")), final_design, width = 12, 
         height = 2.5 * length(plot_data_list)+2.5 )
  

  
  
  
  # ---- 6. Multi-omics Alignment (If RNA exists) ----
  if (!is.null(rna_loc_path) && file.exists(rna_loc_path)) {
    require_optional_package("Morpho")
    cat("6. Spatial Multi-omics Registration (ICP)...\n")
    dat_rna <- fread(rna_loc_path, header = TRUE) %>% filter(tissue == 1 & !is.na(row) & !is.na(col))
    dat_rna$col <- -dat_rna$col 
    
    dat_dna <- data.frame(DNA_Barcode = beads_filt[[1]], xcoord = beads_filt$xcoord, ycoord = beads_filt$ycoord, cluster = beads_filt$cluster)
    
    # Global Joint Scaling
    all_x <- c(dat_dna$xcoord, dat_rna$row); all_y <- c(dat_dna$ycoord, dat_rna$col)
    max_side <- max(max(all_x) - min(all_x), max(all_y) - min(all_y))
    
    dat_dna$x_norm <- (dat_dna$xcoord - min(all_x)) / max_side
    dat_dna$y_norm <- (dat_dna$ycoord - min(all_y)) / max_side
    dat_rna$x_norm <- (dat_rna$row - min(all_x)) / max_side
    dat_rna$y_norm <- (dat_rna$col - min(all_y)) / max_side
    
    # ICP and Label Transfer
    icp_res <- icpmat(x = as.matrix(dat_rna[, c("x_norm", "y_norm")]), y = as.matrix(dat_dna[, c("x_norm", "y_norm")]), iterations = 100)
    knn_result <- FNN::get.knnx(data = as.matrix(dat_dna[, c("x_norm", "y_norm")]), query = icp_res, k = 1)
    
    dat_rna$Mapped_DNA_Cluster <- as.factor(dat_dna$cluster[knn_result$nn.index[, 1]])
    dat_rna$Display_Cluster <- ifelse(knn_result$nn.dist[, 1] > 0.03, "Unmapped", as.character(dat_rna$Mapped_DNA_Cluster))
    dat_rna$aligned_x_norm <- icp_res[,1]
    dat_rna$aligned_y_norm <- icp_res[,2]
  }
  
  cat("Preparing final data objects...\n")
  cluster_means <- lapply(unique_clusters, function(cid) {
    colMeans(counts_log2_plot[beads_filt$cluster == cid, , drop = FALSE], na.rm = TRUE)
  })
  names(cluster_means) <- unique_clusters
  counts_plot <- counts_log2
  
  res_list <- list(
    expt_id = expt_id,
    result_dir = result_dir,
    beads_filt = beads_filt,
    cluster_means = cluster_means,
    counts_plot = counts_plot,
    bins_plot = bins_plot                # required for difference maps and trees
  )
  
  # save RDS for downstream functions
  saveRDS(res_list, file = file.path(result_dir, paste0(expt_id, "_spatial_res.rds")))
  
  cat(sprintf("=== Pipeline Complete for %s ===\n", expt_id))
  
  return(res_list)
}




# ===================== helper functions ===================

# load saved spatial DNA result by expt_id
load_spatial_dna_res <- function(expt_id, base_dir) {
  rds_path <- file.path(base_dir, "results/Spatial_DNA", expt_id, paste0(expt_id, "_spatial_res.rds"))
  if (!file.exists(rds_path)) {
    stop(sprintf("Result file not found: %s\nRun run_unified_spatial_dna() first.", rds_path))
  }
  return(readRDS(rds_path))
}




# ===================== 3. Inter-Sample Subclone Correlation ===================
compare_all_clusters_comprehensive <- function(expt_id_a, expt_id_b, base_dir, gain_loss_cutoff = 0.25) {
  cat(sprintf("\n=== Comprehensive clone comparison: %s vs %s ===\n", expt_id_a, expt_id_b))
  
  # load from disk
  res_a <- load_spatial_dna_res(expt_id_a, base_dir)
  res_b <- load_spatial_dna_res(expt_id_b, base_dir)
  
  means_a <- res_a$cluster_means
  means_b <- res_b$cluster_means
  clusters_a <- names(means_a)
  clusters_b <- names(means_b)
  
  results_list <- list()
  
  for (i in seq_along(clusters_a)) {
    for (j in seq_along(clusters_b)) {
      ca <- clusters_a[i]
      cb <- clusters_b[j]
      
      profile_a <- as.numeric(means_a[[ca]])
      profile_b <- as.numeric(means_b[[cb]])
      
      spearman_rho <- cor(profile_a, profile_b, method = "spearman", use = "complete.obs")
      pearson_r <- cor(profile_a, profile_b, method = "pearson", use = "complete.obs")
      
      state_a <- ifelse(profile_a > gain_loss_cutoff, 1L, ifelse(profile_a < -gain_loss_cutoff, -1L, 0L))
      state_b <- ifelse(profile_b > gain_loss_cutoff, 1L, ifelse(profile_b < -gain_loss_cutoff, -1L, 0L))
      
      union_event <- sum((state_a != 0L) | (state_b != 0L))
      intersect_event <- sum((state_a == state_b) & (state_a != 0L))
      jaccard_val <- ifelse(union_event == 0, NA_real_, intersect_event / union_event)
      
      results_list[[paste0(ca, "_vs_", cb)]] <- data.frame(
        Sample_A_Cluster = paste0(res_a$expt_id, "_C", ca),
        Sample_B_Cluster = paste0(res_b$expt_id, "_C", cb),
        Spearman = spearman_rho,
        Pearson = pearson_r,
        Jaccard = jaccard_val
      )
    }
  }
  
  final_stats <- do.call(rbind, results_list)
  rownames(final_stats) <- NULL
  
  top_match <- final_stats[which.max(final_stats$Spearman), ]
  cat(sprintf("   -> Best evolutionary match: %s and %s\n", top_match$Sample_A_Cluster, top_match$Sample_B_Cluster))
  cat(sprintf("      Spearman: %.3f | Pearson: %.3f | Jaccard: %.3f\n", 
              top_match$Spearman, top_match$Pearson, top_match$Jaccard))
  
  write.csv(final_stats, file.path(res_a$result_dir, paste0("Comprehensive_Stats_", res_a$expt_id, "_vs_", res_b$expt_id, ".csv")), row.names = FALSE)
  
  return(final_stats)
}



# ===================== 4. RNA pseudotime (Monocle2) ===================
run_rna_pseudotime_monocle2 <- function(seurat_obj, output_dir, expt_id, cluster_col = "seurat_clusters", root_cluster = NULL) {
  cat(sprintf("\n=== Building RNA pseudotime trajectory: %s ===\n", expt_id))
  
  expr_matrix <- as(as.matrix(seurat_obj@assays$RNA@counts), 'sparseMatrix')
  cell_ann <- seurat_obj@meta.data
  gene_ann <- data.frame(gene_short_name = rownames(expr_matrix), row.names = rownames(expr_matrix))
  
  pd <- new("AnnotatedDataFrame", data = cell_ann)
  fd <- new("AnnotatedDataFrame", data = gene_ann)
  
  cds <- newCellDataSet(expr_matrix,
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())
  
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  # variable gene selection
  var_genes <- if ("SCT" %in% names(seurat_obj@assays)) VariableFeatures(seurat_obj, assay = "SCT") else VariableFeatures(seurat_obj)
  if (length(var_genes) == 0) stop("No VariableFeatures found in Seurat object.")
  
  cds <- setOrderingFilter(cds, var_genes)
  
  cat("   -> Dimensionality reduction (DDRTree)...\n")
  cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
  
  cat("   -> Initial ordering...\n")
  cds <- orderCells(cds)
  
  # optional manual root cluster
  if (!is.null(root_cluster)) {
    if (root_cluster %in% pData(cds)[[cluster_col]]) {
      cat(sprintf("   -> Setting RNA root cluster: %s\n", root_cluster))
      # dominant state for the root cluster
      state_table <- table(pData(cds)$State[pData(cds)[[cluster_col]] == root_cluster])
      root_state <- names(which.max(state_table))
      
      cds <- orderCells(cds, root_state = root_state)
      cat(sprintf("   -> RNA root state set to %s.\n", root_state))
    } else {
      warning("   -> root_cluster not found; keeping default ordering.")
    }
  }
  
  # visualization
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  p1 <- plot_cell_trajectory(cds, color_by = cluster_col, cell_size = 1.2) +
    labs(title = paste("RNA Trajectory (Root:", ifelse(is.null(root_cluster), "Auto", root_cluster), ")"))
  p2 <- plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 1.2) +
    scale_color_viridis_c(option = "mako")
  
  ggsave(file.path(output_dir, paste0(expt_id, "_RNA_Pseudotime_Calibrated.pdf")), p1 | p2, width = 12, height = 5)
  
  return(cds)
}



# ===================== 5. MEDICC2 export ===================================
# export spatial DNA CNV profiles to MEDICC2 input format
export_for_medicc2 <- function(res_list, out_dir = out_dir, project_name = "Chemo_Evolution") {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("\n=== Exporting MEDICC2 input: %s ===\n", project_name))
  
  medicc_df_list <- list()
  
  # iterate over sample results
  for (res in res_list) {
    
    
    bins <- res$bins_plot
    means <- res$cluster_means
    
    for (c_name in names(means)) {
      # MEDICC2 expects integer absolute copy numbers
      # diploid baseline: absolute CN = 2 * 2^log2_ratio
      # approximate integer CN for spatial data without absolute WGS quantification
      log2_val <- as.numeric(means[[c_name]])
      abs_cn <- round(2 * (2 ^ log2_val))
      
      # cap CN to keep MEDICC2 runtime reasonable
      abs_cn[abs_cn > 8] <- 8 
      abs_cn[abs_cn < 0] <- 0
      
      temp_df <- data.frame(
        sample_id = paste0(res$expt_id, "_C", c_name),
        chrom = bins$chr_clean,
        start = bins$bin_start,
        end = bins$bin_end,
        cn = abs_cn
      )
      medicc_df_list[[length(medicc_df_list) + 1]] <- temp_df
    }
  }
  
  final_medicc_df <- do.call(rbind, medicc_df_list)
  
  # add diploid normal root for MEDICC2
  normal_df <- data.frame(
    sample_id = "Diploid_Root",
    chrom = bins$chr_clean,
    start = bins$bin_start,
    end = bins$bin_end,
    cn = 2
  )
  final_medicc_df <- rbind(normal_df, final_medicc_df)
  
  out_file <- file.path(out_dir, paste0(project_name, "_medicc2_input.tsv"))
  write.table(final_medicc_df, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat(sprintf("   -> MEDICC2 input saved: %s\n", out_file))
  cat("   -> Run MEDICC2 on the server:\n")
  cat(sprintf(
    "      medicc2 %s ./results/MEDICC2_Output/ -j 8 --plot both --total-copy-numbers --normal-name Diploid_Root\n",
    out_file
  ))
  
  return(final_medicc_df)
}


# ==============================================================================
# Usage examples (set RUN_EXAMPLES <- TRUE to execute)
# ==============================================================================

RUN_EXAMPLES <- FALSE

if (isTRUE(RUN_EXAMPLES)) {
  repo_root <- Sys.getenv("AXIS_REPO_ROOT", unset = normalizePath("../..", winslash = "/"))
  base_dir <- repo_root
  resolution <- DEFAULT_SPATIAL_RESOLUTION

  bins_path <- file.path(base_dir, "reference/genomic_bins/GRCh38_1Mb_bins.txt")
  gc_path <- file.path(base_dir, "reference/gc_content/GRCh38_1Mb_gc.txt")
  map_path <- file.path(base_dir, "reference/mappability/hg38_1Mb_map.txt")

  expt_id <- "17_1"
  bulk_control_path <- file.path(base_dir, "data/bulk/human_CRC_17_1Mb_total_reads.txt")

  res_spatial_dna <- run_unified_spatial_dna(
    expt_id = expt_id,
    base_dir = base_dir,
    resolution = resolution,
    bulk_control_path = bulk_control_path,
    bins_path = bins_path,
    gc_path = gc_path,
    map_path = map_path
  )

  # MEDICC2 export example
  # pre_ids <- c("17_1", "19_1", "19_2", "19_3", "20_1", "20_2", "20_3")
  # my_res_list <- lapply(pre_ids, load_spatial_dna_res, base_dir = base_dir)
  # export_for_medicc2(my_res_list, file.path(base_dir, "results/MEDICC2_Input"), "Chemo_Evolution")
  # Then run in shell:
  # medicc2 results/MEDICC2_Input/Chemo_Evolution_medicc2_input.tsv results/MEDICC2_Output/ -j 8 --plot both --total-copy-numbers --normal-name Diploid_Root
}

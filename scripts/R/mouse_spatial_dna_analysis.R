# AXIS-DNA-seq: Mouse spatial DNA CNV analysis (GRCm39)

# install.packages(c("data.table", "Matrix", "ggplot2", "dplyr", "FNN", "patchwork", "RColorBrewer", "viridis", "cluster", "Rtsne", "hexbin"))
rm(list = ls())
gc()

library(data.table)
library(Matrix)
library(ggplot2)
library(dplyr)
library(FNN)            
library(patchwork)      
library(cluster)
library(Rtsne)
library(viridis)

# Shared palettes to keep colors consistent
TRACK_CLUSTER_PALETTE <- c("#D19366", "#EBC6BE", "#DF9A69", "#89669D", "#84A59D", "#C57B86", "#6D98BA",
                           "#E2B475", "#9D8189", "#D19781", "#5D737E", "#CBB6CE")
CNV_REDBLUE_PALETTE <- c(low = "#2166AC", mid = "white", high = "#B2182B")

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

# ===================== 1. Helper Functions =========================================================

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
      axis.line = element_line(color = "black", linewidth = 0.1),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.text = element_text(color = "black", size = 10)
    ) +
    labs(title = title_str, color = "")
  
  if (discrete) {
    p <- p + scale_color_manual(values = c("#4A90E2", "#7ED321", "#FF9500", "#E41A1C", "#377EB8", "#4DAF4A", "#F3E5E2")) 
  } else {
    if (palette_type == "viridis") {
      p <- p + scale_color_viridis_c(option = "magma")
    } else if (palette_type == "redblue") {
      p <- p + scale_color_gradient2(low = CNV_REDBLUE_PALETTE["low"], mid = CNV_REDBLUE_PALETTE["mid"], high = CNV_REDBLUE_PALETTE["high"], midpoint = 0)
    }
  }
  return(p)
}

# ============= 2. Main Analysis Function =================================================================

run_slide_dna_mouse <- function(
    expt_id,
    base_dir = ".",
    resolution = DEFAULT_SPATIAL_RESOLUTION,
    normal_bulk_path = NULL,        # Path to Normal Mouse Bulk (Reference)
    gc_path = NULL,                 # Path to GC annotation file
    map_path = NULL,                # Path to Mappability annotation file
    k_clusters = NULL,              # Allow manual specification of cluster count 'k'
    return_details = FALSE,         # Return detailed objects for cross-slice comparison
    cluster_relabel = NULL          # Standardize cluster labels, e.g., c("1"="Tumor","2"="NonTumor")
) {
  
  # --- Path Setup ---
  spatial_dna_root <- file.path(base_dir, "data/Spatial_DNA")
  data_dir <- file.path(spatial_dna_root, expt_id)
  result_dir <- file.path(base_dir, "results/Spatial_DNA", expt_id)
  
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  plot_dim <- spatial_plot_dim_inch(resolution)
  
  # --- 1. Load Data ---
  cat("1. Loading Data...\n")
  
  # Load Sparse Counts & Beads
  sparse_path <- file.path(data_dir, paste0(expt_id, ".sparse_counts_1Mb.txt"))
  beads_path <- file.path(spatial_dna_root, paste0(expt_id, ".spatial_barcodes_location.csv"))
  sparse_dt <- fread(sparse_path, header = FALSE, col.names = c("row", "col", "val"))
  beads <- fread(beads_path)
  
  # Load GRCm39 Bins
  bins_path <- file.path(base_dir, "reference/genomic_bins/GRCm39_1Mb_bins.txt")
  bins <- fread(bins_path, header = TRUE) 
  
  # Ensure bins are sorted by chromosome and position
  # Mouse chromosomes: 1-19, X, Y
  chr_levels <- c(as.character(1:19), "X", "Y")
  bins$chr_clean <- gsub("chr", "", bins$chr)
  bins$chr_factor <- factor(bins$chr_clean, levels = chr_levels)
  bins <- bins[order(bins$chr_factor, bins$bin_start), ]
  
  # ================= Module 1: Load and align GC and Mappability annotations ==================
  if (!is.null(gc_path) && !is.null(map_path)) {
    cat("   -> Loading GC and Mappability annotations...\n")
    # Read GC (No header, column 8 is gc_ratio)
    gc_dt <- fread(gc_path, header = FALSE)
    gc_dt <- gc_dt[, c(1, 2, 3, 8), with = FALSE]
    colnames(gc_dt) <- c("chr", "bin_start", "bin_end", "gc_ratio")
    
    # Read Mappability (Has header)
    map_dt <- fread(map_path, header = TRUE)
    map_dt <- map_dt[, c("chr", "bin_start", "bin_end", "mappability_score"), with = FALSE]
    
    # Align to main bins
    bins <- left_join(bins, gc_dt, by = c("chr", "bin_start", "bin_end"))
    bins <- left_join(bins, map_dt, by = c("chr", "bin_start", "bin_end"))
    
    # Handle missing values: Impute GC with mean, default Map to 0 (untrusted)
    bins$gc_ratio[is.na(bins$gc_ratio)] <- mean(bins$gc_ratio, na.rm=TRUE)
    bins$mappability_score[is.na(bins$mappability_score)] <- 0
  } else {
    cat("   -> [Warning] GC or Map path not provided. Skipping correction.\n")
    bins$gc_ratio <- NA
    bins$mappability_score <- 1 # Default to perfect mappability
  }
  # ==============================================================================
  
  # Calculate Cumulative Position for Plotting (Manhattan Plot style)
  bins_grouped <- bins %>% group_by(chr_factor) %>% summarize(max_pos = max(bin_end), .groups="drop")
  bins_grouped$chr_shift <- c(0, cumsum(as.numeric(bins_grouped$max_pos))[1:(nrow(bins_grouped)-1)])
  bins <- merge(bins, bins_grouped[,c("chr_factor","chr_shift")], by="chr_factor", all.x=TRUE)
  bins$cum_start <- bins$bin_start + bins$chr_shift
  bins <- bins[order(bins$chr_factor, bins$bin_start), ] # Ensure order persists
  
  counts_mat <- sparseMatrix(i = sparse_dt$row, j = sparse_dt$col, x = sparse_dt$val)
  
  # --- Load Normal Bulk (Dynamic alignment implementation) ---
  if (is.null(normal_bulk_path)) stop("Normal Mouse Bulk path is required for CNV normalization.")
  cat("   -> Loading and aligning Normal Mouse Bulk data...\n")
  
  # 1. Read raw Bulk data (Assuming first 4 cols are chr, start, end, count)
  bulk_raw_df <- fread(normal_bulk_path, header = FALSE)
  colnames(bulk_raw_df)[1:4] <- c("bulk_chr", "bulk_start", "bulk_end", "bulk_count")
  
  # 2. Auto-align: Match bulk data to current reference bins via chromosome and end coordinate
  aligned_bulk <- dplyr::left_join(bins, bulk_raw_df, 
                                   by = c("chr" = "bulk_chr", "bin_end" = "bulk_end"))
  
  # 3. Extract aligned values, impute missing as 0
  norm_bulk_vec <- aligned_bulk$bulk_count
  norm_bulk_vec[is.na(norm_bulk_vec)] <- 0
  
  cat(sprintf("   -> Successfully aligned Bulk data to %d genomic bins!\n", length(norm_bulk_vec)))
  
  # 4. Normalization Processing
  bulk_mean <- mean(norm_bulk_vec[norm_bulk_vec > 0], na.rm = TRUE)
  if(is.na(bulk_mean) || bulk_mean == 0) bulk_mean <- 1 
  
  norm_bulk_baseline <- norm_bulk_vec / bulk_mean
  norm_bulk_baseline[norm_bulk_baseline <= 0] <- 1e-6
  
  # --- 2. Filter & Preprocessing ---
  cat("2. Preprocessing...\n")
  beads$ycoord <- -beads$ycoord
  
  # Filter Coverage
  keep_cells <- rowSums(counts_mat) > 100
  beads_filt <- beads[keep_cells, ]
  counts_filt <- counts_mat[keep_cells, ]
  coords <- as.matrix(beads_filt[, .(xcoord, ycoord)])
  
  # --- 3. Smoothing & Normalization ---
  cat("3. Smoothing & Normalizing...\n")
  
  # KNN Smoothing (XY)
  counts_smo_xy <- smooth_matrix_knn(as.matrix(counts_filt), coords, k=50)
  
  # Normalize by Cell Depth
  counts_norm_depth <- counts_smo_xy / rowMeans(counts_smo_xy)
  
  # Normalize against Normal Mouse Bulk (CNV Calculation)
  counts_cnv <- sweep(counts_norm_depth, 2, norm_bulk_baseline, "/")
  counts_log2 <- log2(counts_cnv + 1e-6)
  
  # ================= Module 2: Core logic for GC Bias & Mappability correction =================
  if (!is.null(gc_path) && !is.null(map_path)) {
    cat("3.5 Running GC Bias Correction (LOESS) & Mappability Filtering...\n")
    
    # Filter reliable bins (Mappability >= 0.8)
    valid_bins_idx <- which(bins$mappability_score >= 0.8)
    invalid_bins_idx <- which(bins$mappability_score < 0.8)
    
    gc_vals <- bins$gc_ratio[valid_bins_idx]
    counts_log2_corrected <- counts_log2
    
    pb <- txtProgressBar(min = 0, max = nrow(counts_log2), style = 3)
    for(i in 1:nrow(counts_log2)) {
      cell_vals <- counts_log2[i, valid_bins_idx]
      
      # Fit GC bias curve for current cell using LOESS
      # tryCatch prevents failure on highly anomalous single cells
      tryCatch({
        loess_fit <- loess(cell_vals ~ gc_vals, span = 0.3)
        corrected_vals <- cell_vals - predict(loess_fit)
        counts_log2_corrected[i, valid_bins_idx] <- corrected_vals
      }, error = function(e) {})
      
      if(i %% 100 == 0) setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Force untrusted regions (low Mappability) to 0 to eliminate false positive spikes
    if (length(invalid_bins_idx) > 0) {
      counts_log2_corrected[, invalid_bins_idx] <- 0
    }
    
    counts_log2 <- counts_log2_corrected
  }
  # ==============================================================================
  
  # Cap extreme values for analysis stability
  counts_log2[counts_log2 > 3] <- 3
  counts_log2[counts_log2 < -3] <- -3
  
  # --- 4. Clustering (PCA -> tSNE -> Kmeans) ---
  cat("4. Clustering...\n")
  
  # Select variable bins
  bin_vars <- apply(counts_log2, 2, var)
  pc_bins <- which(bin_vars > quantile(bin_vars, 0.4)) # Top 60% variable bins
  
  pca_res <- prcomp(counts_log2[, pc_bins], center = TRUE, scale. = FALSE, rank. = 30)
  
  # Smooth PCA
  pc_dist <- pca_res$x[, 1:10]
  pc_smo <- smooth_matrix_knn(pca_res$x[, 1:10], pc_dist, k=50)
  pc_smo_both <- smooth_matrix_knn(pc_smo, coords, k=10)
  
  # tSNE
  set.seed(123) 
  tsne_res <- Rtsne(pc_smo_both, dims = 2, perplexity = 150, verbose = FALSE, check_duplicates=FALSE)
  beads_filt$tSNE_1 <- tsne_res$Y[, 1]
  beads_filt$tSNE_2 <- tsne_res$Y[, 2]
  
  # --- K-means Clustering (Manual specification or auto-detection) ---
  if (!is.null(k_clusters)) {
    # Scenario A: User specified 'k'
    best_k <- k_clusters
    cat(sprintf("   -> User specified Clusters: %d\n", best_k))
  } else {
    # Scenario B: Auto-detect optimal 'k' via Silhouette scoring
    ks <- 2:6
    sil_scores <- sapply(ks, function(k) {
      km <- kmeans(tsne_res$Y, centers = k, nstart = 25)
      mean(silhouette(km$cluster, dist(tsne_res$Y))[, 3])
    })
    best_k <- ks[which.max(sil_scores)]
    cat(sprintf("   -> Optimal Clusters (Auto-detected): %d\n", best_k))
    
    # Plot and save the Silhouette scores to visually verify the optimal K
    cat("   -> Saving Optimal K Silhouette plot...\n")
    sil_df <- data.frame(k = ks, Silhouette_Score = sil_scores)
    p_sil <- ggplot(sil_df, aes(x = k, y = Silhouette_Score)) +
      geom_line(color = "#2166AC", linewidth = 1) +
      geom_point(size = 1.2, color = "#B2182B") +
      theme_bw() +
      labs(title = "Optimal K Selection", x = "Number of Clusters (k)", y = "Mean Silhouette Score") +
      scale_x_continuous(breaks = ks)
    ggsave(file.path(result_dir, paste0(expt_id, "_Silhouette_Scores.pdf")), p_sil, width = 6, height = 4)
    
    
    beads_filt_test <- beads_filt 
    for(i in 2:6){
      set.seed(123) # fixed seed for reproducible clustering
      test_km <- kmeans(tsne_res$Y, centers = i, nstart = 25)
      cluster_test <- as.character(test_km$cluster)
      beads_filt_test$cluster <- as.factor(cluster_test)
      p2 <- plot_spatial_cns(beads_filt_test, resolution, "cluster", "Spatial Clusters", discrete=TRUE)
      ggsave(file.path(result_dir, paste0(expt_id,"_Spatial_Clusters_k",i,".pdf")), p2, width = plot_dim, height = plot_dim)
    }
   
  }
  
  # Execute final K-means
  set.seed(123) # Ensure reproducibility across identical runs
  final_km <- kmeans(tsne_res$Y, centers = best_k, nstart = 25)
  cluster_raw <- as.character(final_km$cluster)
  
  if (!is.null(cluster_relabel)) {
    if (is.null(names(cluster_relabel))) {
      stop("cluster_relabel must be a named vector, e.g., c('1'='Tumor', '2'='NonTumor').")
    }
    cluster_new <- ifelse(
      cluster_raw %in% names(cluster_relabel),
      unname(cluster_relabel[cluster_raw]),
      cluster_raw
    )
    # Keep order stable across samples when using the same relabel dictionary
    beads_filt$cluster <- factor(cluster_new, levels = unique(unname(cluster_relabel)))
    cat("   -> Applied cluster relabel mapping:\n")
    print(cluster_relabel)
  } else {
    beads_filt$cluster <- as.factor(cluster_raw)
  }
  
  # --- Figures A & B (Spatial) ---
  p1 <- plot_spatial_cns(beads_filt, resolution, "cluster", "Spatial Clusters", discrete=TRUE)
  p2 <- plot_spatial_cns(beads_filt, resolution, "tSNE_1", "tSNE Dimension 1", "redblue")
  ggsave(file.path(result_dir, paste0(expt_id,"_Spatial_Clusters.pdf")), p1, width = plot_dim, height = plot_dim)

  
  # --- Figure D: Genome Tracks ---
  cat("5. Generating Genome Profile Tracks (Figure D)...\n")
  
  # 1. Exclude X and Y chromosomes, analyze autosomes only
  cat("   -> Excluding Chromosomes X and Y for cleaner correlation...\n")
  valid_idx <- which(!bins$chr_clean %in% c("X", "Y"))
  bins_plot <- bins[valid_idx, ]
  counts_log2_plot <- counts_log2[, valid_idx]
  
  # 2. Recalculate X-axis coordinate shifts
  bins_grouped_p <- bins_plot %>% 
    group_by(chr_factor) %>% 
    summarize(max_pos = max(bin_end), .groups="drop")
  bins_grouped_p$chr_shift <- c(0, cumsum(as.numeric(bins_grouped_p$max_pos))[1:(nrow(bins_grouped_p)-1)])
  
  bins_plot <- bins_plot %>% 
    select(-chr_shift) %>% 
    left_join(bins_grouped_p[,c("chr_factor","chr_shift")], by="chr_factor")
  bins_plot$cum_start <- bins_plot$bin_start + bins_plot$chr_shift
  bins_plot <- bins_plot[order(bins_plot$chr_factor, bins_plot$bin_start), ]
  
  x_axis_breaks <- (bins_grouped_p$chr_shift + (bins_grouped_p$max_pos / 2))
  x_axis_labels <- bins_grouped_p$chr_factor
  
  # 3. Color definitions
  cluster_colors <- TRACK_CLUSTER_PALETTE
  
  unique_clusters <- sort(unique(beads_filt$cluster))
  plot_data_list <- list()
  
  for (i in seq_along(unique_clusters)) {
    cid <- as.character(unique_clusters[i])
    cells_idx <- which(beads_filt$cluster == cid)
    
    # Calculate spatial statistics
    spatial_mean <- colMeans(counts_log2_plot[cells_idx, ])
    spatial_sd <- apply(counts_log2_plot[cells_idx, ], 2, sd)
    
    # Construct plotting DataFrame
    df_plot <- data.frame(
      x_pos = bins_plot$cum_start,
      y_mean = spatial_mean,
      y_min = spatial_mean - spatial_sd,
      y_max = spatial_mean + spatial_sd,
      cluster = paste0("Cluster ", cid)
    )
    
    # Filter points where y < -1.2 to remove spikes
    df_plot <- df_plot[df_plot$y_mean >= -1.2, ]
    
    title_text <- paste0("Cluster ", cid, " (n=", length(cells_idx), ")")
    
    # =================================================================================
    # [Core] Implement diffused single-bead profile tracks
    # =================================================================================
    
    # Step A: Randomly sample up to 100 beads from the current cluster for the 'diffused' effect
    sample_size <- min(100, length(cells_idx))
    sample_beads <- sample(cells_idx, sample_size)
    
    # Prepare multi-line plotting data (Long Format)
    df_diffused <- as.data.frame(t(counts_log2_plot[sample_beads, , drop=FALSE]))
    colnames(df_diffused) <- paste0("bead_", 1:sample_size)
    df_diffused$x_pos <- bins_plot$cum_start
    
    # Convert to long format
    df_diffused_long <- tidyr::pivot_longer(df_diffused, cols = -x_pos, 
                                            names_to = "bead_id", values_to = "log2_val")
    
    # Apply y >= -1.2 filter again to ensure diffused lines have no spikes
    df_diffused_long <- df_diffused_long[df_diffused_long$log2_val >= -1.2, ]
    
    # Step B: Reconstruct ggplot logic
    p_track <- ggplot() +
      # Core: Draw many semi-transparent thin lines to create a 'cloud' effect
      geom_line(data = df_diffused_long, aes(x = x_pos, y = log2_val, group = bead_id), 
                color = cluster_colors[i], alpha = 0.1, linewidth = 0.2) +
      # Styling
      scale_x_continuous(breaks = x_axis_breaks, labels = x_axis_labels, expand = c(0,0)) +
      scale_y_continuous(limits = c(-1.5, 1.5), oob = scales::squish) +
      labs(title = title_text, y = "Log2 Ratio", x = NULL) +
      theme_bw() +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 10, face = "bold")
      )
    # =================================================================================
    
    # Store in list for downstream assembly
    plot_data_list[[paste0("C", cid)]] <- list(plot = p_track, mean = spatial_mean, df = df_plot)
  }
  
  # Enable X-axis text on the very last plot
  if (length(plot_data_list) > 0) {
    last_idx <- length(plot_data_list)
    plot_data_list[[last_idx]]$plot <- plot_data_list[[last_idx]]$plot + 
      theme(axis.text.x = element_text(size = 8, angle = 0))
  }
  
  # --- 2. Generate Difference Tracks ---
  cat("   -> Generating Difference Tracks (Autosomes Only)...\n")
  diff_plots <- list()
  
  n_clus <- length(unique_clusters)
  if (n_clus > 1) {
    for (k in 1:(n_clus - 1)) {
      c1_name <- paste0("C", unique_clusters[k])
      c2_name <- paste0("C", unique_clusters[k+1])
      
      mean1 <- plot_data_list[[c1_name]]$mean
      mean2 <- plot_data_list[[c2_name]]$mean
      diff_val <- mean1 - mean2
      
      # Ensure row counts match using bins_plot
      df_diff <- data.frame(
        x_pos = bins_plot$cum_start, 
        val = diff_val
      )
      
      p_diff <- ggplot(df_diff, aes(x = x_pos, y = val)) +
        geom_hline(yintercept = 0, color = "grey50") +
        geom_line(color = "grey30", linewidth = 0.4) +
        geom_ribbon(aes(ymin=0, ymax=val, fill = val > 0), alpha=0.6) +
        scale_fill_manual(values = c("TRUE" = cluster_colors[k], "FALSE" = cluster_colors[k+1]), guide="none") +
        scale_x_continuous(breaks = x_axis_breaks, labels = x_axis_labels, expand = c(0,0)) +
        scale_y_continuous(limits = c(-1, 1), oob = scales::squish, breaks=c(-1, 0, 1)) +
        labs(y = "Diff", title = paste0("Difference (", c1_name, " - ", c2_name, ")")) +
        theme_bw() +
        theme(
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
          axis.text.x = if(k == (n_clus-1)) element_text(size=8, angle=0) else element_blank(),
          plot.title = element_text(size = 9, face="italic")
        )
      
      diff_plots[[k]] <- p_diff
    }
  }
  
  # --- 3. Assemble with Patchwork ---
  cat("   -> Assembling Final PDF...\n")
  
  all_grobs <- c(lapply(plot_data_list, function(x) x$plot), diff_plots)
  final_design <- wrap_plots(all_grobs, ncol = 1)
  
  # Calculate dynamic height
  pdf_h <- 2.5 * length(plot_data_list) + 1.5 * length(diff_plots)
  
  save_path <- file.path(result_dir, paste0(expt_id, ".Genome_Tracks_Overlap.pdf"))
  ggsave(save_path, final_design, width = 12, height = pdf_h)
  
  cat("=== Analysis Complete ===\n")
  
  if (!isTRUE(return_details)) {
    return(beads_filt)
  }
  
  cluster_means <- lapply(as.character(unique_clusters), function(cid) {
    plot_data_list[[paste0("C", cid)]]$mean
  })
  names(cluster_means) <- as.character(unique_clusters)
  
  return(list(
    expt_id = expt_id,
    beads_filt = beads_filt,
    bins_plot = bins_plot,
    counts_log2_plot = counts_log2_plot,
    unique_clusters = as.character(unique_clusters),
    cluster_means = cluster_means,
    result_dir = result_dir
  ))
}

# ==================== 3. Adjacent Slice Tumor CNV Similarity ==========================================

calc_cnv_similarity_stats <- function(profile_a, profile_b, gain_loss_cutoff = 0.25, n_boot = 1000, seed = 123) {
  df <- data.frame(
    cnv_a = as.numeric(profile_a),
    cnv_b = as.numeric(profile_b)
  )
  df <- df[is.finite(df$cnv_a) & is.finite(df$cnv_b), ]
  
  if (nrow(df) < 20) stop("Too few valid bins to stably evaluate similarity.")
  
  spearman_rho <- suppressWarnings(cor(df$cnv_a, df$cnv_b, method = "spearman"))
  pearson_r <- suppressWarnings(cor(df$cnv_a, df$cnv_b, method = "pearson"))
  
  mu_a <- mean(df$cnv_a)
  mu_b <- mean(df$cnv_b)
  var_a <- var(df$cnv_a)
  var_b <- var(df$cnv_b)
  cov_ab <- cov(df$cnv_a, df$cnv_b)
  ccc <- (2 * cov_ab) / (var_a + var_b + (mu_a - mu_b)^2)
  
  state_a <- ifelse(df$cnv_a > gain_loss_cutoff, 1L, ifelse(df$cnv_a < -gain_loss_cutoff, -1L, 0L))
  state_b <- ifelse(df$cnv_b > gain_loss_cutoff, 1L, ifelse(df$cnv_b < -gain_loss_cutoff, -1L, 0L))
  union_event <- sum((state_a != 0L) | (state_b != 0L))
  intersect_event <- sum((state_a == state_b) & (state_a != 0L))
  jaccard_gain_loss <- ifelse(union_event == 0, NA_real_, intersect_event / union_event)
  
  # Bland-Altman statistics
  df$ba_mean <- (df$cnv_a + df$cnv_b) / 2
  df$ba_diff <- df$cnv_b - df$cnv_a
  ba_bias <- mean(df$ba_diff)
  ba_sd <- sd(df$ba_diff)
  ba_loa_upper <- ba_bias + 1.96 * ba_sd
  ba_loa_lower <- ba_bias - 1.96 * ba_sd
  
  set.seed(seed)
  n <- nrow(df)
  boot_rho <- replicate(n_boot, {
    idx <- sample.int(n, size = n, replace = TRUE)
    suppressWarnings(cor(df$cnv_a[idx], df$cnv_b[idx], method = "spearman"))
  })
  spearman_ci <- unname(stats::quantile(boot_rho, probs = c(0.025, 0.975), na.rm = TRUE))
  
  list(
    metrics = data.frame(
      spearman_rho = spearman_rho,
      spearman_ci_low = spearman_ci[1],
      spearman_ci_high = spearman_ci[2],
      pearson_r = pearson_r,
      ccc = ccc,
      jaccard_gain_loss = jaccard_gain_loss,
      ba_bias = ba_bias,
      ba_loa_lower = ba_loa_lower,
      ba_loa_upper = ba_loa_upper,
      bins_used = n
    ),
    scatter_df = df
  )
}

compare_adjacent_tumor_clusters <- function(
    slice_a_result,
    slice_b_result,
    tumor_cluster_a,
    tumor_cluster_b,
    label_a = "Slice_A",
    label_b = "Slice_B",
    gain_loss_cutoff = 0.25,
    n_boot = 1000,
    seed = 123,
    out_dir = NULL,
    file_prefix = "Adjacent_Slice_Tumor_CNV_Similarity"
) {
  if (is.null(slice_a_result$cluster_means) || is.null(slice_b_result$cluster_means)) {
    stop("Input object is missing 'cluster_means'. Set return_details = TRUE in run_slide_dna_mouse().")
  }
  
  tumor_cluster_a <- as.character(tumor_cluster_a)
  tumor_cluster_b <- as.character(tumor_cluster_b)
  
  if (!(tumor_cluster_a %in% names(slice_a_result$cluster_means))) {
    stop(sprintf("Cluster '%s' does not exist in Slice A. Available clusters: %s",
                 tumor_cluster_a, paste(names(slice_a_result$cluster_means), collapse = ", ")))
  }
  if (!(tumor_cluster_b %in% names(slice_b_result$cluster_means))) {
    stop(sprintf("Cluster '%s' does not exist in Slice B. Available clusters: %s",
                 tumor_cluster_b, paste(names(slice_b_result$cluster_means), collapse = ", ")))
  }
  
  profile_a <- slice_a_result$cluster_means[[tumor_cluster_a]]
  profile_b <- slice_b_result$cluster_means[[tumor_cluster_b]]
  
  if (length(profile_a) != length(profile_b)) {
    stop("Genomic bin counts between the two slices do not match. Cannot compare directly.")
  }
  
  sim <- calc_cnv_similarity_stats(
    profile_a = profile_a,
    profile_b = profile_b,
    gain_loss_cutoff = gain_loss_cutoff,
    n_boot = n_boot,
    seed = seed
  )
  
  m <- sim$metrics[1, ]
  lim_hex <- range(c(sim$scatter_df$cnv_a, sim$scatter_df$cnv_b), na.rm = TRUE)
  
  p_hexbin <- ggplot(sim$scatter_df, aes(x = cnv_a, y = cnv_b)) +
    geom_point(color = "grey72", alpha = 0.35, size = 0.6) +
    stat_density_2d(
      aes(fill = after_stat(level)),
      geom = "polygon",
      contour = TRUE,
      n = 220,
      bins = 16,
      color = NA,
      alpha = 0.92
    ) +
    scale_fill_gradientn(
      colors = c("#D9D9D9", "#AEB6C2", "#3B4CC0", "#1F3A93", "#D73027", "#FEE08B"),
      name = "Density"
    ) +
    geom_smooth(method = "lm", se = FALSE, color = TRACK_CLUSTER_PALETTE[4], linewidth = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    coord_equal(xlim = lim_hex, ylim = lim_hex, expand = TRUE) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    labs(
      title = sprintf("Hexbin Correlation: %s C%s vs %s C%s", label_a, tumor_cluster_a, label_b, tumor_cluster_b),
      subtitle = sprintf("Spearman=%.3f (95%%CI %.3f~%.3f), Pearson=%.3f, CCC=%.3f, Jaccard=%.3f",
                         m$spearman_rho, m$spearman_ci_low, m$spearman_ci_high,
                         m$pearson_r, m$ccc, m$jaccard_gain_loss),
      x = paste0(label_a, " CNV (log2 ratio)"),
      y = paste0(label_b, " CNV (log2 ratio)")
    )
  
  p_bland_altman <- ggplot(sim$scatter_df, aes(x = ba_mean, y = ba_diff)) +
    geom_point(color = TRACK_CLUSTER_PALETTE[5], alpha = 0.22, size = 0.65) +
    geom_hline(yintercept = m$ba_bias, color = TRACK_CLUSTER_PALETTE[4], linewidth = 0.8) +
    geom_hline(yintercept = m$ba_loa_upper, color = "grey45", linetype = "dashed", linewidth = 0.6) +
    geom_hline(yintercept = m$ba_loa_lower, color = "grey45", linetype = "dashed", linewidth = 0.6) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    labs(
      title = "Bland-Altman Plot",
      subtitle = sprintf("Bias=%.3f, LoA=[%.3f, %.3f]", m$ba_bias, m$ba_loa_lower, m$ba_loa_upper),
      x = "Mean CNV of two slices",
      y = paste0(label_b, " - ", label_a, " (log2 ratio)")
    )
  
  heat_df <- rbind(
    data.frame(sample = paste0(label_a, "_C", tumor_cluster_a), bin = seq_along(profile_a), log2 = profile_a),
    data.frame(sample = paste0(label_b, "_C", tumor_cluster_b), bin = seq_along(profile_b), log2 = profile_b)
  )
  heat_df$sample <- factor(heat_df$sample, levels = c(paste0(label_a, "_C", tumor_cluster_a), paste0(label_b, "_C", tumor_cluster_b)))
  
  p_heat <- ggplot(heat_df, aes(x = bin, y = sample, fill = log2)) +
    geom_tile() +
    scale_fill_gradient2(
      low = CNV_REDBLUE_PALETTE["low"],
      mid = CNV_REDBLUE_PALETTE["mid"],
      high = CNV_REDBLUE_PALETTE["high"],
      midpoint = 0,
      name = "CNV"
    ) +
    labs(x = "Genomic Bins (Autosomes)", y = NULL) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  p_top <- p_hexbin | p_bland_altman
  p_final <- p_top / p_heat + plot_layout(heights = c(2.0, 1.0))
  
  if (is.null(out_dir)) {
    out_dir <- slice_a_result$result_dir
  }
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_path <- file.path(out_dir, paste0(file_prefix, ".pdf"))
    ggsave(out_path, p_final, width = 8.5, height = 7.2)
  } else {
    out_path <- NA_character_
  }
  
  return(list(
    metrics = sim$metrics,
    plot_hexbin = p_hexbin,
    plot_bland_altman = p_bland_altman,
    plot = p_final,
    output_pdf = out_path
  ))
}

# ==================== 4. Local Run Template ====================================================
# Usage Instructions:
# 1) Modify the base directory and expt_id below.
# 2) Set run_local to TRUE.
# 3) Source this script to automatically run Tumor C1 similarity for slices 13_1 vs 13_2.

run_local <- FALSE

if (isTRUE(run_local)) {
  repo_root <- Sys.getenv("AXIS_REPO_ROOT", unset = normalizePath("../..", winslash = "/"))
  base_dir <- repo_root
  normal_bulk_path <- file.path(base_dir, "data/bulk/mouse_blood_1Mb_total_reads.txt")
  gc_path <- file.path(base_dir, "reference/gc_content/GRCm39_1Mb_gc.txt")
  map_path <- file.path(base_dir, "reference/mappability/GRCm39_1Mb_map.txt")

  slice_12_2 <- run_slide_dna_mouse(
    expt_id = "12_2",
    base_dir = base_dir,
    normal_bulk_path = normal_bulk_path,
    gc_path = gc_path,
    map_path = map_path,
    return_details = TRUE
  )
  
  # Slice 12_4: Inverse mapping if cluster 2 is the actual Tumor
  slice_12_4 <- run_slide_dna_mouse(
    expt_id = "12_4",
    base_dir = base_dir,
    normal_bulk_path = normal_bulk_path,
    gc_path = gc_path,
    map_path = map_path,
    k_clusters = 2,
    return_details = TRUE,
    cluster_relabel = c("1" = "Tumor", "2" = "NonTumor")
  )
  
  sim_res <- compare_adjacent_tumor_clusters(
    slice_a_result = slice_12_2,
    slice_b_result = slice_12_4,
    tumor_cluster_a = "Tumor",
    tumor_cluster_b = "Tumor",
    label_a = "12_2",
    label_b = "12_4",
    file_prefix = "12_2_vs_12_4_Tumor"
  )
  
  print(sim_res$metrics)
  cat("Similarity PDF:", sim_res$output_pdf, "\n")
}

# AXIS-DNA-seq: MEDICC2 export and phylogeny for merged multi-sample CNV clusters.

prepare_bins_plot <- function(bins) {
  bins <- as.data.frame(bins)
  bins$chr_clean <- gsub("^chr", "", bins$chr)
  valid_idx <- which(!bins$chr_clean %in% c("X", "Y"))
  bins_plot <- bins[valid_idx, , drop = FALSE]
  bins_plot$chr_factor <- factor(bins_plot$chr_clean, levels = as.character(1:22))
  bins_plot[order(bins_plot$chr_factor, bins_plot$bin_start), , drop = FALSE]
}


build_sample_cluster_log2_profiles <- function(
    beads_filt,
    counts_norm_raw,
    target_bins,
    bulk_norm_vec = NULL,
    min_spots = 5
) {
  combos <- unique(beads_filt[, c("sample_id", "cluster")])
  profiles <- list()
  meta <- list()

  raw_mat <- as.matrix(counts_norm_raw[, target_bins, drop = FALSE])

  ref_vec <- if (!is.null(bulk_norm_vec)) {
    bulk_norm_vec[target_bins]
  } else {
    apply(raw_mat, 2, median, na.rm = TRUE)
  }
  ref_vec[!is.finite(ref_vec) | ref_vec <= 0] <- 1e-6

  for (i in seq_len(nrow(combos))) {
    sid <- as.character(combos$sample_id[i])
    cl <- as.character(combos$cluster[i])
    cells <- which(beads_filt$sample_id == sid & as.character(beads_filt$cluster) == cl)
    if (length(cells) < min_spots) {
      cat(sprintf("   -> Skip %s_C%s: only %d spots (< %d)\n", sid, cl, length(cells), min_spots))
      next
    }
    prof <- if (length(cells) > 1) colMeans(raw_mat[cells, , drop = FALSE]) else raw_mat[cells, ]
    ratio <- prof / ref_vec
    log2_val <- log2(ratio + 1e-6)
    log2_val[log2_val > 3] <- 3
    log2_val[log2_val < -3] <- -3
    clone_id <- paste0(sid, "_C", cl)
    profiles[[clone_id]] <- log2_val
    meta[[clone_id]] <- data.frame(
      clone_id = clone_id,
      sample_id = sid,
      cluster = cl,
      n_spots = length(cells),
      stringsAsFactors = FALSE
    )
  }

  list(
    profiles = profiles,
    meta = if (length(meta) > 0) do.call(rbind, meta) else data.frame()
  )
}


export_merged_for_medicc2 <- function(
    profiles,
    bins_plot,
    out_dir,
    project_name
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("\n=== Export MEDICC2 input: %s ===\n", project_name))

  if (length(profiles) < 2) {
    stop("Need at least 2 sample-cluster profiles for MEDICC2 (got ", length(profiles), ")")
  }

  medicc_df_list <- list()
  for (clone_id in names(profiles)) {
    log2_val <- as.numeric(profiles[[clone_id]])
    abs_cn <- round(2 * (2 ^ log2_val))
    abs_cn[abs_cn > 8] <- 8
    abs_cn[abs_cn < 0] <- 0
    medicc_df_list[[clone_id]] <- data.frame(
      sample_id = clone_id,
      chrom = bins_plot$chr_clean,
      start = bins_plot$bin_start,
      end = bins_plot$bin_end,
      cn = abs_cn,
      stringsAsFactors = FALSE
    )
  }

  final_medicc_df <- do.call(rbind, medicc_df_list)
  normal_df <- data.frame(
    sample_id = "Diploid_Root",
    chrom = bins_plot$chr_clean,
    start = bins_plot$bin_start,
    end = bins_plot$bin_end,
    cn = 2,
    stringsAsFactors = FALSE
  )
  final_medicc_df <- rbind(normal_df, final_medicc_df)

  out_file <- file.path(out_dir, paste0(project_name, "_medicc2_input.tsv"))
  write.table(final_medicc_df, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("   -> MEDICC2 input: %s (%d clones)\n", out_file, length(profiles)))
  list(input_file = out_file, table = final_medicc_df)
}


run_medicc2_cli <- function(input_file, out_dir, medicc2_bin = "medicc2") {
  medicc2_path <- Sys.which(medicc2_bin)
  if (!nzchar(medicc2_path)) {
    cat("   -> medicc2 not found in PATH; export only.\n")
    cat(sprintf("      Run manually: %s %s %s\n", medicc2_bin, input_file, out_dir))
    return(NULL)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cmd <- paste(
    shQuote(medicc2_path), shQuote(input_file), shQuote(out_dir),
    "-j 8 --plot both --total-copy-numbers --normal-name Diploid_Root"
  )
  cat("   -> Running:", cmd, "\n")
  status <- system(cmd)
  if (status != 0) {
    warning("medicc2 exited with status ", status)
    return(NULL)
  }
  tree_files <- list.files(out_dir, pattern = "final_tree\\.new$", recursive = TRUE, full.names = TRUE)
  if (length(tree_files) == 0) {
    tree_files <- list.files(out_dir, pattern = "\\.new$", recursive = TRUE, full.names = TRUE)
  }
  if (length(tree_files) == 0) return(NULL)
  tree_files[[1]]
}


plot_medicc2_tree <- function(tree_file, clone_meta, fig_dir, project_name) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    cat("   -> Package 'ape' not installed; skip tree plot.\n")
    return(invisible(NULL))
  }
  tree <- ape::read.tree(tree_file)
  pdf(file.path(fig_dir, paste0(project_name, "_MEDICC2_tree.pdf")), width = 14, height = 10)
  plot(tree, type = "phylogram", direction = "rightwards", edge.width = 1.5, cex = 0.7,
       main = paste0(project_name, " — MEDICC2 CNV Evolution"))
  axisPhylo(backward = FALSE)
  dev.off()

  if (requireNamespace("ggtree", quietly = TRUE) && nrow(clone_meta) > 0) {
    metadata <- data.frame(Tip = tree$tip.label, stringsAsFactors = FALSE)
    metadata$sample_id <- vapply(metadata$Tip, function(tip) {
      if (tip == "Diploid_Root") return("Diploid_Root")
      m <- clone_meta$sample_id[match(tip, clone_meta$clone_id)]
      if (is.na(m)) sub("_C[0-9]+$", "", tip) else m
    }, character(1))

    p <- ggtree::ggtree(tree, size = 0.8, color = "#4A4A4A") %<+% metadata +
      ggtree::geom_tippoint(ggplot2::aes(color = sample_id), size = 3) +
      ggtree::geom_tiplab(ggplot2::aes(label = Tip), size = 2.5, offset = 0.02) +
      ggplot2::theme_tree2() +
      ggplot2::labs(title = paste0(project_name, " — MEDICC2 (colored by sample)"))
    ggplot2::ggsave(
      file.path(fig_dir, paste0(project_name, "_MEDICC2_tree_ggtree.pdf")),
      p, width = 14, height = 10
    )
  }
  invisible(tree)
}


run_merged_medicc2_pipeline <- function(
    beads_filt,
    counts_norm_raw,
    bins,
    bulk_norm_vec = NULL,
    out_dir,
    merged_id,
    min_spots = 5,
    run_medicc2 = TRUE,
    medicc2_bin = "medicc2"
) {
  medicc2_dir <- file.path(out_dir, "MEDICC2")
  input_dir <- file.path(medicc2_dir, "input")
  output_dir <- file.path(medicc2_dir, "output")
  fig_dir <- file.path(out_dir, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  bins_plot <- prepare_bins_plot(bins)
  target_bins <- match(
    paste(bins_plot$chr, bins_plot$bin_start, bins_plot$bin_end),
    paste(bins$chr, bins$bin_start, bins$bin_end)
  )
  target_bins <- target_bins[!is.na(target_bins)]

  prof_res <- build_sample_cluster_log2_profiles(
    beads_filt = beads_filt,
    counts_norm_raw = counts_norm_raw,
    target_bins = target_bins,
    bulk_norm_vec = bulk_norm_vec,
    min_spots = min_spots
  )

  if (nrow(prof_res$meta) < 2) {
    warning("Too few sample-cluster profiles for MEDICC2; skipping.")
    return(NULL)
  }

  fwrite(prof_res$meta, file.path(out_dir, paste0(merged_id, "_sample_cluster_profiles.csv")))

  export_res <- export_merged_for_medicc2(
    profiles = prof_res$profiles,
    bins_plot = bins_plot,
    out_dir = input_dir,
    project_name = merged_id
  )

  tree_file <- NULL
  if (isTRUE(run_medicc2)) {
    tree_file <- run_medicc2_cli(export_res$input_file, output_dir, medicc2_bin = medicc2_bin)
  }

  if (!is.null(tree_file) && file.exists(tree_file)) {
    plot_medicc2_tree(tree_file, prof_res$meta, fig_dir, merged_id)
    cat(sprintf("   -> Tree saved from: %s\n", tree_file))
  }

  list(
    profiles = prof_res$profiles,
    meta = prof_res$meta,
    medicc2_input = export_res$input_file,
    tree_file = tree_file
  )
}

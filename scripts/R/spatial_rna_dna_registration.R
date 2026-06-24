# AXIS-DNA-seq: Spatial RNA-DNA registration and QC.
#
# Section 0: barcode grid assignment from tissue image coordinates
# Section 1-3: SpotClean / Seurat RNA preprocessing (optional)
# Section 4: Procrustes alignment between adjacent DNA and RNA slices
rm(list = ls())
gc()  # free memory

# Load packages for preprocessing
library(dplyr)
library(tidyr)

repo_root <- Sys.getenv("AXIS_REPO_ROOT", unset = normalizePath("../..", winslash = "/"))
setwd(repo_root)
source(file.path(repo_root, "scripts/R/spatial_rna_helpers.R"))

sample_id <- "test10" #test10, 11_1, 11_2, HCC, LSIL; 12_2, 12_4, 13_1
species <- "human"  #human; mouse
resolution <- "80"  # 50 or 80

# output dir
output_dir <- paste0("./results/Spatial_RNA/",sample_id)
ensure_dir(output_dir)


#### Fetch slide information based on per-sample barcode location ####
# Raw tissue outline: comma-separated rowxcol tokens (e.g. position_1_3.txt)
filename <- paste0("./data/Spatial_RNA/", sample_id, "/position_", sample_id, ".txt")
location <- parse_raw_tissue_position(filename)

barcode_index <- resolve_spatial_rna_barcode_index(
  sample_id,
  repo_root = repo_root,
  resolution = as.integer(resolution)
)

# tissue barcode index
slide_info <- fetch_side_infor(barcode_index, location, 1)

# save result
filename <- paste0(output_dir,"/position_",sample_id,".txt")
write.table(slide_info, file = filename, sep = "\t", row.names = F, quote = F)



#### 1: inferCNV analysis ####
rm(list = ls())
gc()

#### Step 1: R environment and data loaders ####
library(SpotClean)
library(Seurat)
library(SummarizedExperiment)
library(infercnv)
library(tidyverse)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(Matrix)
library(dplyr)
library(tidyr)
library(OpenImageR)
library(grid)

sample_id <- "13_1" #test10, 11_1, 11_2, HCC, LSIL; 12_2, 12_4, 13_1
species <- "mouse"  #human; mouse
k <- 2 # equal to the subclone number of AXIS-seq 

repo_root <- Sys.getenv("AXIS_REPO_ROOT", unset = normalizePath("../..", winslash = "/"))
setwd(repo_root)
source(file.path(repo_root, "scripts/R/spatial_rna_helpers.R"))

resolution <- DEFAULT_SPATIAL_RESOLUTION  # 50 or 80

# output dir
output_dir <- paste0("./results/Spatial_RNA/",sample_id)
output_inferCNV <- paste0(output_dir,"/1.inferCNV/")
ensure_dir(output_inferCNV)


#### Step 2: Prepare data ####
# inferCNV requires raw counts
# read cellRanger raw matrix (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
data_dir <- resolve_cellranger_raw_matrix_dir(sample_id, repo_root = repo_root)
cellRanger_result <- read10xRaw(data_dir)
# Create a new slide object
slide_info <- read.table(file =  paste0(output_dir,"/position_",sample_id,".txt"),sep = "\t",header = T)
slide_obj <- createSlide(cellRanger_result, slide_info)
 
# data in tissue object
keep_cols <- metadata(slide_obj)$slide$tissue == 1
slide_obj_subset <- slide_obj[, keep_cols]
metadata(slide_obj_subset)$slide <- metadata(slide_obj_subset)$slide[keep_cols, ]

# Label all spots as one group (e.g. "Mixed") when identities are unknown
annotations <- data.frame(cell_id = colnames(slide_obj_subset), group = "Mixed")
write.table(annotations, file = paste0(output_inferCNV,sample_id,"_infercnv_annotations.txt"), sep = "\t", quote = F, row.names = F, col.names = F)

# Create inferCNV Object
counts_matrix <- assay(slide_obj_subset, "raw")
dim(counts_matrix)

if(species == "human"){
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts_matrix,
    annotations_file = paste0(output_inferCNV,sample_id,"_infercnv_annotations.txt"),
    gene_order_file = "./reference/gene_order_file/gene_ordering_file_from_ref_human.txt", # species-matched gene order file
    ref_group_names = NULL,  # NULL: no reference group; set to e.g. "CON" to specify a normal reference
    min_max_counts_per_cell = c(1, +Inf)  # default min_max_counts_per_cell = 100 filters low-count spots
  )
}else if(species == "mouse"){
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts_matrix,
    annotations_file = paste0(output_inferCNV,sample_id,"_infercnv_annotations.txt"),
    gene_order_file = "./reference/gene_order_file/gene_ordering_file_from_ref_mouse.txt", # species-matched gene order file
    ref_group_names = NULL,  # NULL: no reference group; set to e.g. "CON" to specify a normal reference
    min_max_counts_per_cell = c(1, +Inf)  # default min_max_counts_per_cell = 100 filters low-count spots
  )
}


#### Step 3: Run analysis ####
infercnv_run <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1, 
  out_dir = paste0(output_inferCNV,"infercnv_output_basic"), 
  cluster_by_groups = FALSE, # single "Mixed" group: cluster globally (FALSE)
  analysis_mode = "subclusters", # find subclones within the group
  min_cells_per_gene = 1,     # filter genes expressed in too few spots
  denoise = TRUE,
  HMM = TRUE,
  write_expr_matrix = T,
  num_threads = 6  # number of threads
)
save(infercnv_run,file = paste0(output_inferCNV,"_infercnv_run_result.Rdata"))


#### Step 4: Visualize results ####
# map inferCNV subclone labels back to spatial coordinates
grouping_file <- paste0(output_inferCNV,"infercnv_output_basic/infercnv.observation_groupings.txt")  # verify exact filename in output directory
subclones <- read.table(grouping_file, header = TRUE, sep = "", check.names = TRUE)
# Visualization
raw_inferCNV_p <- inferCNV_spatila_cluster(slide_obj_subset,subclones,resolution)
ggsave(raw_inferCNV_p, filename = paste0(output_inferCNV,"/raw_inferCNV_observation_subclones.pdf"), width = 12, height = 12)


#### Step 5: Re-cluster CNV matrix with manual k ####
load(file = paste0(output_inferCNV,"_infercnv_run_result.Rdata"))
# extract expression matrix
m_cnv_matrix <- infercnv_run@expr.data
# Calculate the distance between spots, simulated the clustering process within inferCNV 
dist_mat <- dist(t(m_cnv_matrix)) # transpose: spot-to-spot distance
hc <- hclust(dist_mat, method = "ward.D2") # hierarchical clustering

# Visualization
# source(file.path(repo_root, "scripts/R/spatial_rna_helpers.R"))
K_inferCNV_p <- inferCNV_spatila_K(slide_obj_subset,hc,resolution,k)
ggsave(K_inferCNV_p, filename = paste0(output_inferCNV,"/reClustered_inferCNV_subclones.pdf"), width = 6.2, height = 6.2)


#### 2: decontamination ####
rm(list = ls())
gc()

#### Step 1: R environment ####
library(SpotClean)
library(Seurat)
library(SummarizedExperiment)
library(tidyverse)
library(data.table)
library(ggplot2)
library(Matrix)
library(dplyr)
library(tidyr)
library(OpenImageR)
library(grid)

sample_id <- "test10" #test10, 11_1, 11_2; 12_2
species <- "human"  #human; mouse

repo_root <- Sys.getenv("AXIS_REPO_ROOT", unset = normalizePath("../..", winslash = "/"))
setwd(repo_root)
source(file.path(repo_root, "scripts/R/spatial_rna_helpers.R"))

resolution <- DEFAULT_SPATIAL_RESOLUTION  # 50 or 80

# output dir
output_dir <- paste0("./results/Spatial_RNA/",sample_id)
output_decontaminate <- paste0(output_dir,"/2.Decontaminate/")
ensure_dir(output_decontaminate)


#### Step 2: Decontaminate data ####

# Decontamination (e.g. SpotClean) requires raw unprocessed data
# SpotClean models bleeding rate by comparing tissue vs blank spot expression.

# read cellRanger raw matrix (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
data_dir <- resolve_cellranger_raw_matrix_dir(sample_id, repo_root = repo_root)
cellRanger_result <- read10xRaw(data_dir)
# read slide position information
slide_info <- read.table(file =  paste0(output_dir,"/position_",sample_id,".txt"),sep = "\t",header = T)
# Create a new slide object
slide_obj <- createSlide(cellRanger_result, slide_info)

# Decontaminated
decont_obj <- spotclean(slide_obj, tol=10, candidate_radius=20)

# Visualization of contamination rate
decont_rate_plot <- contamination_heatmap(decont_obj, resolution)
ggsave(decont_rate_plot, filename = paste0(output_decontaminate,"/Contamination_heatmap.pdf"), width = 8.6, height = 8.6)

# extract decontaminated matrix and map spatial indices (row x col)
decont_matrix <- data.frame(assay(decont_obj),check.names = F)
slide_info$index <- paste0(slide_info$row, "x", slide_info$col)
colnames(decont_matrix) <- slide_info[match(colnames(decont_matrix), slide_info$barcode), "index"]
# transpose and save backup for faster reload
data_filtered <- t(decont_matrix)
filename <- paste0(output_decontaminate,"/Decontaminated_filtered_matrix.tsv")
write.table(cbind(Spot_ID = rownames(data_filtered), data_filtered), file = filename, sep = "\t", col.names = T, row.names = F, quote = F)

# re_read
my_data <- read.table(paste0(output_decontaminate,"/Decontaminated_filtered_matrix.tsv"), sep = "\t", header = T, stringsAsFactors = FALSE, check.names = F)
names(my_data)[1] = "X"
count <- rowSums(my_data[, 2:ncol(my_data)])
data_filtered_binary <- my_data[, 2:ncol(my_data)] %>% mutate_all(as.logical)
gene_count <- rowSums(data_filtered_binary)

# UMI Count 
df <- data.frame(number = 1, c = count)
region <- max(count)
UMI_plot <- umi_plot(df, region)
filename <- paste0(output_decontaminate,"/Decontaminated_UMI_counts.pdf")
ggsave(UMI_plot, filename=filename, width = 8.6, height = 8.6)

# Gene Count
df <- data.frame(number = 1, c = gene_count)
region <- max(gene_count)
Gene_plot <- gene_plot(df, region)
filename <- paste0(output_decontaminate,"/Decontaminated_Gene_counts.pdf")
ggsave(Gene_plot, filename=filename, width = 8.6, height = 8.6)

# UMI heatmap, adjust the limits for scale_color_gradientn, select the limit to be close to the maximum number
test <- my_data %>% separate(X, c("A", "B"),  sep = "x")
color_manual <- c("#252A62","#692F7C", "#B43970", "#D96558", "#EFA143", "#F6C63C")
UMI_heatmap_plot <- umi_heatmap(test, color_manual, max(count), resolution)
filename <- paste0(output_decontaminate,"/Decontaminated_UMI_count_heatmap.pdf")
ggsave(UMI_heatmap_plot, filename=filename, width = 8.6, height = 8.6)

# Gene heatmap, adjust the limits for scale_color_gradientn, select the limit to be close to the maximum number
Gene_heatmap_plot <- gene_heatmap(test, color_manual, max(gene_count), resolution)
filename <- paste0(output_decontaminate,"/Decontaminated_Gene_count_heatmap.pdf")
ggsave(Gene_heatmap_plot, filename=filename, width = 8.6, height = 8.6)



#### Step 3: basic QC ####
# Create Decontaminated Seurate object 
matrix_data <- Matrix(as.matrix(decont_matrix), sparse = TRUE)
seurat_obj <- CreateSeuratObject(matrix_data, min.cells = 10)
filename <- paste0(output_decontaminate,"/Decontaminated_seurat_obj.Rdata")
save(seurat_obj,file = filename)

# mitochondrial gene fraction
# human: "^MT-"; mouse: "^mt-"
if(species == "human"){
  seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")
}else if (species == "mouse"){
  seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^mt-", col.name = "percent.mt")
}

# filter spots with fewer than 200 detected genes
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200)
filename <- paste0(output_decontaminate,"/Decontaminated_QCed_seurat_obj.Rdata")
save(seurat_obj,file = filename)




#### 3: Basic analysis ####
rm(list = ls())
gc()

#### Step 1: R environment ####
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(rhdf5)
library(Matrix)
library(sctransform)
library(plyr)
library(gridExtra)
library(magrittr)
library(tidyr)
library(raster)
library(OpenImageR)
library(ggpubr)
library(grid)
library(wesanderson)
library(glmGamPoi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

sample_id <- "test10" #test10, 11_1, 11_2; 12_2
species <- "human"  #human; mouse

repo_root <- Sys.getenv("AXIS_REPO_ROOT", unset = normalizePath("../..", winslash = "/"))
setwd(repo_root)
source(file.path(repo_root, "scripts/R/spatial_rna_helpers.R"))

resolution <- DEFAULT_SPATIAL_RESOLUTION  # 50 or 80

# output dir
output_dir <- paste0("./results/Spatial_RNA/",sample_id)
output_Basic <- paste0(output_dir,"/3.Basic_analysis/")
ensure_dir(output_Basic)


#### Step 2: Clustering and marker analysis ####
# load data
load(file = paste0(output_dir,"/2.Decontaminate/Decontaminated_QCed_seurat_obj.Rdata"))

qc_metrics <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
filename <- paste0(output_Basic,"/QC_metrics.pdf")
ggsave(qc_metrics, filename=filename, width = 8.6, height = 8.6)

plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Feature_cor <- plot1 + plot2
filename <- paste0(output_Basic,"/Fature_cor.pdf")
ggsave(Feature_cor, filename=filename, width = 8.6, height = 8.6)

# normalize
seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)

# cluster
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)

# visualization
umap_plot <- DimPlot(seurat_obj, label = TRUE, cols = c( "#4A90E2",  "#7ED321",   "#FF9500","#E2B4B4", "#F6C8A8", "#E6CECF", "#D9B9A3", "#C69D9D", "#BCAEAA", "#F3E5E2")) + 
  xlab("UMAP_1") + ylab("UMAP_2")





#### 4: DNA-RNA alignment ####

## Procrustes registration of adjacent DNA (AXIS-seq) and RNA slices
## stratified sampling, coverage map, displacement field, smoothed distance, label agreement QC

## ----- install dependencies (once) -----
# install.packages("MASS")     # Procrustes registration
# install.packages("dplyr")    # data manipulation
# install.packages("ggplot2")  # visualization
# install.packages("fields")   # Tps smoothing
# install.packages("patchwork")# patchwork layouts

library(MASS)
library(dplyr)
library(ggplot2)
library(fields)
library(patchwork)
library(arrow)

grid_res <- as.integer(resolution)
if (is.na(grid_res) || grid_res <= 0) {
  grid_res <- DEFAULT_SPATIAL_RESOLUTION
}
grid_n <- grid_res * grid_res
spot_point_size <- spatial_grid_point_size(grid_res)
spatial_coord <- function(p) {
  p + coord_fixed(ratio = 1, xlim = c(0.5, grid_res + 0.5), ylim = c(grid_res + 0.5, 0.5), expand = FALSE)
}

## ===== Step 1: build spatial grid coordinates =====
# all spot (x, y) coordinates on an grid_res x grid_res array
grid_xy <- expand.grid(x = 1:grid_res, y = 1:grid_res)

# flatten cluster matrix column-wise and bind to grid
# requires top_cluster and bottom_cluster in environment:
#         grid_res x grid_res integer matrix; 0 = background, 1/2/3 = subclone
top_vec    <- as.vector(top_cluster)    # length grid_n
bottom_vec <- as.vector(bottom_cluster)

dat_top <- cbind(grid_xy, cluster = top_vec)
dat_bot <- cbind(grid_xy, cluster = bottom_vec)

## ===== Step 2: stratified subclone sampling =====
# exclude background (cluster == 0)
top_org <- dat_top %>% filter(cluster != 0)
bot_org <- dat_bot %>% filter(cluster != 0)

# stratified sampling per subclone (not global random)
# equal fraction per subclone to avoid large-clone bias
set.seed(123)

# sample equal fractions from top and bottom per subclone
bot_frac <- 0.3   # 30% of spots per subclone for registration; adjust as needed

bot_s <- bot_org %>%
  group_by(cluster) %>%
  sample_frac(bot_frac) %>%
  ungroup()

# stratified sampling on top slice, matched to bottom
top_s <- top_org %>%
  group_by(cluster) %>%
  sample_frac(bot_frac) %>%
  ungroup()

# use the smaller matched spot count
n_match <- min(nrow(top_s), nrow(bot_s))
set.seed(456)
top_s <- top_s %>% sample_n(n_match)
bot_s <- bot_s %>% sample_n(n_match)

X <- as.matrix(bot_s[, c("x", "y")])   # moving (to be transformed)
Y <- as.matrix(top_s[, c("x", "y")])   # fixed (reference frame)

cat(sprintf("Registration sample size: %d\n", n_match))

## ===== Step 3: Procrustes transform bottom -> top =====
proc  <- procrustes(Y, X, scale = TRUE)  # transform X into Y space

Rmat  <- proc$rotation * proc$scale      # 2x2 rotation + scale matrix
trans <- proc$translation                 # 1x2 translation vector

# apply transform from bottom to top coordinates
transform_xy <- function(xy) {
  sweep(as.matrix(xy) %*% Rmat, 2, trans, "+")
}

## ===== Step 4: map all grid spots to top coordinates =====
bot_all_xy_tr <- as.data.frame(transform_xy(grid_xy))
colnames(bot_all_xy_tr) <- c("x_t", "y_t")

## ===== Step 5: round coordinates and filter out-of-bounds =====
bot_all_xy_tr$xr <- round(bot_all_xy_tr$x_t)
bot_all_xy_tr$yr <- round(bot_all_xy_tr$y_t)

valid        <- with(bot_all_xy_tr, xr >= 1 & xr <= grid_res & yr >= 1 & yr <= grid_res)
bot_valid_idx <- which(valid)

# linear index into top matrix (column-major)
index_top <- with(bot_all_xy_tr[bot_valid_idx, ], (xr - 1) * grid_res + yr)

## ===== Step 6: transfer cluster labels =====
top_vec_full <- as.vector(top_cluster)
shared_label <- rep(NA_integer_, grid_n)
shared_label[bot_valid_idx] <- top_vec_full[index_top]

# build result data frame
bottom_shared <- data.frame(
  x_bottom        = grid_xy$x[bot_valid_idx],
  y_bottom        = grid_xy$y[bot_valid_idx],
  top_x           = bot_all_xy_tr$xr[bot_valid_idx],
  top_y           = bot_all_xy_tr$yr[bot_valid_idx],
  top_cluster     = shared_label[bot_valid_idx],
  bottom_cluster  = bottom_vec[bot_valid_idx],
  # displacement vectors
  dx              = bot_all_xy_tr$xr[bot_valid_idx] - grid_xy$x[bot_valid_idx],
  dy              = bot_all_xy_tr$yr[bot_valid_idx] - grid_xy$y[bot_valid_idx],
  # Euclidean residual after rounding
  dist            = sqrt(
    (bot_all_xy_tr$x_t[bot_valid_idx] - bot_all_xy_tr$xr[bot_valid_idx])^2 +
      (bot_all_xy_tr$y_t[bot_valid_idx] - bot_all_xy_tr$yr[bot_valid_idx])^2
  )
)

## ===== Step 7: quantitative QC =====

# --- 7a. mapping coverage ---
total_bot_tissue <- sum(bottom_vec != 0)
mapped_count     <- sum(bottom_vec[bot_valid_idx] != 0)
mapping_rate     <- mapped_count / total_bot_tissue * 100

cat(sprintf("Total tissue spots (bottom): %d\n", total_bot_tissue))
cat(sprintf("Successfully mapped spots: %d\n", mapped_count))
cat(sprintf("Mapping coverage:         %.1f%%\n", mapping_rate))

# --- 7b. global registration RMSE (spot units) ---
rmse <- sqrt(mean(
  (bot_all_xy_tr$x_t - bot_all_xy_tr$xr)^2 +
    (bot_all_xy_tr$y_t - bot_all_xy_tr$yr)^2,
  na.rm = TRUE
))
cat(sprintf("Global registration RMSE: %.4f spots\n", rmse))

# --- 7c. subclone label agreement (bottom vs mapped top) ---
# computed on tissue spots with valid labels on both slices
valid_tissue <- bottom_shared %>%
  filter(bottom_cluster != 0, !is.na(top_cluster), top_cluster != 0)

agreement <- mean(valid_tissue$bottom_cluster == valid_tissue$top_cluster)
cat(sprintf("Subclone label agreement: %.1f%% (%d comparable spots)\n",
            agreement * 100, nrow(valid_tissue)))

## ===== Step 8: visualization =====

subclone_colors <- c("1" = "#4C9BE8", "2" = "#73C26F", "3" = "#F0A84F",
                     "Unmapped" = "#D94F5C")

# ---------- Plot 1: AXIS-seq subclones (bottom) ---
p1 <- spatial_coord(
  dat_bot %>%
    filter(cluster != 0) %>%
    mutate(cluster = factor(cluster)) %>%
    ggplot(aes(x, y, color = cluster)) +
    geom_point(size = spot_point_size, shape = 16) +
    scale_color_manual(values = subclone_colors, name = "Subclones") +
    theme_void() +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 11)) +
    labs(title = "AXIS-seq Subclones")
)

# ---------- Plot 2: mapping coverage ---
# mapping status for all bottom tissue spots
coverage_df <- dat_bot %>%
  filter(cluster != 0) %>%
  mutate(
    lin_idx = (x - 1) * grid_res + y,
    mapped  = lin_idx %in% bot_valid_idx,
    status  = ifelse(mapped, "Mapped", "Unmapped")
  )

p2 <- spatial_coord(
  ggplot(coverage_df, aes(x, y, color = status)) +
    geom_point(size = spot_point_size, shape = 16) +
    scale_color_manual(values = c("Mapped" = "#73C26F", "Unmapped" = "#D94F5C"),
                       name = "Status") +
    annotate("text", x = 1, y = grid_res,
             label = sprintf("Mapping rate: %.1f%%", mapping_rate),
             hjust = 0, vjust = 1, size = 3.2, color = "gray30") +
    theme_void() +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 11)) +
    labs(title = "Mapping Coverage")
)

# ---------- Plot 3: mapped subclones ---
mapped_vis <- data.frame(
  x       = grid_xy$x,
  y       = grid_xy$y,
  cluster = factor(ifelse(is.na(shared_label) & bottom_vec != 0, "Unmapped",
                          ifelse(bottom_vec == 0, NA, as.character(shared_label))))
) %>% filter(!is.na(cluster))

p3 <- spatial_coord(
  ggplot(mapped_vis, aes(x, y, color = cluster)) +
    geom_point(size = spot_point_size, shape = 16) +
    scale_color_manual(
      values = c(subclone_colors, "NA" = "grey80"),
      name = "Subclones",
      na.value = "grey80"
    ) +
    theme_void() +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 11)) +
    labs(title = "Mapped Subclones")
)

# ---------- Plot 4: displacement vector field ---
# subsample arrows every arrow_step spots
arrow_step <- 4
arrow_df <- bottom_shared %>%
  filter(bottom_cluster != 0) %>%
  filter(x_bottom %% arrow_step == 0, y_bottom %% arrow_step == 0)

p4 <- ggplot(arrow_df, aes(x = x_bottom, y = y_bottom,
                           xend = x_bottom + dx, yend = y_bottom + dy)) +
  geom_segment(arrow = arrow::arrow(length = unit(0.08, "cm"), type = "closed"),
               color = "steelblue", alpha = 0.7, linewidth = 0.4) +
  coord_fixed() + theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  labs(title = "Displacement Vectors\n(bottom → top)")

# ---------- Plot 5: smoothed distance heatmap ---
# Tps smoothing of registration residuals
smooth_df <- bottom_shared %>% filter(bottom_cluster != 0)

if (nrow(smooth_df) > 10) {
  fit_tps   <- Tps(as.matrix(smooth_df[, c("x_bottom", "y_bottom")]), smooth_df$dist)
  pred_grid <- expand.grid(x_bottom = 1:grid_res, y_bottom = 1:grid_res)
  pred_grid$smooth_dist <- predict(fit_tps, as.matrix(pred_grid))
  
  p5 <- ggplot(pred_grid, aes(x_bottom, y_bottom, fill = smooth_dist)) +
    geom_tile() +
    scale_fill_viridis_c(option = "plasma", name = "Distance\n(smoothed)") +
    coord_fixed() + theme_void() +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 11)) +
    labs(title = "QC: Mapping Distance\n(Smoothed)")
} else {
  # fallback to scatter plot when too few points
  p5 <- ggplot(smooth_df, aes(x_bottom, y_bottom, color = dist)) +
    geom_point(size = spot_point_size, shape = 16) +
    scale_color_viridis_c(option = "plasma", name = "Distance") +
    coord_fixed() + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 11)) +
    labs(title = "QC: Mapping Distance")
}

# ---------- Plot 6: subclone label confusion heatmap ---
confusion_df <- valid_tissue %>%
  count(bottom_cluster, top_cluster) %>%
  group_by(bottom_cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    bottom_cluster = factor(bottom_cluster),
    top_cluster    = factor(top_cluster)
  )

p6 <- ggplot(confusion_df,
             aes(x = top_cluster, y = bottom_cluster, fill = prop)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.0f%%", prop * 100)), size = 3.5) +
  scale_fill_gradient(low = "#f7f7f7", high = "#2171b5",
                      name = "Proportion", limits = c(0, 1)) +
  labs(
    title    = "Subclone Label Agreement",
    subtitle = sprintf("Overall: %.1f%%", agreement * 100),
    x        = "Top (RNA) Cluster",
    y        = "Bottom (DNA) Subclone"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 11),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"),
    panel.grid    = element_blank()
  )

## ===== Step 9: assemble output =====
# 3-column layout: row1 p1-p3, row2 p4-p6
final_plot <- (p1 | p2 | p3) / (p4 | p5 | p6) +
  plot_annotation(
    title   = "DNA–RNA Spatial Spot Mapping",
    subtitle = sprintf(
      "RMSE = %.3f spots  |  Mapping rate = %.1f%%  |  Subclone agreement = %.1f%%",
      rmse, mapping_rate, agreement * 100
    ),
    theme = theme(
      plot.title    = element_text(hjust = 0.5, size = 13, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9,  color = "gray40")
    )
  )

# save figure
ggsave("DNA_RNA_Mapping_QC.pdf", final_plot, width = 15, height = 10)

## ===== main result object ===
# bottom_shared columns:
#   x_bottom / y_bottom    — original bottom coordinates
#   top_x / top_y          — rounded top coordinates
#   bottom_cluster         — original bottom subclone
#   top_cluster            — mapped top cluster
#   dx / dy                — displacement (top - bottom)
#   dist                   — registration residual (spot units)






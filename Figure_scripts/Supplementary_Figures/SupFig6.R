############################################################
## DIABIMMUNE Karelia Bray-Curtis trajectory analysis
##
## 1) Load metadata and MetaPhlAn species-level profiles
## 2) Filter species by prevalence and relative abundance
## 3) Restrict samples to the first 1000 days of life
## 4) Select subjects with sufficient longitudinal samples
## 5) Calculate Bray-Curtis distance to each subject's own baseline
## 6) Cluster longitudinal trajectories using kmlShape
## 7) Plot Bray-Curtis distance trajectories by cluster
############################################################

###########################
## 0. Load packages
###########################
library(dplyr)
library(tidyr)
library(vegan)
library(traj)
library(ggplot2)

###########################
## 1. Set project root and load input data 
###########################
set.seed(1)
setwd("~/ddong/TEDDY_project/revision/other_data/three_cohort/")

load("./DIABIMMUNE_Karelia_metadata.RData")

metaphlan <- fread("./diabimmune_karelia_metaphlan_table.txt")

metaphlan <- metaphlan %>%
  separate(
    ID,
    c("kingdom", "phylum", "class", "order", "family", "genus", "species", "strain"),
    sep = "\\|",
    remove = FALSE
  )

###########################
## 2. Extract species-level profiles
###########################
Species <- metaphlan %>%
  filter(is.na(strain) & !is.na(species))

Species <- Species[, -c(1:7, 9)]

###########################
## 3. Define Bray-Curtis trajectory clustering pipeline
###########################
run_bray_traj_pipeline <- function(
    Species,
    metadata,
    min_samples    = 4,
    min_RA         = 0.01,
    min_RA_percent = 10,
    nclusters      = 3
) {
  
  ###########################
  ## 3.1 Species-level QC
  ###########################
  sp_df <- as.data.frame(Species)
  
  if (!"species" %in% colnames(sp_df)) {
    stop("The first column of Species must be named 'species'.")
  }
  
  rownames(sp_df) <- sp_df$species
  sp_df$species <- NULL
  
  number_pass <- rowSums(sp_df >= min_RA)
  
  select_species <- names(number_pass)[
    number_pass / ncol(sp_df) * 100 >= min_RA_percent
  ]
  
  sp_df_QCed <- sp_df[select_species, , drop = FALSE]
  
  cat("Number of species retained after QC:", length(select_species), "\n")
  
  sp_mat_full <- t(as.matrix(sp_df_QCed))
  rownames(sp_mat_full) <- as.character(rownames(sp_mat_full))
  
  ###########################
  ## 3.2 Prepare metadata
  ###########################
  meta <- metadata %>%
    mutate(
      gid_wgs   = as.character(gid_wgs),
      subjectID = as.character(subjectID)
    ) %>%
    filter(gid_wgs %in% rownames(sp_mat_full)) %>%
    filter(!is.na(subjectID), !is.na(age_at_collection)) %>%
    filter(age_at_collection <= 1000)
  
  meta <- meta %>%
    group_by(subjectID) %>%
    filter(n() >= min_samples) %>%
    ungroup()
  
  cat(
    "Number of subjects with at least ", min_samples,
    " samples within the first 1000 days: ",
    n_distinct(meta$subjectID), "\n",
    sep = ""
  )
  
  if (nrow(meta) == 0) {
    stop("No available samples after filtering by age and minimum sample count.")
  }
  
  ###########################
  ## 3.3 Prepare species abundance matrix
  ###########################
  sp_mat <- sp_mat_full[meta$gid_wgs, , drop = FALSE]
  
  row_sum <- rowSums(sp_mat)
  keep <- row_sum > 0
  
  sp_mat <- sp_mat[keep, , drop = FALSE]
  meta <- meta[keep, ]
  
  cat("Number of samples retained after removing all-zero samples:", nrow(meta), "\n")
  
  ###########################
  ## 3.4 Calculate Bray-Curtis distance to baseline
  ###########################
  bc <- vegdist(sp_mat, method = "bray")
  bc_mat <- as.matrix(bc)
  
  baseline_df <- meta %>%
    group_by(subjectID) %>%
    slice_min(age_at_collection, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(subjectID, baseline_gid = gid_wgs)
  
  meta2 <- meta %>%
    left_join(baseline_df, by = "subjectID")
  
  meta2$bc_dists_init <- mapply(
    function(sample_gid, base_gid) {
      bc_mat[sample_gid, base_gid]
    },
    sample_gid = meta2$gid_wgs,
    base_gid   = meta2$baseline_gid
  )
  
  ###########################
  ## 3.5 Create wide trajectory matrix
  ###########################
  meta2 <- meta2 %>%
    filter(gid_wgs != baseline_gid)
  
  traj_long <- meta2 %>%
    select(subjectID, age_at_collection, bc_dists_init) %>%
    rename(
      ID  = subjectID,
      Day = age_at_collection
    )
  
  traj_wide <- traj_long %>%
    pivot_wider(
      id_cols = ID,
      names_from = Day,
      values_from = bc_dists_init
    )
  
  day_cols <- colnames(traj_wide)[-1]
  sorted_day_cols <- as.character(sort(as.numeric(day_cols)))
  
  traj_wide <- traj_wide[, c("ID", sorted_day_cols)]
  time_points <- as.numeric(sorted_day_cols)
  
  ###########################
  ## 3.6 Cluster trajectories
  ###########################
  step1 <- Step1Measures(traj_wide, Time = time_points, ID = TRUE)
  step2 <- Step2Selection(step1)
  step3 <- Step3Clusters(step2, nclusters = nclusters)
  
  cluster.df <- step3$partition %>%
    rename(
      subjectID    = ID,
      traj_cluster = Cluster
    )
  
  ###########################
  ## 3.7 Merge cluster assignments
  ###########################
  meta_with_cluster <- meta2 %>%
    left_join(cluster.df, by = "subjectID")
  
  ###########################
  ## 3.8 Generate trajectory plot
  ###########################
  plot_df <- meta_with_cluster %>%
    filter(
      !is.na(traj_cluster),
      !is.na(bc_dists_init),
      !is.na(age_at_collection)
    ) %>%
    mutate(traj_cluster = factor(traj_cluster))
  
  p_traj <- ggplot(
    plot_df,
    aes(x = age_at_collection, y = bc_dists_init)
  ) +
    # geom_line(
    #   aes(group = subjectID, color = traj_cluster),
    #   alpha = 0.12,
    #   linewidth = 0.25
    # ) +
    geom_smooth(
      aes(color = traj_cluster, fill = traj_cluster),
      method = "loess",
      se = TRUE,
      linewidth = 1,
      alpha = 0.25
    ) +
    labs(
      x = "Time (Days)",
      y = "Bray-Curtis distance from the initial timepoint",
      color = "Trajectory cluster",
      fill = "Trajectory cluster"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank()
    )
  
  ###########################
  ## 3.9 Return outputs
  ###########################
  invisible(list(
    meta_with_cluster = meta_with_cluster,
    traj_wide         = traj_wide,
    time_points       = time_points,
    step1             = step1,
    step2             = step2,
    step3             = step3,
    cluster_df        = cluster.df,
    plot_df           = plot_df,
    plot              = p_traj
  ))
}

###########################
## 4. Run trajectory clustering
###########################
res1000 <- run_bray_traj_pipeline(
  Species     = Species,
  metadata    = metadata,
  min_samples = 4,
  nclusters   = 3
)

###########################
## 5. Customize plot colors
###########################
cluster_n <- res1000$cluster_df %>%
  count(traj_cluster)

cluster_labels <- setNames(
  paste0(cluster_n$traj_cluster, " (n = ", cluster_n$n, ")"),
  cluster_n$traj_cluster
)

p_traj <- res1000$plot +
  scale_color_manual(
    values = c("#4c956c", "#457b9d", "#e26d5c"),
    labels = cluster_labels
  ) +
  scale_fill_manual(
    values = c("#4c956c", "#457b9d", "#e26d5c"),
    labels = cluster_labels
  ) +
  theme_cowplot()

p_traj

###########################
## 6. Save plot
###########################
ggsave(
  filename = "~/ddong/TEDDY_project/figures/supplement/SupFig6_DIABIMMUNE_trajectory.pdf",
  plot = p_traj,
  width = 6,
  height = 4
)
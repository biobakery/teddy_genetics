############################################################
## Supplementary Figure 1
## Determination of the optimal number of microbiome
## maturational patterns (k) for trajectory clustering.
##
## a) Elbow method
## b) Gap statistic
## c) Calinski-Harabasz index
############################################################

###########################
## 0. Load packages
###########################

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(traj)
  library(cluster)
  library(ggplot2)
  library(cowplot)
}))

set.seed(1)

###########################
## 1. Set input/output paths
###########################

genetics_file <- "~/ddong/TEDDY_project/data/genetics.bed/filter/teddy_out8.qc.PCA_PC1-PC5.tsv"
metadata_file <- "~/ddong/TEDDY_project/output/traj_cluster_results.Rdata"
output_dir <- "~/ddong/TEDDY_project/figures/supplement/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

###########################
## 2. Load and preprocess data
###########################

load(metadata_file)
metadata <- metadata_cluster

geneticpc <- fread(genetics_file) %>%
  rename(subject_id = 1)

allphenotype <- metadata %>%
  inner_join(geneticpc, by = "subject_id") %>%
  filter(
    Days <= 800,
    (is.na(T1D) | T1D == FALSE),
    (is.na(AB) | AB == FALSE),
    !is.na(bc_dists_init)
  ) %>%
  group_by(subject_id) %>%
  filter(n() > 3) %>%
  ungroup()

allphenotype.df <- allphenotype %>%
  dcast(subject_id ~ Days,
        value.var = "bc_dists_init") %>%
  rename(ID = subject_id)

###########################
## 3. Generate trajectory feature matrix
###########################

step1 <- Step1Measures(
  allphenotype.df,
  Time = as.numeric(colnames(allphenotype.df)[-1]),
  ID = TRUE
)

step2 <- Step2Selection(step1)

get_measures <- function(s2, s1) {
  
  cand <- list(
    s2$selection,
    s2$measures,
    s2$data,
    s1$measures
  )
  
  for (m in cand) {
    
    if (is.null(m)) next
    
    df <- as.data.frame(m)
    
    num <- df[, vapply(df, is.numeric, logical(1)),
              drop = FALSE]
    
    idlike <- vapply(
      num,
      function(x)
        all(x == round(x), na.rm = TRUE) &&
        !any(duplicated(x)),
      logical(1)
    )
    
    if (any(idlike))
      num <- num[, !idlike, drop = FALSE]
    
    if (ncol(num) >= 2 && nrow(num) > 10)
      return(as.matrix(num))
  }
  
  stop("Cannot locate trajectory measures.")
}

M <- get_measures(step2, step1)

M <- scale(M)
M <- M[complete.cases(M), , drop = FALSE]

n <- nrow(M)

###########################
## 4. Calculate clustering diagnostics
###########################

k_max <- 10
ks <- 1:k_max

wss <- numeric(length(ks))
ch <- rep(NA_real_, length(ks))

for (i in seq_along(ks)) {
  
  k <- ks[i]
  
  km <- kmeans(
    M,
    centers = k,
    nstart = 25,
    iter.max = 100
  )
  
  wss[i] <- km$tot.withinss
  
  if (k > 1) {
    ch[i] <- (km$betweenss / (k - 1)) /
      (km$tot.withinss / (n - k))
  }
}

###########################
## 5. Calculate Gap statistic
###########################

gap_kmax <- 8

gap <- clusGap(
  M,
  FUN = kmeans,
  nstart = 25,
  K.max = gap_kmax,
  B = 100,
  iter.max = 100
)

###########################
## 6. Prepare summary tables
###########################

gap.df <- data.frame(
  k = seq_len(gap_kmax),
  gap = gap$Tab[, "gap"],
  SE = gap$Tab[, "SE.sim"]
)

elbow.df <- data.frame(
  k = ks,
  wss = wss
)

ch.df <- data.frame(
  k = ks,
  ch = ch
) %>%
  filter(!is.na(ch))

###########################
## 7. Determine optimal k
###########################

## Elbow (maximum curvature)

elbow_best_k <-
  elbow.df$k[
    which.max(diff(diff(elbow.df$wss))) + 1
  ]

## Gap statistic (Tibshirani 1-SE rule)

gap_best_k <- maxSE(
  gap$Tab[, "gap"],
  gap$Tab[, "SE.sim"],
  method = "Tibs2001SEmax"
)

## Calinski-Harabasz

ch_best_k <-
  ch.df$k[
    which.max(ch.df$ch)
  ]

best_k.df <- data.frame(
  method = c(
    "Elbow",
    "Gap statistic",
    "Calinski-Harabasz"
  ),
  best_k = c(
    elbow_best_k,
    gap_best_k,
    ch_best_k
  )
)

###########################
## 8. Generate figure panels
###########################

p_a <- ggplot(elbow.df, aes(k, wss)) +
  geom_line(color = "#3B7BB4") +
  geom_point(color = "#3B7BB4") +
  geom_vline(
    xintercept = elbow_best_k,
    linetype = "dashed",
    colour = "grey60"
  ) +
  scale_x_continuous(breaks = ks) +
  labs(
    title = paste0("Elbow Method: k = ", elbow_best_k),
    x = "Number of clusters k",
    y = "Total Within Sum of Square"
  ) +
  theme_cowplot()

p_b <- ggplot(gap.df, aes(k, gap)) +
  geom_line(color = "#3B7BB4") +
  geom_point(color = "#3B7BB4") +
  geom_errorbar(
    aes(ymin = gap - SE, ymax = gap + SE),
    width = 0.2,
    color = "#3B7BB4"
  ) +
  geom_vline(
    xintercept = gap_best_k,
    linetype = "dashed",
    colour = "#8Fb8d8"
  ) +
  scale_x_continuous(breaks = seq_len(gap_kmax)) +
  labs(
    title = paste0("Gap Statistic: k = ", gap_best_k),
    x = "Number of clusters k",
    y = "Gap statistic (k)"
  ) +
  theme_cowplot()

p_c <- ggplot(ch.df, aes(k, ch)) +
  geom_line(color = "black") +
  geom_point(color = "black") +
  geom_vline(
    xintercept = ch_best_k,
    linetype = "dashed",
    colour = "grey60"
  ) +
  scale_x_continuous(breaks = seq(2, k_max, 2)) +
  labs(
    title = paste0("Calinski-Harabasz: k = ", ch_best_k),
    x = "Number of clusters K",
    y = "Calinski-Harabasz Index"
  ) +
  theme_cowplot()
###########################
## 9. Assemble figure
###########################

fig <- plot_grid(
  plot_grid(
    p_a,
    p_b,
    labels = c("a", "b"),
    nrow = 1
  ),
  plot_grid(
    p_c,
    NULL,
    labels = c("c", ""),
    nrow = 1,
    rel_widths = c(1, 1)
  ),
  nrow = 2
)

###########################
## 10. Save outputs
###########################

out_pdf <- file.path(
  output_dir,
  "SupFig3.pdf"
)

ggsave(
  out_pdf,
  fig,
  width = 10,
  height = 7
)

fwrite(
  elbow.df,
  file.path(output_dir,
            "SupFig3_elbow_wss.tsv"),
  sep = "\t"
)

fwrite(
  gap.df,
  file.path(output_dir,
            "SupFig3_gap_statistic.tsv"),
  sep = "\t"
)

fwrite(
  ch.df,
  file.path(output_dir,
            "SupFig3_calinski_harabasz.tsv"),
  sep = "\t"
)

fwrite(
  best_k.df,
  file.path(output_dir,
            "SupFig3_best_k_summary.tsv"),
  sep = "\t"
)

message("Wrote: ", out_pdf)
print(best_k.df)
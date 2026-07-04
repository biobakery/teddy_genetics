############################################################
## Supplementary Figure 5
## Sensitivity analyses of microbiome maturation trajectories
## using alternative cluster numbers (k = 2 and k = 4).
############################################################

###########################
## 0. Load packages
###########################
library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)
library(traj)
library(nnet)
library(survival)
library(survminer)
library(ggplot2)
library(cowplot)
library(patchwork)



###########################
## 1. Set input/output paths
###########################
set.seed(2025)
setwd("~/ddong/TEDDY_project/")

genetics_file <- "./data/genetics.bed/filter/teddy.qc.PCA_PC1-PC5.tsv"
metadata_file <- "./output/traj_cluster_results.Rdata"
output_dir    <- "./figures/supplement/"

###########################
## 2. Load data and build clustering input
###########################

load(metadata_file)   # -> metadata_cluster (long), Cluster_results (subject-level)

geneticpc <- fread(genetics_file) %>%
  rename(subject_id = 1)

allphenotype <- metadata_cluster %>%
  inner_join(geneticpc, by = "subject_id") %>%
  filter(
    Days <= 800,
    (is.na(T1D) | T1D == FALSE),
    (is.na(AB)  | AB  == FALSE),
    !is.na(bc_dists_init)
  ) %>%
  group_by(subject_id) %>%
  filter(n() > 4) %>%
  ungroup()

allphenotype.df <- allphenotype %>%
  dcast(subject_id ~ Days, value.var = "bc_dists_init") %>%
  rename(ID = subject_id)

###########################
## 3. Trajectory measures (computed once, re-partitioned per k)
###########################

step1 <- Step1Measures(
  allphenotype.df,
  Time = as.numeric(colnames(allphenotype.df)[-1]),
  ID   = TRUE
)
step2 <- Step2Selection(step1)

###########################
## 4. Per-k analysis function
###########################

palettes <- list(
  `2` = c("#457b9d", "#e26d5c"),
  `4` = c("#4c956c", "#457b9d", "#e26d5c", "#e9c46a")
)

ps_covariates <- c("Clinical_Center", "Sex", "FDR", "delivery",
                   "PC1", "PC2", "PC3", "PC4", "PC5", "first_day",
                   "Antibiotics", "Probiotic",
                   "time_to_brstfed_stop", "age_startsolid")

analyse_k <- function(k) {
  
  ## 4.1 Re-partition into k clusters
  step3 <- Step3Clusters(step2, nclusters = k)
  
  cluster.df <- step3$partition %>%
    rename(subject_id = ID) %>%
    mutate(Cluster = factor(paste0("C", Cluster)))
  
  ## 4.2 Panel a input: long trajectory + new cluster labels
  traj.df <- metadata_cluster %>%
    inner_join(cluster.df, by = "subject_id")
  
  p_traj <- ggplot(traj.df,
                   aes(x = Days, y = bc_dists_init, colour = Cluster)) +
    geom_smooth(method = "loess") +
    scale_colour_manual(values = palettes[[as.character(k)]]) +
    labs(title = paste0("k = ", k),
         x = "Days", y = "Bray-Curtis to baseline") +
    theme_cowplot()
  
  ## 4.3 Subject-level covariates + new cluster (drop primary Cluster)
  surv.df <- Cluster_results %>%
    select(-any_of("Cluster")) %>%
    inner_join(cluster.df, by = "subject_id")
  
  model_vars <- c("IA_T1D_time", "AB_T1D", ps_covariates)
  surv.df <- surv.df %>%
    filter(if_all(all_of(model_vars), ~ !is.na(.)))
  
  ## 4.4 Inverse-probability (propensity) weights, matching Fig2b's design
  ps_model <- multinom(
    as.formula(paste("Cluster ~", paste(ps_covariates, collapse = " + "))),
    data = surv.df, trace = FALSE
  )
  
  ps_mat <- predict(ps_model, type = "probs")
  if (is.null(dim(ps_mat))) {                       # k == 2 => vector
    ps_mat <- cbind(1 - ps_mat, ps_mat)
    colnames(ps_mat) <- levels(surv.df$Cluster)
  }
  assigned <- as.character(surv.df$Cluster)
  p_assigned <- ps_mat[cbind(seq_len(nrow(ps_mat)),
                             match(assigned, colnames(ps_mat)))]
  surv.df$weight <- 1 / p_assigned                  # IPTW
  
  ## 4.5 Weighted Cox: cumulative incidence (stratified) + HR (covariate)
  cox_strata <- coxph(Surv(IA_T1D_time, AB_T1D) ~ strata(Cluster),
                      weights = weight, data = surv.df)
  
  p_cuminc <- ggsurvplot(
    survfit(cox_strata), data = surv.df, fun = "event",
    palette = palettes[[as.character(k)]],
    conf.int = TRUE, xlim = c(0, 1825), ylim = c(0, 1),
    xlab = "Time", ylab = "Cumulative event rate",
    title = paste0("k = ", k), ggtheme = theme_cowplot()
  )$plot
  
  cox_hr <- coxph(Surv(IA_T1D_time, AB_T1D) ~ Cluster,
                  weights = weight, data = surv.df)
  hr.df <- data.frame(
    k        = k,
    term     = rownames(summary(cox_hr)$conf.int),
    HR       = summary(cox_hr)$conf.int[, "exp(coef)"],
    CI_low   = summary(cox_hr)$conf.int[, "lower .95"],
    CI_high  = summary(cox_hr)$conf.int[, "upper .95"],
    p        = summary(cox_hr)$coefficients[, "Pr(>|z|)"],
    row.names = NULL
  )
  
  list(traj = p_traj, cuminc = p_cuminc, hr = hr.df,
       clusters = cluster.df)
}

###########################
## 5. Run for k = 2 and k = 4
###########################

res2 <- analyse_k(2)
res4 <- analyse_k(4)

###########################
## 6. Assemble figure (a: trajectories, b: cumulative incidence)
###########################

fig <- (res2$traj   | res4$traj) /
  (res2$cuminc | res4$cuminc) +
  plot_annotation(tag_levels = list(c("a", "", "b", "")))

###########################
## 7. Save outputs
###########################
out_pdf <- file.path(output_dir, "SupFig5.pdf")
ggsave(out_pdf, fig, width = 11, height = 9)
print(hr.df)

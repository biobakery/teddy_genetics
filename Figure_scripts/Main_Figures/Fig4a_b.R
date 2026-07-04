# ============================================================
# Figure 4a–4b: Genetic PC3 modification of maturation-pattern risk
# (a) Cross-pattern marginal HR vs Early Matured at the SAME PC3 value
# (b) Forest plot of Late Matured/Early Plateaued vs Early Matured within PC3 strata 
# ============================================================

## Set project root
setwd("~/ddong/TEDDY_project/")

###########################
## 0. Load packages
###########################
library(data.table)
library(dplyr)
library(survival)
library(ggplot2)
library(cowplot)

###########################
## 1. Load data
###########################
load("./output/traj_cluster_results.Rdata")          # Load clustering results
metadata <- fread("./data/metadata/metadata_subject.txt")     # subject-level covariates

## Merge clustering results with metadata
Cluster_results_meta <- Cluster_results %>%
  left_join(metadata)

###########################
## 2. Harmonize factor coding and reference levels
###########################
## Cluster 1 = Early Matured (reference)
ref_cluster <- "1"

Cluster_results_meta <- Cluster_results_meta %>%
  mutate(
    Cluster = factor(Cluster),
    Cluster = relevel(Cluster, ref = ref_cluster),
    Sex = factor(Sex),
    Clinical_Center = factor(Clinical_Center),
    FDR = factor(FDR),
    delivery_class = factor(delivery_class),
    Antibiotics = factor(Antibiotics),
    Probiotic = factor(Probiotic)
  )

############################################################
## Figure 4a
## Cross-pattern marginal HR vs Early Matured at SAME PC3
############################################################

###########################
## 3. Fit Cox model with PC × Cluster interactions
###########################
fit <- coxph(
  Surv(AB_T1D_time_to0, AB_T1D) ~
    Clinical_Center + Sex + FDR + delivery_class +
    Antibiotics + Probiotic + time_to_brstfed_stop + age_startsolid +
    PC1 + PC2 + PC3 + PC4 + PC5 + Cluster +
    PC1:Cluster + PC2:Cluster + PC3:Cluster + PC4:Cluster + PC5:Cluster,
  data = Cluster_results_meta
)

###########################
## 4. Build prediction grid across PC3
###########################
pc3_rng <- range(Cluster_results_meta$PC3, na.rm = TRUE)
pc3_seq <- seq(pc3_rng[1], pc3_rng[2], by = 0.05)

grid <- expand.grid(
  PC3 = pc3_seq,
  Cluster = levels(Cluster_results_meta$Cluster)
)

###########################
## 5. Fix other covariates at mean / modal values
###########################
ref_row <- Cluster_results_meta %>%
  summarise(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    PC4 = mean(PC4, na.rm = TRUE),
    PC5 = mean(PC5, na.rm = TRUE),
    
    Sex = names(sort(table(Sex), decreasing = TRUE))[1],
    Clinical_Center = names(sort(table(Clinical_Center), decreasing = TRUE))[1],
    FDR = names(sort(table(FDR), decreasing = TRUE))[1],
    delivery_class = names(sort(table(delivery_class), decreasing = TRUE))[1],
    Antibiotics = names(sort(table(Antibiotics), decreasing = TRUE))[1],
    Probiotic = names(sort(table(Probiotic), decreasing = TRUE))[1],
    
    time_to_brstfed_stop = mean(time_to_brstfed_stop, na.rm = TRUE),
    age_startsolid = mean(age_startsolid, na.rm = TRUE)
  )

grid <- cbind(grid, ref_row[rep(1, nrow(grid)), ])

## Ensure factor levels align with model frame
grid <- grid %>%
  mutate(
    Cluster = factor(Cluster, levels = levels(Cluster_results_meta$Cluster)),
    Sex = factor(Sex, levels = levels(Cluster_results_meta$Sex)),
    Clinical_Center = factor(Clinical_Center, levels = levels(Cluster_results_meta$Clinical_Center)),
    FDR = factor(FDR, levels = levels(Cluster_results_meta$FDR)),
    delivery_class = factor(delivery_class, levels = levels(Cluster_results_meta$delivery_class)),
    Antibiotics = factor(Antibiotics, levels = levels(Cluster_results_meta$Antibiotics)),
    Probiotic = factor(Probiotic, levels = levels(Cluster_results_meta$Probiotic))
  )

###########################
## 6. Predict linear predictor and compute marginal HR vs Cluster 1 at SAME PC3
###########################
grid$lp <- as.numeric(predict(fit, newdata = grid, type = "lp"))

ref_early <- grid %>%
  filter(Cluster == ref_cluster) %>%
  select(PC3, lp_ref = lp)

grid_plot <- grid %>%
  left_join(ref_early, by = "PC3") %>%
  mutate(HR = exp(lp - lp_ref))

###########################
## 7. Plot Fig4a
###########################
cluster_labels <- c("1" = "Early Matured",
                    "2" = "Late Matured",
                    "3" = "Early Plateaued")

cluster_colors <- c("1" = "#4c956c",
                    "2" = "#457b9d",
                    "3" = "#e26d5c")

pdf("./figures/main/Fig4a_MarginPC3_crosspattern.pdf", width = 6, height = 4)
ggplot(grid_plot, aes(x = PC3, y = HR, color = Cluster)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = cluster_colors, labels = cluster_labels) +
  labs(
    x = "Genetic PC3 (z-score)",
    y = "Risk ratio per unit increase\n(Seroconversion or T1D)"
  ) +
  theme_cowplot()
dev.off()

############################################################
## Figure 4b
## Forest plot: Late Matured/Early Plateaued vs Early Matured within PC3 strata (Q1/Q2)
############################################################

###########################
## 8. Define PC3 strata (median split: Q1/Q2)
###########################
q_pc3 <- as.numeric(stats::quantile(Cluster_results_meta$PC3, probs = c(0, 0.5, 1), na.rm = TRUE))

Cluster_results_meta <- Cluster_results_meta %>%
  mutate(
    PC3_half = cut(
      PC3,
      breaks = q_pc3,
      labels = c("Q1", "Q2"),
      include.lowest = TRUE
    ),
    PC3_half = factor(PC3_half, levels = c("Q1", "Q2")),
    Cluster = relevel(factor(Cluster), ref = "1")
  )

###########################
## 9. Fit Cox model within each stratum and extract Late Matured/Early Plateaued effects
###########################
fit_cox_stratum <- function(dat, stratum_label) {
  dat2 <- dat %>% filter(PC3_half == stratum_label)
  
  m <- coxph(
    Surv(AB_T1D_time_to0, AB_T1D) ~
      Cluster +
      Clinical_Center + Sex + FDR +
      delivery_class +
      time_to_brstfed_stop + age_startsolid +
      PC1 + PC2 + PC3 + PC4 + PC5 +
      Antibiotics + Probiotic,
    data = dat2
  )
  
  s  <- summary(m)
  ci <- as.data.frame(s$conf.int)
  ci$term <- rownames(ci)
  
  out <- ci %>%
    filter(term %in% c("Cluster2", "Cluster3")) %>%
    transmute(
      Side  = ifelse(stratum_label == "Q1", "Low PC3 Level", "High PC3 Level"),
      term  = term,
      HR    = as.numeric(`exp(coef)`),
      Lower = as.numeric(`lower .95`),
      Upper = as.numeric(`upper .95`)
    )
  
  out
}

res_q1 <- fit_cox_stratum(Cluster_results_meta, "Q1")
res_q2 <- fit_cox_stratum(Cluster_results_meta, "Q2")
res    <- bind_rows(res_q1, res_q2)

###########################
## 10. Map cluster terms to figure comparison labels and set display order
###########################
res <- res %>%
  mutate(
    Comparison = case_when(
      term == "Cluster2" ~ "Late Matured vs Early Matured",
      term == "Cluster3" ~ "Early Plateaued vs Early Matured",
      TRUE ~ term
    ),
    Side = factor(Side, levels = c("High PC3 Level", "Low PC3 Level")),
    Comparison = factor(
      Comparison,
      levels = c("Late Matured vs Early Matured", "Early Plateaued vs Early Matured")
    )
  )

###########################
## 11. Forest plot 
###########################

pdf("./figures/main/Fig4b_PC3_Q1Q2_forest.pdf", width = 5, height = 3)
ggplot(res, aes(x = HR, y = Side, fill = Side)) +
  scale_fill_manual(values = c("#978164", "#498985")) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(shape = 21, size = 3) +
  facet_wrap(~ Comparison, scales = "free_y", nrow = 2) +
  scale_x_continuous(
    trans  = "log10",
    breaks = c(0.5, 1, 2, 5),
    labels = c("0.5", "1", "2", "5")
  ) +
  labs(x = "Risk Ratio for Seroconversion or T1D (95% CI)", y = NULL) +
  theme_cowplot() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
dev.off()
fwrite(as.data.table(res), "./output/Fig4b_PC3_Q1Q2_forest_numbers.tsv", sep = "\t")

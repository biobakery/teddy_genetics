############################################################
## Figure 5b: Infant & Toddler examples (function-based)
## Representative genetic PC × time associations 
############################################################

## Set project root
setwd("~/ddong/TEDDY_project/")

###########################
## 0. Load packages
###########################
library(tidyverse)
library(reshape2)
library(ggplot2)
library(viridis)
library(scales)
library(lmerTest)
library(broom.mixed)

###########################
## 1. Utility functions
###########################

load_shared_inputs <- function(genetics_path, metadata_path) {
  snps <- read.csv(genetics_path, sep = "\t", header = TRUE, check.names = FALSE)
  colnames(snps)[1] <- "subject_id"
  
  metadata <- read.csv(metadata_path, sep = "\t", header = TRUE)
  
  list(snps = snps, metadata = metadata)
}

load_microbiome_log10 <- function(micro_path) {
  microbiome_melt <- read.csv(micro_path, sep = "\t", header = TRUE) %>%
    mutate(sample_id = factor(sample_id)) %>%
    melt() %>%
    dplyr::filter(value > 0)
  
  eps <- min(microbiome_melt$value)
  
  microbiome <- read.csv(micro_path, sep = "\t", header = TRUE) %>%
    mutate(sample_id = factor(sample_id)) %>%
    melt() %>%
    mutate(value = log10(value + eps)) %>%
    dcast(sample_id ~ variable, value.var = "value")
  
  list(microbiome = microbiome, eps = eps)
}

make_df_sig <- function(model_stats, p_method = "bonferroni") {
  model_stats %>%
    filter(!Predictor %in% c("Country", "Day", "Probiotic")) %>%
    group_by(Predictor) %>%
    mutate(P_adj = p.adjust(P, method = p_method)) %>%
    filter(P_adj < 0.05) %>%
    filter(grepl("PC", Predictor)) %>%
    mutate(
      Type = case_when(
        !grepl(":", Predictor) ~ "genetics",
        grepl(":", Predictor)  ~ "genetics:time"
      ),
      Coe = Coefficient,
      Coefficient = case_when(
        Coefficient > 0 ~  1,
        Coefficient < 0 ~ -1
      ),
      value = Coefficient * -log10(P_adj),
      MicroFeat = gsub("lmer\\(", "", Model),
      MicroFeat = gsub(" .*", "", MicroFeat),
      MicroFeat = gsub("_noname_unclassified", "_?_?", MicroFeat),
      MicroFeat = gsub("_unclassified", "_?", MicroFeat),
      pcn = gsub("\\:.*", "", Predictor)
    ) %>%
    ungroup()
}

make_ann_sig_interaction <- function(df_sig) {
  df_sig %>%
    filter(grepl(":", Predictor)) %>%
    mutate(
      pcn = gsub("\\:.*", "", Predictor),
      P_adj = paste0("q=", scientific(P_adj))
    ) %>%
    select(MicroFeat, pcn, P_adj, Coefficient)
}

augment_main_effect <- function(model_text, dat, augment_engine = c("broom", "broom.mixed")) {
  augment_engine <- match.arg(augment_engine)
  
  model <- eval(parse(text = model_text),
                envir = list2env(list(dat = dat), parent = environment()))
  
  aug <- if (augment_engine == "broom.mixed") {
    broom.mixed::augment(model)
  } else {
    broom::augment(model)
  }
  
  aug %>%
    rename(PC = 2) %>%
    mutate(
      Model = model_text,
      MicroFeat = gsub("lmer\\(| .*", "", Model),
      MicroFeat = gsub("_noname_unclassified", "_?_?", MicroFeat),
      MicroFeat = gsub("_unclassified", "_?", MicroFeat),
      pcn = gsub("*.\\~  | \\*.*", "", Model),
      pcn = gsub(".*P", "P", pcn),
      Quartile = ntile(Day, 4)
    ) %>%
    select(Subject, MicroFeat, pcn, Quartile, .fitted, PC)
}

augment_interaction <- function(model_text, dat) {
  model <- eval(parse(text = model_text),
                envir = list2env(list(dat = dat), parent = environment()))
  
  broom.mixed::augment(model) %>%
    rename(RelAb = 1, PC = 2) %>%
    mutate(
      Model = model_text,
      MicroFeat = gsub("lmer\\(| .*", "", Model),
      MicroFeat = gsub("_noname_unclassified", "_?_?", MicroFeat),
      MicroFeat = gsub("_unclassified", "_?", MicroFeat),
      pcn = gsub("*.\\~  | \\*.*", "", Model),
      pcn = gsub(".*P", "P", pcn),
      Quartile = ntile(PC, 4)
    ) %>%
    select(Subject, MicroFeat, pcn, Quartile, Day, .fitted, RelAb, PC)
}

plot_interaction <- function(combined, ann_sig, pdf_out, width, height) {
  pdf(pdf_out, width, height)
  
  print(ggplot(combined, aes(x = Day, y = .fitted)) +
          facet_wrap(~ MicroFeat + pcn, scales = "free", ncol = 5) +
          geom_smooth(aes(colour = factor(Quartile), fill = factor(Quartile)), method = "lm") +
          scale_colour_viridis(discrete = TRUE) +
          scale_fill_viridis(discrete = TRUE) +
          theme_bw() +
          theme(
            plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 12.5, face = "bold"),
            axis.text = element_text(size = 8.75),
            legend.title = element_text(size = 12.5, face = "bold", vjust = 0.75),
            legend.text = element_text(size = 12.5),
            strip.text = element_text(size = 8.75, face = "bold"),
            legend.position = "top"
          ) +
          labs(
            title = "Interaction effect(s)",
            colour = "PC Quartile",
            fill = "PC Quartile",
            y = expression(paste(bold(log[10] * " abundance (%)")))
          ) +
          geom_text(
            data = ann_sig,
            aes(x = Inf, y = Inf, label = P_adj),
            hjust = 1.5,
            vjust = 1.5
          ))
  
  dev.off()
}

###########################
## 2. Panel runner
###########################
run_panel <- function(panel_name,
                      micro_path,
                      genetics_path,
                      metadata_path,
                      model_stats_path,
                      p_method,
                      outdir,
                      pdf_out,
                      pdf_w,
                      pdf_h,
                      main_effect_augment_engine = c("broom", "broom.mixed")) {
  
  main_effect_augment_engine <- match.arg(main_effect_augment_engine)
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  shared <- load_shared_inputs(genetics_path, metadata_path)
  micro  <- load_microbiome_log10(micro_path)
  
  df_subjects <- merge(shared$metadata, shared$snps, by = "subject_id")
  
  dat <- merge(df_subjects, micro$microbiome, by = "sample_id") %>%
    rename(Subject = subject_id, Day = Days) %>%
    mutate(Subject = factor(Subject))
  
  model_stats <- read.csv(model_stats_path, sep = "\t", header = TRUE)
  
  df_sig <- make_df_sig(model_stats, p_method = p_method)
  
  df_pc <- df_sig %>% filter(!grepl(":", Type))
  df_pcxday <- df_sig %>% filter(grepl(":", Type))
  
  combined_main <- list()
  if (nrow(df_pc) > 0) {
    for (m in df_pc$Model) {
      combined_main <- rbind(combined_main, augment_main_effect(m, dat, augment_engine = main_effect_augment_engine))
    }
  }
  
  combined_int <- list()
  if (nrow(df_pcxday) > 0) {
    for (m in df_pcxday$Model) {
      combined_int <- rbind(combined_int, augment_interaction(m, dat))
    }
  }
  
  ann_sig_int <- make_ann_sig_interaction(df_sig)
  
  plot_interaction(combined_int, ann_sig_int, pdf_out, pdf_w, pdf_h)
  
  invisible(list(
    panel = panel_name,
    dat = dat,
    df_sig = df_sig,
    combined_main = combined_main,
    combined_int = combined_int,
    ann_sig_int = ann_sig_int
  ))
}

###########################
## 3. Run infant & toddler
###########################
res_infant <- run_panel(
  panel_name = "infant",
  micro_path = "./output/microbiome_infant.txt",
  genetics_path = "./data/genetics.bed/filter/teddy.qc.PCA_PC1-PC5.tsv",
  metadata_path = "./data/metadata/metadata_sample.txt",
  model_stats_path = "./output/model_anova_total_infant.txt",
  p_method = "bonferroni",
  outdir = "./output/micro_infant",
  pdf_out = "./figures/main/Fig5b_infant_interaction.pdf",
  pdf_w = 14,
  pdf_h = 10,
  main_effect_augment_engine = "broom.mixed"
)

res_toddler <- run_panel(
  panel_name = "toddler",
  micro_path = "./output/microbiome_toddler.txt",
  genetics_path = "./data/genetics.bed/filter/teddy.qc.PCA_PC1-PC5.tsv",
  metadata_path = "./data/metadata/metadata_sample.txt",
  model_stats_path = "./output/model_anova_total_toddler.txt",
  p_method = "bonferroni",
  outdir = "./output/micro_toddler",
  pdf_out = "./figures/main/Fig5b_toddler_interaction.pdf",
  pdf_w = 10,
  pdf_h = 5,
  main_effect_augment_engine = "broom.mixed"
)

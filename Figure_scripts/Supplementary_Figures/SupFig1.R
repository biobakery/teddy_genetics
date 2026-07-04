############################################################
## Supplementary Figure 1
## Per-participant distribution of metagenomic stool samples
## across age.
############################################################

###########################
## 0. Load packages
###########################
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)

###########################
## 1. Set input/output paths
###########################
set.seed(1)
setwd("~/ddong/TEDDY_project")

metadata_file <- "./output/preprocess/metadata_microbiome.Rdata"
output_dir    <- "./figures/supplement/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

###########################
## 2. Load sample-level metadata
##    (one row = one shotgun-sequenced sample)
###########################

load(metadata_file)   # -> object `metadata`

samples <- metadata %>%
  select(subject_id, sample_id, mgx_age, Country) %>%
  filter(!is.na(mgx_age)) %>%
  distinct()

###########################
## 3. Order participants by sampling density
###########################

subject_order <- samples %>%
  group_by(subject_id) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  arrange(n_samples, subject_id)

samples <- samples %>%
  mutate(subject_id = factor(subject_id,
                             levels = subject_order$subject_id))

###########################
## 4. Report cohort-level summary
###########################

summary.df <- data.frame(
  n_samples          = nrow(samples),
  n_participants     = nlevels(samples$subject_id),
  mean_per_subject   = mean(subject_order$n_samples),
  sd_per_subject     = sd(subject_order$n_samples)
)

message("Samples: ", summary.df$n_samples,
        " | Participants: ", summary.df$n_participants,
        " | mean/subject: ", round(summary.df$mean_per_subject, 1),
        " +/- ", round(summary.df$sd_per_subject, 1))

###########################
## 5. Generate figure
###########################

p <- ggplot(samples,
            aes(x = mgx_age, y = subject_id, colour = Country)) +
  geom_point(size = 0.35, alpha = 0.7) +
  scale_colour_manual(values = c("#556fa1", "#c04d4d", "#e6ad5b", "#a26162")) +
  scale_x_continuous(breaks = seq(0, 2000, 500)) +
  labs(
    x = "Age (days)",
    y = "Participant (sorted by sampling density)",
    colour = "Country"
  ) +
  theme_cowplot() +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

###########################
## 6. Save outputs
###########################
out_pdf <- file.path(output_dir, "SupFig1.pdf")
ggsave(out_pdf, p, width = 7, height = 8)
############################################################
## Supplementary Figure 4
## Longitudinal trajectories of selected metabolic EC groups
##    Panel A: LOESS-smoothed abundance trajectories
##    Panel B: Wilcoxon effect sizes by age bin
############################################################


############################################################
## 1. Load  packages
############################################################

library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(mgcv)
library(emmeans)
library(tibble)
library(ggplot2)
library(cowplot)
library(stringr)
library(scales)
library(viridis)
library(rstatix)
library(patchwork)


############################################################
## 2. Set working directory
############################################################

setwd("~/ddong/TEDDY_project/")


############################################################
## 3. Load EC abundance table and EC annotation
############################################################

ec.df <- fread("./output/ec_all_qc_1e6.txt")
ec_name <- fread("./output/ec_name.txt")

############################################################
## 4. Load manually selected ECs
############################################################

ec_select <- fread(
  "./data/external/select_EC.csv",
  header = TRUE
)

ec_select <- ec_select %>%
  left_join(ec_name)

ec_select <- ec_select %>%
  mutate(feature = paste("X", EC, sep = ""))


############################################################
## 5. Define developmental age bins
############################################################

binnames <- c(
  sprintf("[%d,%d]", 100, 200),
  sprintf("(%d,%d]", seq(200, 700, by = 100),
          seq(300, 800, by = 100))
)


############################################################
## 6. Load MaAsLin results for EC associations
############################################################

data <- list()
all_EC.df <- data.frame()

setwd("./output/maaslin_EC/")

## Comparisons using Cluster 1 as reference
for (i in seq_along(binnames)) {
  b_n <- binnames[i]
  
  data[[i]] <- fread(
    paste("./", paste("Cluster", b_n, sep = "_"),
          "/all_results.tsv", sep = "")
  )
  
  data[[i]] <- data[[i]] %>%
    filter(metadata == "Cluster") %>%
    mutate(
      Type = paste("Cluster1_Cluster", value, sep = ""),
      bin = binnames[i]
    )
  
  all_EC.df <- rbind(all_EC.df, data[[i]])
}

## Comparisons using Cluster 2 as reference
data2 <- list()

for (i in seq_along(binnames)) {
  b_n <- binnames[i]
  
  data2[[i]] <- fread(
    paste("./", paste("Cluster_base2_", b_n, sep = "_"),
          "/all_results.tsv", sep = "")
  )
  
  data2[[i]] <- data2[[i]] %>%
    filter(metadata == "Cluster") %>%
    mutate(
      Type = paste("Cluster2_Cluster", value, sep = ""),
      bin = binnames[i]
    )
  
  all_EC.df <- rbind(all_EC.df, data2[[i]])
}

## Check available pairwise comparison types
table(all_EC.df$Type)


############################################################
## 7. Add EC annotation to MaAsLin results
############################################################

ec_name <- ec_name %>%
  mutate(feature = paste("X", EC, sep = ""))

all_EC.df <- all_EC.df %>%
  left_join(ec_name)

## Remove redundant comparison: Cluster2 versus Cluster1
all_EC.df <- all_EC.df %>%
  filter(Type != "Cluster2_Cluster1")

## Optional: save combined MaAsLin EC results
# save(all_EC.df, file = "~/Desktop//allEC_result.Rdata")


############################################################
## 8. Select ECs with significant MaAsLin associations
############################################################

sig <- all_EC.df %>%
  filter(qval < 0.25)

alldata_select_EC <- all_EC.df %>%
  filter(EC %in% sig$EC)

alldata_select_EC$bin <- factor(
  alldata_select_EC$bin,
  levels = binnames
)

alldata_select_EC <- alldata_select_EC %>%
  mutate(
    sig = case_when(
      qval < 0.001 ~ "***",
      qval < 0.05  ~ "**",
      qval < 0.25  ~ "*"
    )
  )


############################################################
## 9. Load KEGG Orthology annotation
############################################################

KO <- fread("./data/external/KO_03202022.csv")

KO <- KO %>%
  mutate(
    EC_number = gsub("EC", "", ec),
    feature = paste("X", EC_number, sep = "")
  )


############################################################
## 10. Annotate significant ECs with KEGG information
############################################################

sig_EC <- alldata_select_EC %>%
  select(feature, EC_name, EC) %>%
  unique()

sig_EC <- sig_EC %>%
  left_join(KO)

sig_EC <- sig_EC %>%
  select(
    cat1, cat2, cat3, cat4, ko,
    feature, EC, ec, EC_number,
    one_of(colnames(sig_EC))
  )


############################################################
## 11. Select broad KEGG categories of interest
############################################################

KO_select <- KO %>%
  filter(
    cat2 %in% c(
      "09105 Amino acid metabolism",
      "09106 Metabolism of other amino acids",
      "09108 Metabolism of cofactors and vitamins"
    ) |
      cat3 %in% c("00052 Galactose metabolism [PATH:ko00052]")
  )


############################################################
## 12. Reconstruct significant EC table with KEGG annotation
############################################################

sig_EC <- alldata_select_EC %>%
  select(feature, EC_name, EC) %>%
  unique()

sig_EC <- sig_EC %>%
  left_join(KO)

sig_EC <- sig_EC %>%
  select(
    cat1, cat2, cat3, cat4, ko,
    feature, EC, ec, EC_number,
    one_of(colnames(sig_EC))
  )


############################################################
## 13. Add manually curated EC group information
############################################################

sig_EC <- sig_EC %>%
  left_join(ec_select)

sig_EC_select <- sig_EC %>%
  filter(
    cat2 %in% c(
      "09105 Amino acid metabolism",
      "09106 Metabolism of other amino acids",
      "09108 Metabolism of cofactors and vitamins"
    ) |
      cat3 %in% c("00052 Galactose metabolism [PATH:ko00052]") |
      !is.na(EC_group_top)
  )

## Remove duplicated EC-pathway combinations
sig_EC_select <- sig_EC_select %>%
  mutate(paste = paste(cat1, cat2, cat3, EC, sep = "")) %>%
  filter(!duplicated(paste))


############################################################
## 14. Define KEGG pathways selected for visualization
############################################################

select_pathway <- c(
  "00260 Glycine, serine and threonine metabolism [PATH:ko00260]",
  "00270 Cysteine and methionine metabolism [PATH:ko00270]",
  "00280 Valine, leucine and isoleucine degradation [PATH:ko00280]",
  "00290 Valine, leucine and isoleucine biosynthesis [PATH:ko00290]",
  "00300 Lysine biosynthesis [PATH:ko00300]",
  "00310 Lysine degradation [PATH:ko00310]",
  "00340 Histidine metabolism [PATH:ko00340]",
  "00360 Phenylalanine metabolism [PATH:ko00360]",
  "00380 Tryptophan metabolism [PATH:ko00380]",
  "00400 Phenylalanine, tyrosine and tryptophan biosynthesis [PATH:ko00400]",
  "00670 One carbon pool by folate [PATH:ko00670]",
  "00730 Thiamine metabolism [PATH:ko00730]",
  "00740 Riboflavin metabolism [PATH:ko00740]",
  "00750 Vitamin B6 metabolism [PATH:ko00750]",
  "00760 Nicotinate and nicotinamide metabolism [PATH:ko00760]",
  "00770 Pantothenate and CoA biosynthesis [PATH:ko00770]",
  "00780 Biotin metabolism [PATH:ko00780]",
  "00790 Folate biosynthesis [PATH:ko00790]",
  "00052 Galactose metabolism [PATH:ko00052]"
)

KO_select <- KO %>%
  filter(cat3 %in% select_pathway)


############################################################
## 15. Convert EC abundance table to long format
############################################################

ec_abundance_melt.df <- melt(ec.df, "sample_id")

colnames(ec_abundance_melt.df) <- c(
  "sample_id",
  "feature",
  "abundance"
)


############################################################
## 16. Load trajectory cluster assignments
############################################################

load("~/Downloads/traj_cluster_results.Rdata")


############################################################
## 17. Assign samples to developmental age bins
############################################################

metadata_cluster <- metadata_cluster %>%
  mutate(
    bins = cut(
      Day,
      breaks = seq(100, 800, by = 100),
      include.lowest = TRUE,
      right = FALSE
    )
  )

cluster_bin.df <- metadata_cluster %>%
  dplyr::select(subject_id, sample_id, bins, Cluster, Days)


############################################################
## 18. Merge EC abundance with trajectory and KEGG information
############################################################

ec_abundance_melt.df <- ec_abundance_melt.df %>%
  inner_join(cluster_bin.df)

ec_abundance_melt.df <- ec_abundance_melt.df %>%
  left_join(ec_select)

ec_abundance_melt.df <- ec_abundance_melt.df %>%
  mutate(hit = 1)

ec_abundance_melt.df <- ec_abundance_melt.df %>%
  inner_join(KO_select)


############################################################
## 19. Define pathway subgroups
############################################################

esentail_biosynthesis <- c(
  "00290 Valine, leucine and isoleucine biosynthesis [PATH:ko00290]",
  "00300 Lysine biosynthesis [PATH:ko00300]",
  "00400 Phenylalanine, tyrosine and tryptophan biosynthesis [PATH:ko00400]"
)

metabolism <- setdiff(
  select_pathway[1:10],
  esentail_biosynthesis
)


############################################################
## 20. Assign ECs to top-level metabolic groups
############################################################

ec_abundance_melt.df <- ec_abundance_melt.df %>%
  mutate(
    EC_top = case_when(
      cat3 %in% c(
        "00290 Valine, leucine and isoleucine biosynthesis [PATH:ko00290]"
      ) ~ "BACC biosynthesis",
      
      cat3 %in% c(
        "00400 Phenylalanine, tyrosine and tryptophan biosynthesis [PATH:ko00400]"
      ) ~ "AAA biosynthesis",
      
      cat3 %in% select_pathway[11:18] ~ "Vitamin B metabolism",
      
      cat3 %in% select_pathway[19] ~ "Galactose metabolism"
    )
  )


############################################################
## 21. Scale EC abundance within each EC feature
############################################################

ec_abundance_melt.df <- ec_abundance_melt.df %>%
  group_by(feature) %>%
  mutate(scale = scale(abundance))


############################################################
## 22. Summarize EC abundance by metabolic group
############################################################

ec_abundance_melt_sum_top.df <- ec_abundance_melt.df %>%
  ungroup() %>%
  mutate(paste = paste(EC_number, bins, sample_id, EC_top)) %>%
  filter(!duplicated(paste)) %>%
  group_by(sample_id, subject_id, bins, Cluster, EC_top, Days) %>%
  summarise(
    sum_group = sum(abundance),
    hit = sum(hit),
    .groups = "drop"
  )


############################################################
## 23. Prepare plotting and statistical testing dataframe
############################################################

df <- ec_abundance_melt_sum_top.df %>%
  filter(!is.na(EC_top)) %>%
  mutate(
    bins = cut(Days, breaks = seq(100, 800, by = 100), right = FALSE),
    bin_start = as.numeric(
      gsub("\\[|\\)|\\,", "", sub(",.*", "", as.character(bins)))
    ),
    bin_mid = bin_start + 50,
    bin_end = bin_start + 100
  )


############################################################
## 24. Perform pairwise Wilcoxon tests
############################################################

pw_raw <- df %>%
  group_by(EC_top, bins) %>%
  filter(n_distinct(Cluster) > 1) %>%
  rstatix::pairwise_wilcox_test(
    sum_group ~ Cluster,
    p.adjust.method = "none"
  )


############################################################
## 25. Apply global FDR correction
############################################################

pw <- pw_raw %>%
  ungroup() %>%
  mutate(
    p.adj = p.adjust(p, method = "fdr"),
    p.adj.signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      p.adj < 0.10  ~ "+",
      p.adj < 0.25  ~ "-",
      TRUE          ~ "ns"
    )
  )


############################################################
## 26. Calculate Wilcoxon effect sizes
############################################################

eff <- df %>%
  group_by(EC_top, bins) %>%
  filter(n_distinct(Cluster) > 1) %>%
  wilcox_effsize(sum_group ~ Cluster)


############################################################
## 27. Combine p-values and effect sizes
############################################################

stat_df <- left_join(
  pw,
  eff,
  by = c("EC_top", "bins", "group1", "group2")
) %>%
  mutate(
    pair = paste0(group1, "-", group2),
    bin_start = as.numeric(
      gsub("\\[|\\)|\\,", "", sub(",.*", "", as.character(bins)))
    ),
    bin_mid = bin_start + 50,
    bin_end = bin_start + 100,
    sig_bin = p.adj < 0.25
  )


############################################################
## 28. Identify age bins to highlight in Panel A
############################################################

highlight_bins <- stat_df %>%
  filter(p.adj < 0.25) %>%
  group_by(EC_top, bins, bin_start, bin_end) %>%
  summarise(max_r = max(effsize), .groups = "drop") %>%
  group_by(EC_top) %>%
  slice_max(max_r, n = 3, with_ties = FALSE) %>%
  mutate(ymin = -Inf, ymax = Inf)


############################################################
## 29. Define significance color mapping for Panel B
############################################################

signif_colors <- c(
  "***" = "#a50026",
  "**"  = "#d73027",
  "*"   = "#fc8d59",
  "+"   = "#fee08b",
  "-"   = "#d9ef8b",
  "ns"  = "#f0f0f0"
)


############################################################
## 30. Set facet order for metabolic groups
############################################################

ec_top_levels <- c(
  "BACC biosynthesis",
  "Galactose metabolism",
  "Vitamin B metabolism",
  "AAA biosynthesis"
)

df$EC_top <- factor(df$EC_top, levels = ec_top_levels)
highlight_bins$EC_top <- factor(highlight_bins$EC_top, levels = ec_top_levels)
stat_df$EC_top <- factor(stat_df$EC_top, levels = ec_top_levels)


############################################################
## 31. Generate Panel A
############################################################

p1 <- ggplot(df, aes(x = Days, y = sum_group, color = Cluster)) +
  geom_rect(
    data = highlight_bins,
    aes(xmin = bin_start, xmax = bin_end, ymin = ymin, ymax = ymax),
    fill = "grey85",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE,
    linewidth = 1.2,
    alpha = 0.25
  ) +
  facet_wrap(~ EC_top, scales = "free_y", ncol = 4) +
  scale_color_manual(
    values = c(
      "1" = "#4c956c",
      "2" = "#457b9d",
      "3" = "#e26d5c"
    )
  ) +
  labs(
    x = "Days",
    y = "Sum Group Abundance"
  ) +
  theme_cowplot() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "bottom"
  )


############################################################
## 32. Generate Panel B
############################################################

p2 <- ggplot(stat_df, aes(x = bin_mid, y = effsize)) +
  geom_rect(
    data = stat_df %>% filter(sig_bin),
    aes(xmin = bin_mid - 50, xmax = bin_mid + 50,
        ymin = -Inf, ymax = Inf),
    fill = "grey90",
    alpha = 0.4,
    inherit.aes = FALSE
  ) +
  geom_line(
    aes(group = pair),
    color = "black",
    linewidth = 0.7
  ) +
  geom_point(
    aes(fill = p.adj.signif),
    shape = 21,
    color = "black",
    size = 4,
    stroke = 0.3
  ) +
  facet_wrap(pair ~ EC_top, ncol = 4) +
  scale_fill_manual(
    values = signif_colors,
    name = "FDR Significance",
    breaks = c("***", "**", "*", "+", "-", "ns"),
    labels = c(
      "*** (<0.001)",
      "** (<0.01)",
      "* (<0.05)",
      "+ (<0.10)",
      "- (<0.25)",
      "ns (≥0.25)"
    )
  ) +
  guides(color = "none") +
  ylim(-0.05, 0.25) +
  labs(
    x = "Days (bin midpoint)",
    y = "Effect size (r)"
  ) +
  theme_cowplot() +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom"
  )


############################################################
## 33. Export combined figure using cowplot
############################################################

pdf("~/Desktop/SupFig4_loess.pdf", 14, 9)

plot_grid(
  p1,
  p2,
  labels = c("A", "B"),
  ncol = 1,
  rel_heights = c(1, 1.1),
  align = "v",
  axis = "lr"
)

dev.off()

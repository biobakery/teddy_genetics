############################################################
## Supplementary Figure 11
## SNP-level fgsea enrichment curves for genetic PC3.
############################################################

setwd("~/ddong/TEDDY_project/")

###########################
## 0. Load packages
###########################
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(fgsea)
library(cowplot)

###########################
## 1. Load GO gene sets
###########################
go_file <- "./data/GO/hg19_geneticPCsnp20_GO.txt"

go_raw <- fread(go_file, header = FALSE, fill = TRUE)

gene_sets <- melt(go_raw, id.vars = "V1") %>%
  dplyr::select(V1, value) %>%
  dplyr::rename(Term = V1, Gene = value) %>%
  filter(Gene != "")

###########################
## 2. Load SNP-to-gene annotation
###########################
map <- fread("./output/map.txt") %>%
  dplyr::rename(snp = rs) %>%
  filter(!duplicated(snp))

###########################
## 3. Load genetic PC SNP loadings
###########################
loadings <- read.table(
  "./genetics.bed/filter/teddy.qc.excl_outliers.eigenvec.var",
  header = TRUE
) %>%
  dplyr::select(2, 5:9) %>%
  dplyr::rename(snp = VAR) %>%
  mutate(
    `PC3-PC2` = PC3 - PC2,
    `PC3-PC4` = PC3 - PC4,
    `PC3-PC5` = PC3 - PC5
  )

###########################
## 4. Convert GO gene sets
## from gene-level to SNP-level
###########################
map_select <- map %>%
  filter(snp %in% loadings$snp) %>%
  dplyr::select(snp, hgnc_symbol)

gene_sets_snp <- gene_sets %>%
  left_join(map_select, by = c("Gene" = "hgnc_symbol")) %>%
  filter(!is.na(snp)) %>%
  distinct(Term, snp)

grouped_list_snp <- split(gene_sets_snp$snp, gene_sets_snp$Term)

###########################
## 5. Reshape SNP loading matrix
###########################
ml <- loadings %>%
  pivot_longer(
    cols = -snp,
    names_to = "PC",
    values_to = "loading"
  )

###########################
## 6. Perform SNP-level fgsea
###########################
pcs_to_plot <- c("PC3")

fgseaRes <- list()
ranked_stats <- list()

for (pc_i in pcs_to_plot) {
  
  ## Rank SNPs by signed loading
  snps <- ml %>%
    filter(PC == pc_i) %>%
    select(snp, loading) %>%
    distinct() %>%
    arrange(desc(loading))
  
  FCsnplist <- snps$loading
  names(FCsnplist) <- snps$snp
  
  ranked_stats[[pc_i]] <- FCsnplist
  
  ## Run fgsea
  fgseaRes[[pc_i]] <- fgsea(
    pathways = grouped_list_snp,
    stats = FCsnplist,
    minSize = 15,
    maxSize = 500
  ) %>%
    arrange(padj) %>%
    mutate(PC = pc_i)
}

fgsea_all <- bind_rows(fgseaRes)

###########################
## 7. Select significant pathways
###########################
sig_pathways <- fgsea_all %>%
  filter(padj < 0.1) %>%
  arrange(PC, padj)

## Plot only the top six pathways
sig_pathways_top <- sig_pathways %>%
  group_by(PC) %>%
  slice_head(n = 6) %>%
  ungroup()

###########################
## 8. Draw fgsea enrichment curves
###########################
plot_list <- list()

for (i in seq_len(nrow(sig_pathways_top))) {
  
  pc_i <- sig_pathways_top$PC[i]
  pathway_i <- sig_pathways_top$pathway[i]
  
  plot_list[[i]] <- plotEnrichment(
    pathway = grouped_list_snp[[pathway_i]],
    stats = ranked_stats[[pc_i]]
  ) +
    labs(
      title = paste0(
        pc_i, ": ", pathway_i,
        "\nFDR = ", signif(sig_pathways_top$padj[i], 3),
        ", ES = ", signif(sig_pathways_top$ES[i], 3),
        ", NES = ", signif(sig_pathways_top$NES[i], 3)
      ),
      x = "SNPs ranked by genetic PC loading",
      y = "Running enrichment score"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 8),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8)
    )
}

###########################
## 9. Save supplementary figure
###########################
pdf(
  "./figures/supplement/SupFig11enrichment.pdf",
  width = 12,
  height = 10
)

plot_grid(
  plotlist = plot_list,
  ncol = 3
)

dev.off()
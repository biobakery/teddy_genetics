#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plink_gsea.R [-i <input> -m <map> -o <output>]

Options:
   -i input
   -m map
   -o output prefix

' -> doc

opts <- docopt(doc)

#

library(dplyr)
library(reshape2)

loadings_raw <- read.csv(opts$i, header=T, check.names=F, sep="\t")

loadings <- loadings_raw %>%
	select(2, 5:9) %>%
	rename(snp=1)

map <- data.table::fread(opts$m) %>%
	as.data.frame()

##############################
# gene set enrichment analysis

library(fgsea)

# gene set

genes1 <- loadings %>%
	reshape2::melt(id="snp", variable.name="PC", value.name="Value") %>%
	merge(., map, by="snp") %>%
	dplyr::select(name, snp)

genes2 <- genes1 %>%
	group_by(name) %>%
	group_nest()

genes3 <- unlist(genes2$data, recursive=FALSE)

names(genes3) <- genes2$name

genes4 <- as.list(genes3)

# run gsea for PC

gsea <- list()

for (i in colnames(loadings)[-1]) {
		
ranks_pc <- loadings %>%
	select(snp, i) %>%
	tibble::deframe()

ranks_pc <- sort(ranks_pc, decreasing=TRUE)

res <- as.data.frame(fgsea(genes4, ranks_pc)) %>%
	filter(padj < 0.05) %>%
	arrange(ES) %>%
	select(-leadingEdge) %>%
	mutate(PC = i)

gsea <- rbind(gsea, res)

}

write.table(gsea, paste0(opts$o, "_gsea_all.tsv"), sep="\t", row.names=F, quote=F)

topUp <- gsea %>% 
    filter(ES > 0) %>%
    mutate(ES_padj = ES/padj) %>%
    group_by(PC) %>%
    top_n(10, wt=ES_padj) %>%
    arrange(PC, desc(ES_padj)) %>%
    ungroup() %>%
    as.data.frame()

topDown <- gsea %>% 
    filter(ES < 0) %>%
    mutate(ES_padj = -ES/padj) %>%
    group_by(PC) %>%
    top_n(10, wt=ES_padj) %>%
    arrange(PC, desc(ES_padj)) %>%
    ungroup() %>%
    as.data.frame()

topBoth <- rbind(topUp, topDown) %>%
	select(8, 1) %>%
	arrange(PC, pathway)
	
write.csv(topBoth, paste0(opts$o, "_gsea_top.tsv"), sep="\t", row.names=F, quote=F)

#####
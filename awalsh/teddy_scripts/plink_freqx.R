opts$i <- "species_pvalues_toddlers.tsv"
opts$m <- "days_549-1095.metaphlan2.major.tsv"
opts$d <- "metadata.tsv"
opts$g <- "../../../plink_output/genetics.SNPs.tsv"
opts$a <- "toddler"


df <- fread("toddlers.frqx") %>%
	select(2, 5, 6, 7) %>%
	rename(HOM_A1=2, HET=3, HOM_A2=4) %>%
	reshape2::melt(id="SNP") %>%
	group_by(SNP) %>%
	mutate(freqx = 100 * value / sum(value) ) %>%
	filter(value == min(value)) %>%
	ungroup() %>%
	filter(freqx > 1) %>%
	select(SNP) %>%
	distinct()
#!/usr/bin/env Rscript

require(docopt)
'Usage:
   plot_results_snps.R [-i <input> -o <output>]

Options:
   -i model stats
   -p p-value adjustment method
   -m microbiome features
   -g genetic features
   -d metadata
   -o output prefix
   
' -> doc

opts <- docopt(doc)

###

library(tidyverse)
library(reshape2)
library(ggrepel)
library(lmerTest)
library(scales)
library(cowplot)

# flags: stats (fixed + interaction), microbiome, genetics, metadata, output (directory)

####################
# input for models #
####################

model_stats <- data.table::fread(opts$i, sep="\t", header=T)

####################
# p-value adjustment

head(model_stats)

levels(as.factor(model_stats$Type))

# fixed effect

fe <- data.table::fread("~/Desktop/teddy2/test/mphlan_infants_model_stats_genetics_fixed_effects.tsv")

sig_fe <- fe %>%
	top_n(P_adj, n=-25) %>%
	mutate(SNP = gsub("snps_", "", SNP)) %>%
	distinct() %>%
	mutate(Model = paste0("lmer(", MicroFeat, " ~ ", SNP, " * ", "Day + Country + ( 1 + Day | Subject ), data=data, na.action=na.omit)")) %>%
	as.data.frame()

# interaction effect

int <- data.table::fread("~/Desktop/teddy2/test/mphlan_infants_model_stats_genetics_interaction_effects.tsv")

sig_int <- int %>%
	filter(P_adj < 0.05) %>%
	top_n(P_adj, n=-25) %>%
	mutate(SNP = gsub("snps_", "", SNP)) %>%
	distinct() %>%
	mutate(Model = paste0("lmer(", MicroFeat, " ~ ", SNP, " * ", "Day + Country + ( 1 + Day | Subject ), data=data, na.action=na.omit)"))

# list significant features

sig_bugs <- c(sig_fe$MicroFeat, sig_int$MicroFeat) %>%
	as.data.frame() %>%
	distinct() %>%
	rename(MicroFeat = 1)

sig_snps <- c(sig_fe$SNP, sig_int$SNP) %>%
	as.data.frame() %>%
	distinct() %>%
	rename(SNP = 1)

####################
# input for models #
####################

# genetics

genetics <- read.csv("gen_micro_test.qc.genetics_updated.tsv", header=T, sep="\t", check.names=F) %>%
	mutate(RsID = gsub("-", "_", RsID)) %>%
	filter(RsID %in% levels(as.factor(sig_snps$SNP))) %>%
	column_to_rownames("RsID") %>%
	t(.) %>%
	as.data.frame() %>%
	setNames(paste0('genetics.', names(.))) %>%
	rownames_to_column("Subject")

# metaphlan2

mphlan <- read.csv("days_549-1095_metaphlan2_major.tsv", header=T, sep="\t", row.names=1) %>%
	t(.) %>%
	as.data.frame() %>%
	rownames_to_column("Species") %>%
	mutate(Species = gsub("microbiome\\.", "", Species)) %>%
	filter(Species %in% levels(as.factor(sig_bugs$MicroFeat))) %>%
	column_to_rownames("Species") %>%
	t(.) %>%
	as.data.frame() %>%
	setNames(paste0('microbiome.', names(.))) %>%
	rownames_to_column("sample_id")

# metadata

metadata <- read.csv("metadata.tsv", header=T, sep="\t")

df <- merge(genetics, metadata, by="Subject") %>%
	merge(., mphlan, by="sample_id")

# list

micro_list <- colnames(select(df, contains("microbiome."))) %>%
	gsub("microbiome\\.", "", .) %>%
	gsub("\\.", "_", .)

genet_list <- colnames(select(df, contains("genetics."))) %>%
	gsub(".*\\.", "", .)

# model input
	
data <- df %>%
	setNames(gsub("microbiome\\.|genetics\\.", "", names(.))) %>%
	setNames(gsub("\\.", "_", names(.)))

########################
# combinations of nucl #
########################

nucl_1 <- combn(c("A", "T", "C", "G"), 2, simplify=TRUE) %>%
	t(.) %>%
	as.data.frame() %>%
	mutate(combo = paste0(V1, V2)) %>%
	mutate(combo_rev = stringi::stri_reverse(combo)) %>%
	mutate(combination = paste0(combo, "/", combo_rev)) %>%
	select(3:5) %>%
	melt(id="combination") %>%
	select(3, 1) %>%
	rename(Allele=1)
	
nucl_2 <- expand.grid(c("A", "T", "C", "G"), c("A", "T", "C", "G")) %>%
	mutate(combination = case_when(Var1 == Var2 ~ paste0(Var1, Var2))) %>%
	drop_na() %>%
	mutate(Allele = combination) %>%
	select(4, 3)

nucl <- rbind(nucl_1, nucl_2)

#################
# fixed effects #
#################

sig_micro <- fe %>%
	filter(!grepl(":", Type)) %>%
	filter(P_adj < 0.05) %>%
	select(MicroFeat) %>%
	distinct()

##################
# Manhattan plot #
##################

DF <- fe %>%
	filter(!grepl(":", Type)) %>%
	merge(., sig_micro, by="MicroFeat") %>%
	as.data.frame()

# positions

rsid_map <- read.csv("gen_micro_test.qc.map",
		header=F,
		sep="\t") %>%
	rename(Chr=1, RsID=2, Morgans=3, BP=4)
	
rsid_list <- read.csv("gen_micro_test.qc.RsID_list.tsv",
	header=T,
	sep="\t")

map <- merge(rsid_map, rsid_list, by="RsID") %>%
	mutate(SNP = gsub("-", "_", RsID)) %>%
	mutate(SNP = paste0("snps_", SNP)) %>%
	arrange(Chr, BP) %>%
	mutate(BP = cumsum(as.numeric(BP)))

# merge

DF2 <- merge(DF, map, by="SNP")

# set x-axis labels to the centre of each chromosome

axis.set <- map %>% 
  group_by(Chr) %>% 
  summarize(center = (max(BP) + min(BP)) / 2)

# set the y-axis limits

ylim <- abs(floor(log10(min(map$P)))) + 2 

# the significance threshold

adjusted <- DF2 %>%
	filter(P_adj < 0.05)
sig <- max(adjusted$P)

# the number of chromosomes

nCHR <- nlevels(as.factor(map$Chr))

# annotation file

ann <- DF2 %>%
	filter(P_adj < 0.05) %>%
	group_by(MicroFeat) %>%
	top_n(P, n=-5)

# generate the plot

manhplot <- ggplot(DF2, aes(x = BP, y = -log10(P), 
                                 color = as.factor(Chr), size = -log10(P))) +
  facet_wrap(~MicroFeat, scales="fixed", ncol=5) +
  geom_point(alpha = 0.75) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$Chr, 
    breaks = axis.set$center,
    guide = guide_axis(check.overlap = TRUE)) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  scale_color_manual(values = rep(c("orangered", "midnightblue"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log10(p)") + 
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, size = 5, vjust = 0.5),
    axis.text.y = element_text(size=10),
    axis.title = element_text(size=12.5, face="bold"),
    strip.text = element_text(size=12.5, face="bold")) +
  labs(x="Chromosome", y=expression(bold(-log[10]*(p[unadjusted])))) +
  geom_label_repel(data=ann, aes(x = BP, y = -log10(P), label=RsID), colour="black")

mh <- ceiling( nrow(sig_micro) / 5 ) * 4

if (nrow( sig_micro ) / 5 > 5){ mw = 25 } else { mw = nrow(sig_micro) / 5 * 25 }

ggsave(manhplot, file="infants_manhattan_plot.png", dpi=300, width=mw, height=mh)

###############################
# fixed effects estimate plot #
###############################

map <- data.table::fread("~/Desktop/teddy2/test/gen_micro_test.qc.rsid_to_gene.map.txt") %>%
	rename(SNP=snp)

models_fe <- levels(as.factor(sig_fe$Model))

combined_fe <- list()

for (i in models_fe) {

model <- eval(parse(text=i))

slope <- as.data.frame(coef(model)$Subject) %>%
	tibble::rownames_to_column("Subject") %>%
	mutate(Model = i)

tidy_mod <- broom.mixed::tidy(model) %>%
	as.data.frame() %>%
	mutate(Allele = str_sub(term, start=-2)) %>%
	mutate(SNP = str_sub(term, end=-3)) %>%
	filter(SNP %in% levels(as.factor(sig_snps$SNP))) %>%
	select(10, 9, 4:8) %>%
	mutate(MicroFeat = i) %>%
	mutate(MicroFeat = gsub("lmer\\(", "", MicroFeat)) %>%
	mutate(MicroFeat = gsub(" .*", "", MicroFeat))

combined_fe <- rbind(combined_fe, tidy_mod)	

}

df <- combined_fe %>%
	merge(., rbind(nucl_1, nucl_2), by="Allele") %>%
	merge(., map, by="SNP") %>%
	mutate(SNP_Gene = paste0(SNP, " (", name, ")"))

fe_estimate <- ggplot(data=df, aes(x=SNP_Gene, y=estimate)) +
	facet_grid(MicroFeat ~ ., scales="free", space="free", switch="both") +
	geom_hline(yintercept=0, linetype="dashed", colour="gray") +
	geom_pointrange(aes(ymin=(estimate-std.error), ymax=(estimate+std.error), colour=combination)) +
	theme_bw() +
	labs(title="Estimated fixed effects (top 25)",
		colour="Allele",
		x="SNP (closest gene)",
		y="Estimate") +
	scale_colour_brewer(palette="Dark2") +
	scale_fill_brewer(palette="Dark2") +
	scale_x_discrete(position="top") +
	coord_flip() +
	theme(
		panel.grid=element_blank(),
		plot.title = element_text(size=12.5, face="bold", hjust=0.5),
		axis.title = element_text(size=12.5, face="bold"),
		axis.text = element_text(size=8.75),
		legend.title = element_text(size=12.5, face="bold", vjust=0.75),
		legend.text = element_text(size=12.5),
		legend.position = "top",
		strip.text.y.left = element_text(size=8.75, angle=0, face="bold"),
		strip.placement = "inside")

ggsave(fe_estimate, file="infants_top_25_fixed_effects.png", dpi=300, height=7.5, width=11.25)

################
# interactions #
################

if (nrow(sig_int) > 25){
	
	sig_int_top <- sig_int %>%
		top_n(P_adj, n=-25)
	
	models <- levels(as.factor(sig_int_top$Model))
	
} else {
	
	models <- levels(as.factor(sig_int$Model))
	
}

#########################
# interaction line plot #
#########################

combined_int <- list()

for (i in models) {

model <- eval(parse(text=i))

slope <- as.data.frame(coef(model)$Subject) %>%
	tibble::rownames_to_column("Subject") %>%
	mutate(Model = i)

aug <- broom.mixed::augment(model) %>%
	select(-matches(".rownames"), everything()) %>%
	rename(RelAb=1, Allele=2) %>%
	mutate(Model = i) %>%
	mutate(MicroFeat = gsub("lmer\\(| .*", "", Model)) %>%
	mutate(MicroFeat = gsub("_noname_unclassified", "_?_?", MicroFeat)) %>%
	mutate(MicroFeat = gsub("_unclassified", "_?", MicroFeat)) %>%
	mutate(SNP = gsub(".* \\~ ", "", Model)) %>%
	mutate(SNP = gsub(" .*", "", SNP)) %>%
	select(Subject, MicroFeat, SNP, Allele, Day, .fitted)

combined_int <- rbind(combined_int, aug)	

}

combined_int2 <- merge(nucl, combined_int, by="Allele")

ann_int <- sig_int_top %>%
	 mutate(P_adj = paste0("q=", scientific(P_adj))) %>%
	 select(MicroFeat, SNP, P_adj) %>%
	 arrange(MicroFeat, P_adj)

interaction <- ggplot(combined_int2, aes(x=Day, y=.fitted)) +
	facet_wrap(~ MicroFeat + SNP, scales="free", ncol=5) +
	geom_smooth(aes(colour=combination, fill=combination), method="lm") +
	theme_bw() +
	theme(
		plot.title = element_text(size=12.5, face="bold", hjust=0.5),
		axis.title = element_text(size=12.5, face="bold"),
		axis.text = element_text(size=8.75),
		legend.title = element_text(size=12.5, face="bold", vjust=0.75),
		legend.text = element_text(size=12.5),
		strip.text = element_text(size=8.75, face="bold"),
		legend.position = "top") +
	labs(title="Interaction effect(s)", 
		colour="Allele",
		fill="Allele") +
	scale_colour_brewer(palette="Dark2") +
	scale_fill_brewer(palette="Dark2") +
	geom_text(data=ann_int, aes(x=Inf, y=Inf, label=P_adj), 
		hjust=1.5,
		vjust=1.5)

ggsave(interaction, file="_top_25_interaction_effect_line.png", dpi=300, height=15, width=15)

#############################
# interaction slope boxplot #
#############################

library(purrr)
library(broom)

slope_df <- combined_int %>%
  split(., list(.$Subject, .$MicroFeat, .$SNP), drop = TRUE) %>%
  map(~lm(.fitted ~ Day, data = .x)) %>% 
  map_df(tidy, .id = "Subject.MicroFeat.SNP") %>%
  filter(term == "Day") %>%
  mutate(SNP = gsub(".*\\.", "", Subject.MicroFeat.SNP)) %>%
  mutate(Subject.MicroFeat = str_sub(sub(sprintf("^((?:[^.]*.){2}).*", 2), "\\1", Subject.MicroFeat.SNP), end=-2)) %>%
  mutate(Subject = gsub("\\..*", "", Subject.MicroFeat)) %>%
  mutate(MicroFeat = gsub(".*\\.", "", Subject.MicroFeat)) %>%
  select(SNP, MicroFeat, Subject, estimate) %>%
  as.data.frame()

slope_map <- combined_int %>%
	select(SNP, MicroFeat, Subject, Allele)

DF <- merge(slope_df, slope_map, by=c("SNP", "MicroFeat", "Subject")) %>%
	select(Subject, SNP, Allele, MicroFeat, estimate) %>%
	drop_na() %>%
	distinct() %>%
	merge(., nucl, by="Allele")

slope_box <- ggplot(data=DF, aes(x=Allele, y=estimate)) +
	facet_wrap(MicroFeat~SNP, scales="free", ncol=5) +
	geom_jitter(aes(colour=combination, fill=combination), pch=21, alpha=0.25) +
	geom_boxplot(aes(colour=combination), fill=NA, outlier.shape=NA) +
	theme_bw() +
	theme(
		panel.grid=element_blank(),
		plot.title = element_text(size=12.5, face="bold", hjust=0.5),
		axis.title = element_text(size=12.5, face="bold"),
		axis.text = element_text(size=8.75),
		legend.title = element_text(size=12.5, face="bold", vjust=0.75),
		legend.text = element_text(size=12.5),
		strip.text = element_text(size=8.75, face="bold"),
		legend.position = "top"
		) +
	labs(title="Slope per subject", 
		colour="Allele", fill="Allele",
		y="slope") +
	scale_colour_brewer(palette="Dark2") +
	scale_fill_brewer(palette="Dark2")

ggsave(slope_box, file="_top_25_interaction_effects_box.png", dpi=300, height=15, width=15)

##################################
# comnbine the interaction plots #
##################################

fig_ab <- plot_grid(interaction, slope_box, 
	nrow=1, align="h", scale=0.9,
	labels=c("a", "b"), label_x=0.25, label_size = 17.5, label_colour="firebrick")

ggsave(fig_ab, file="_top_25_interaction_effects_cowplot.png", dpi=300, height=15, width=30)

###
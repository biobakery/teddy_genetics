#!/usr/bin/env Rscript

require(docopt)

'Usage:
   run_models.R [-i <input> -c <cores> -o <output>]

Options:
   -i input
   -c number of cpus [default: 1]
   -o model output

' -> doc

opts <- docopt(doc)

###########################
# make a list of the models

library(dplyr)

df <- read.csv(opts$i, sep="\t", header=T)

micro_list <- colnames(select(df, contains("microbiome."))) %>%
	gsub("microbiome\\.", "", .) %>%
	gsub("\\.", "_", .)

genet_list <- colnames(select(df, contains("genetics."))) %>%
	gsub(".*\\.", "", .)

data <- df %>%
	setNames(gsub("microbiome\\.|genetics\\.", "", names(.))) %>%
	setNames(gsub("\\.", "_", names(.)))

#

formulae <- list()

for (i in micro_list) {
	
	for (j in genet_list) {
		
		formula <- data.frame(formula = paste0("lmer(", i, " ~ ", j, " * ", "Days + Country + ( 1 + Days | Subject ), data=data)"))
		
		formulae <- rbind(formulae, formula)
	}
}

formulae <- levels(as.factor(formulae$formula))

###################################
# loop through the list in parallel

library(doParallel)

registerDoParallel(cores=opts$c)

system.time(model_anova_total <- foreach(i=formulae, 
			      .combine=rbind, 
			      .packages=c("lmerTest", "tidyverse"), 
			      .errorhandling = 'remove'
			      )	%dopar% {

	model <- eval(parse( text= i ))

	model_coef <- as.data.frame(coef(summary(model))[ , "Estimate"]) %>%
		rename(Coefficient=1) %>%
		tibble::rownames_to_column("Predictor") %>%
		filter(!grepl("Intercept|Country", Predictor))
	
	model_anova <- anova(model) %>%
		tibble::rownames_to_column("Predictor") %>%
		select(Predictor, `Pr(>F)`) %>%
		rename(P = `Pr(>F)`) %>%
		mutate(Type = gsub(".*\\:", "PC\\:", Predictor)) %>%
		mutate(Predictor = factor(Predictor)) %>%
		mutate(Model = i) %>%
		mutate(Coefficient = c(.subset2(model_coef, 2)[1], .subset2(model_coef, 2)[2], NA, .subset2(model_coef, 2)[3]))

	})

write.table(model_anova_total, opts$o, sep="\t", quote=F, row.names=F)

#####
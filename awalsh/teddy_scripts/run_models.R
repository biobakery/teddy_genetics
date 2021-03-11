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
	gsub("genetics\\.", "genetics_", .) %>%
	gsub("-", "_", .) %>%
	gsub("\\.", "_", .)

data <- df %>%
	setNames(gsub("microbiome\\.", "", names(.))) %>%
	setNames(gsub("genetics\\.", "genetics_", names(.))) %>%
	setNames(gsub("\\.|-", "_", names(.)))

#

formulae <- list()

for (i in micro_list) {
	
	for (j in genet_list) {
		
		formula <- data.frame(formula = paste0("lmer(", i, " ~ ", j, " * ", "Days + Country + ( 1 + Days | Subject ), data=data, na.action=na.omit)"))
		
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
			      .packages=c("lmerTest", "tidyverse", "performance", "DHARMa"),
			      .errorhandling = 'remove'
			      )	%dopar% {

	model <- eval(parse( text= i ))

	model_coef <- as.data.frame(coef(summary(model))[ , "Estimate"]) %>%
		rename(Coefficient=1) %>%
		rownames_to_column("Predictor") %>%
		filter(!grepl("Intercept|Country", Predictor))
	
	R2_conditional <- (data.frame(r2(model)$R2_conditional))[1,1]
	
	simulationOutput <- simulateResiduals(fittedModel = model, refit = F)
	p_uni <- testUniformity(simulationOutput, plot=F)$p.value
	p_out <- testOutliers(simulationOutput, plot=F, type="binomial")$p.value
	p_dis <- testDispersion(simulationOutput, plot=F)$p.value
	
	model_anova <- anova(model) %>%
		rownames_to_column("Predictor") %>%
		select(Predictor, `Pr(>F)`) %>%
		rename(P = `Pr(>F)`) %>%
		mutate(Model = i) %>%
		mutate(Coefficient = c(.subset2(model_coef, 2)[1], .subset2(model_coef, 2)[2], NA, .subset2(model_coef, 2)[3])) %>%
		mutate(R2_conditional=R2_conditional, Uniformity=p_uni, Outliers=p_out, Dispersion=p_dis) %>%
		as.data.frame()

	})

head(model_anova_total)

write.table(model_anova_total, opts$o, sep="\t", quote=F, row.names=F)

#####
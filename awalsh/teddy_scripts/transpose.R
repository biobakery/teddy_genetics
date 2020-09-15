#!/usr/bin/env Rscript

require(docopt, quietly=TRUE)
'Usage:
   transpose.R [-i <input> -o <output>]

Options:
   -i input (.tsv file)
   -o output (.tsv file)

' -> doc

opts <- docopt(doc)

###

suppressPackageStartupMessages(library(tidyverse))

df <- read.csv(opts$i, sep="\t", header=T, row.names=1, check.names=F) %>%
   t(.) %>%
   as.data.frame() %>%
   rownames_to_column("Column1")

write.table(df, opts$o, sep="\t", quote=F, row.names=F)

###
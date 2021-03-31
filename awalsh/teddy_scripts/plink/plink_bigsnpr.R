#!/usr/bin/env Rscript

require(docopt)

'Usage:
   plink_bigsnpr.R [-i <input> -o <output>]

Options:
   -i input (.bed)
   -o output (.bed)
    
' -> doc

opts <- docopt(doc)

###

library(bigsnpr)

tmp <- snp_readBed(opts$i, backingfile = tempfile())

bigSNP_raw <- snp_attach(tmp)

auto_SVD_out <- bed_autoSVD(bed(opts$i))

bigSNP_filtered <- snp_subset(bigSNP_raw, ind.col = attr(auto_SVD_out, "subset"))

snp_writeBed(snp_attach(bigSNP_filtered), opts$o)

###

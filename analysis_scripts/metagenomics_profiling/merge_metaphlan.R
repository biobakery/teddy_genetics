###############################################################################
## Merge per-sample MetaPhlAn2 taxonomic profiles into one species x sample
## matrix. Generalises the repo's two-file demo to all samples: full_join by
## bug name, fill NA -> 0, sort rows (bugs) and columns (samples).
##
## Usage:  Rscript merge_metaphlan.R <profiles_dir> <output.tsv>
##   <profiles_dir> : folder of per-sample *_taxonomic_profile.tsv
##   <output.tsv>   : merged matrix (default merged_metaphlan_table.tsv)
###############################################################################
suppressWarnings(suppressMessages({
  library(dplyr)
  library(purrr)
}))

args    <- commandArgs(trailingOnly = TRUE)
in_dir  <- ifelse(length(args) >= 1, args[1], "metaphlan_outputs")
out_tsv <- ifelse(length(args) >= 2, args[2], "merged_metaphlan_table.tsv")

files <- list.files(in_dir, pattern = "\\.tsv$", full.names = TRUE)
stopifnot(length(files) > 0)

## Read one profile as (bugname, <sample>); sample name = file basename
read_one <- function(f) {
  s <- tools::file_path_sans_ext(basename(f))
  read.table(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
             colClasses = "character", comment.char = "") %>%
    setNames(c("bugname", s))
}

merged <- reduce(map(files, read_one), full_join, by = "bugname")
merged[is.na(merged)] <- 0

## sort rows by bug, columns by sample
sample_cols <- sort(setdiff(names(merged), "bugname"))
merged <- merged[order(merged$bugname), c("bugname", sample_cols)]

write.table(merged, out_tsv, sep = "\t",
            quote = FALSE, row.names = FALSE, na = "0")

message("Wrote: ", out_tsv, "  (", nrow(merged), " taxa x ",
        length(sample_cols), " samples)")

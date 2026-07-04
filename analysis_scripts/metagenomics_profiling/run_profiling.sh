#!/bin/bash
###############################################################################
## Run MetaPhlAn2 + HUMAnN2 on paired-end shotgun stool reads (bioBakery wmgx).
##
## One command runs BOTH: HUMAnN2 does the functional profiling and internally
## calls MetaPhlAn2 for the taxonomic prescreen. Read QC + human-read removal
## already happened upstream, so quality control is bypassed here.
##
## Usage:  run_profiling.sh <input_dir> <output_dir>
##   <input_dir>  : decompressed paired-end fastq named <sample>_1.fastq / _2.fastq
##   <output_dir> : where per-sample HUMAnN2 / MetaPhlAn2 outputs are written
###############################################################################
set -euo pipefail

INPUT_DIR=${1:?usage: run_profiling.sh <input_dir> <output_dir>}
OUTPUT_DIR=${2:?usage: run_profiling.sh <input_dir> <output_dir>}

## If inputs are still .bz2, decompress first:
##   bzip2 --decompress "$INPUT_DIR"/*.bz2

biobakery_workflows wmgx \
  --input  "$INPUT_DIR" \
  --output "$OUTPUT_DIR" \
  --pair-identifier _1 \
  --bypass-quality-control \
  --bypass-strain-profiling \
  --threads 8 \
  --input-extension fastq \
  --remove-intermediate-output
## On a SLURM cluster add:  --grid slurm --grid-jobs 2

###############################################################################
## Per-sample outputs produced:
##   HUMAnN2:    <sample>_pathabundance.tsv  <sample>_pathcoverage.tsv
##               <sample>_genefamilies.tsv   <sample>_ecs.tsv  (regrouped EC)
##   MetaPhlAn2: <sample>_taxonomic_profile.tsv
## -> merge with merge_metaphlan.R (species) and merge_humann2.sh (function).
###############################################################################

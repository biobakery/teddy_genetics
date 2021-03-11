#!/usr/bin/env python

from anadama2 import Workflow
import os
import subprocess
import re

class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end

workflow = Workflow(version="0.1", description="A workflow to run PC models for TEDDY2")

workflow.add_argument("subset-microbiome", desc="Should the samples be subset into two age groups?", choices=["yes", "no"], default="yes")
workflow.add_argument("min-a", desc="Minimum age (days) for the first subset", type=int, default=0)
workflow.add_argument("max-a", desc="Maximum age (days) for the first subset", type=int, default=548)
workflow.add_argument("min-b", desc="Minimum age (days) for the second subset", type=int, default=549)
workflow.add_argument("max-b", desc="Maximum age (days) for the second subset", type=int, default=1095)
workflow.add_argument("abundance", desc="Abundance threshold for microbiome features", type=float, default=1)
workflow.add_argument("prevalence", desc="Prevalence threshold for microbiome features", type=float, choices=[Range(0.0, 1.0)], default=0.25)
workflow.add_argument("transformation", desc="Transformation method for microbiome measurements", choices=["abs", "log", "log10", "sqrt"], required=True)
workflow.add_argument("genetics", desc="Genetics data", required=True)
workflow.add_argument("cores", desc="Number of processors to use for modeling", type=int, default=1)
workflow.add_argument("padj", desc="Method for p-value adjustment", choices=["bonferroni", "fdr"], default="bonferroni")
workflow.add_argument("group-by", desc="Group by Type or Predictor", choices=["Type", "Predictor"], default="Type")
workflow.add_argument("humann-dir", desc="HUMAnN2 directory", required=True)
workflow.add_argument("mpa", desc="MetaPhlAn2 table", required=True)

args = workflow.parse_args()

# set species targets

directory_species = "./PC_models/Species" 
if not os.path.exists(directory_species):
	os.makedirs(directory_species)

if args.subset_microbiome == "yes":
	prefix_species = [directory_species + "/days_" + str(args.min_a) + "-" + str(args.max_a) + ".metaphlan2.tsv", directory_species + "/days_" + str(args.min_b) + "-" + str(args.max_b) + ".metaphlan2.tsv"]
else:
	prefix_species = [directory_species + "/days_" + str(args.min_a) + "-" + str(args.max_b) + ".metaphlan2.tsv"]

targets_species = []
for i in prefix_species:
   	targets_species.append(i.replace(".tsv", ".major.model_output.tsv"))

depends_species = ["metaphlan2.tsv", "metadata.tsv"]

workflow.add_task(
	"anadama2_models_PCs.py " \
	"--microbiome metaphlan2.tsv " \
	"--microbiome-feature Species " \
	"--metadata metadata.tsv " \
	"--subset-microbiome [subset_microbiome] " \
	"--min-a [min_a] " \
	"--max-a [max_a] " \
	"--min-b [min_b] " \
	"--max-b [max_b] " \
	"--abundance [abundance] " \
	"--prevalence [prevalence] " \
	"--transformation [transformation] " \
	"--genetics [genetics] " \
	"--cores [cores] " \
	"--padj [padj] " \
	"--group-by [group_by] " \
	"-o [output]",
	depends=depends_species,
	targets=targets_species,
	name="Species x PC models",
	subset_microbiome=args.subset_microbiome,
	min_a=args.min_a,
	max_a=args.max_a,
	min_b=args.min_b,
	max_b=args.max_b,
	abundance=args.abundance,
	prevalence=args.prevalence,
	transformation=args.transformation,
	genetics=args.genetics,
	cores=args.cores,
	padj=args.padj,
	group_by=args.group_by,
	output=args.output
	)

# set EC targets

directory_ec = "./PC_models/EC" 
if not os.path.exists(directory_ec):
	os.makedirs(directory_ec)

if args.subset_microbiome == "yes":
	prefix_ec = [directory_ec + "/days_" + str(args.min_a) + "-" + str(args.max_a) + ".humann2_ec.tsv", directory_ec + "/days_" + str(args.min_b) + "-" + str(args.max_b) + ".humann2_ec.tsv"]
else:
	prefix_ec = [directory_ec + "/days_" + str(args.min_a) + "-" + str(args.max_b) + ".humann2_ec.tsv"]

targets_ec = []
for i in prefix_ec:
   	targets_ec.append(i.replace(".tsv", ".major.model_output.tsv"))

depends_ec = ["humann2_ec.tsv", "metadata.tsv", "metaphlan2.tsv", args.humann_dir + "map_level4ec_name.txt", args.humann_dir + "humann2_ec_stratified.tsv"]

ec_names = args.humann_dir + "map_level4ec_name.txt"
ec_strat = args.humann_dir + "humann2_ec_stratified.tsv"

workflow.add_task(
	"anadama2_models_PCs.py " \
	"--microbiome humann2_ec.tsv " \
	"--microbiome-feature EC " \
	"--metadata metadata.tsv " \
	"--subset-microbiome [subset_microbiome] " \
	"--min-a [min_a] " \
	"--max-a [max_a] " \
	"--min-b [min_b] " \
	"--max-b [max_b] " \
	"--abundance [abundance] " \
	"--prevalence [prevalence] " \
	"--transformation [transformation] " \
	"--genetics [genetics] " \
	"--cores [cores] " \
	"--padj [padj] " \
	"--group-by [group_by] " \
	"--path-names [path_names] " \
	"--path-strat [path_strat] " \
	"--mpa metaphlan2.tsv " \
	"-o [output]",
	depends=depends_species,
	targets=targets_species,
	name="EC x PC models",
	subset_microbiome=args.subset_microbiome,
	min_a=args.min_a,
	max_a=args.max_a,
	min_b=args.min_b,
	max_b=args.max_b,
	abundance=args.abundance,
	prevalence=args.prevalence,
	transformation=args.transformation,
	genetics=args.genetics,
	cores=args.cores,
	padj=args.padj,
	group_by=args.group_by,
	path_names=ec_names,
	path_strat=ec_strat,
	output=args.output
	)

# set Pfam targets

directory_pfam = "./PC_models/Pfam" 
if not os.path.exists(directory_pfam):
	os.makedirs(directory_pfam)

depends_pfam = ["humann2_pfam.tsv", "metadata.tsv", "metaphlan2.tsv"]

if args.subset_microbiome == "yes":
	prefix_pfam = [directory_pfam + "/days_" + str(args.min_a) + "-" + str(args.max_a) + ".humann2_pfam.tsv", directory_pfam + "/days_" + str(args.min_b) + "-" + str(args.max_b) + ".humann2_pfam.tsv"]
else:
	prefix_pfam = [directory_pfam + "/days_" + str(args.min_a) + "-" + str(args.max_b) + ".humann2_pfam.tsv"]

targets_pfam = []
for i in prefix_pfam:
   	targets_pfam.append(i.replace(".tsv", ".major.model_output.tsv"))

depends_pfam = ["humann2_pfam.tsv", "metadata.tsv", "metaphlan2.tsv", args.humann_dir + "map_pfam_name.txt", args.humann_dir + "humann2_pfam_stratified.tsv"]

pfam_names = args.humann_dir + "map_pfam_name.txt"
pfam_strat = args.humann_dir + "humann2_pfam_stratified.tsv"

workflow.add_task(
	"anadama2_models_PCs.py " \
	"--microbiome humann2_pfam.tsv " \
	"--microbiome-feature Pfam " \
	"--metadata metadata.tsv " \
	"--subset-microbiome [subset_microbiome] " \
	"--min-a [min_a] " \
	"--max-a [max_a] " \
	"--min-b [min_b] " \
	"--max-b [max_b] " \
	"--abundance [abundance] " \
	"--prevalence [prevalence] " \
	"--transformation [transformation] " \
	"--genetics [genetics] " \
	"--cores [cores] " \
	"--padj [padj] " \
	"--group-by [group_by] " \
	"--path-names [path_names] " \
	"--path-strat [path_strat] " \
	"--mpa metaphlan2.tsv " \
	"-o [output]",
	depends=depends_pfam,
	targets=targets_pfam,
	name="Pfam x PC models",
	subset_microbiome=args.subset_microbiome,
	min_a=args.min_a,
	max_a=args.max_a,
	min_b=args.min_b,
	max_b=args.max_b,
	abundance=args.abundance,
	prevalence=args.prevalence,
	transformation=args.transformation,
	genetics=args.genetics,
	cores=args.cores,
	padj=args.padj,
	group_by=args.group_by,
	path_names=pfam_names,
	path_strat=pfam_strat,
	output=args.output
	)

# set stability targets

directory_bc = "./PC_models/Bray-Curtis" 
if not os.path.exists(directory_bc):
	os.makedirs(directory_bc)

if args.subset_microbiome == "yes":
	prefix_bc = [directory_bc + "/days_" + str(args.min_a) + "-" + str(args.max_a) + ".bc_dists.tsv", directory_bc + "/days_" + str(args.min_b) + "-" + str(args.max_b) + ".bc_dists.tsv"]
else:
	prefix_bc = [directory_bc + "/days_" + str(args.min_a) + "-" + str(args.max_b) + ".bc_dists.tsv"]

targets_bc = []
for i in prefix_bc:
   	targets_bc.append(i.replace(".tsv", ".major.model_output.tsv"))

depends_bc = ["bc_dists.tsv", "metadata.tsv"]

workflow.add_task(
	"anadama2_models_PCs.py " \
	"--microbiome bc_dists.tsv " \
	"--microbiome-feature Bray-Curtis " \
	"--metadata metadata.tsv " \
	"--subset-microbiome [subset_microbiome] " \
	"--min-a [min_a] " \
	"--max-a [max_a] " \
	"--min-b [min_b] " \
	"--max-b [max_b] " \
	"--abundance [abundance] " \
	"--prevalence [prevalence] " \
	"--transformation [transformation] " \
	"--genetics [genetics] " \
	"--cores [cores] " \
	"--padj [padj] " \
	"--group-by [group_by] " \
	"-o [output]",
	depends=depends_species,
	targets=targets_species,
	name="Bray-Curtis x PC models",
	subset_microbiome=args.subset_microbiome,
	min_a=args.min_a,
	max_a=args.max_a,
	min_b=args.min_b,
	max_b=args.max_b,
	abundance=args.abundance,
	prevalence=args.prevalence,
	transformation="abs",
	genetics=args.genetics,
	cores=args.cores,
	padj=args.padj,
	group_by=args.group_by,
	output=args.output
	)

# run the workflow

title = "# Microbiome x PC models #"

print("#" * (len(title) - 1) + "\n" + title + "\n" + "#" * (len(title) - 1))
workflow.go()

for i in range(25):
	print(args.output)

#

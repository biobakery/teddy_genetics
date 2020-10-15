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

workflow = Workflow(version="0.1", description="A workflow to run linear mixed effect models for TEDDY2")

workflow.add_argument("microbiome", desc="Microbiome data", required=True)
workflow.add_argument("pathway", desc="Add this flag if the microbiome features are pathways")
workflow.add_argument("microbiome-feature", desc="Type of microbiome feature", choices=["Bray-Curtis", "EC", "Pfam", "Species"], required=True)
workflow.add_argument("metadata", desc="Metadata for samples", required=True)
workflow.add_argument("subset-microbiome", desc="Should the samples be subset into two age groups?", choices=["yes", "no"], default="yes")
workflow.add_argument("min-a", desc="Minimum age (days) for the first subset", type=int, default=0)
workflow.add_argument("max-a", desc="Maximum age (days) for the first subset", type=int, default=548)
workflow.add_argument("min-b", desc="Minimum age (days) for the second subset", type=int, default=549)
workflow.add_argument("max-b", desc="Maximum age (days) for the second subset", type=int, default=1095)
workflow.add_argument("abundance", desc="Abundance threshold for microbiome features", type=float, default=0)
workflow.add_argument("prevalence", desc="Prevalence threshold for microbiome features", type=float, choices=[Range(0.0, 1.0)], default=0.1)
workflow.add_argument("transformation", desc="Transformation method for microbiome measurements", choices=["abs", "log", "log10", "sqrt"], required=True)
workflow.add_argument("genetics", desc="Genetics data", required=True)
workflow.add_argument("genetics-feature", desc="Genetics data", choices=["PC", "SNP"], required=True)
workflow.add_argument("subset-genetics", desc="Should the genetics data be subset into chunks?", choices=["yes", "no"], default="no")
workflow.add_argument("size", desc="Size of the chunks into which the genetics data are subsetted", type=int, default=250)
workflow.add_argument("threads", desc="Number of threads to use when splitting the genetics data into subsets", type=int, default=1)
workflow.add_argument("system", desc="Where are the models run from?", choices=["server", "cluster"], default="server")
workflow.add_argument("cores", desc="Number of processors to use for modeling", type=int, default=1)
workflow.add_argument("time", desc="Limit on runtime (minutes) off jobs on cluster", type=int, default=720)
workflow.add_argument("mem", desc="Amount of memory (MB) to allocate to processors on cluster", type=int, default=1000)
workflow.add_argument("partition", desc="Partition to which jobs are submitted on cluster", default="serial_requeue")
workflow.add_argument("padj", desc="Method for p-value adjustment", choices=["bonferroni", "fdr"], default="bonferroni")
workflow.add_argument("group-by", desc="Group by Type or Predictor", choices=["Type", "Predictor"], default="Type")
workflow.add_argument("path-names", desc="Path to directory containing HUMAnN2 pathway names files")
workflow.add_argument("path-strat", desc="HUMAnN2 stratified pathways table")
workflow.add_argument("mpa", desc="MetaPhlAn2 table")

args = workflow.parse_args()

# subset microbiome

if args.subset_microbiome == "no":

	out_files_0 = ["days_0-1095." + args.microbiome]

	workflow.add_task(
		"samples_0-3_y.R --microbiome [depends[0]] --metadata [depends[1]]",
		depends=[args.microbiome, args.metadata],
		targets=out_files_0
		)

	if args.microbiome_feature == "EC" or "Pfam":

		out_files_0_a = ["days_0-1095." + args.mpa]

		workflow.add_task(
			"samples_0-3_y.R --microbiome [depends[0]] --metadata [depends[1]]",
			depends=[args.mpa, args.metadata],
			targets=out_files_0_a
			)

else:

	out_files_0 = ["days_" + str(args.min_a) + "-" + str(args.max_a) + "." + args.microbiome, "days_" + str(args.min_b) + "-" + str(args.max_b) + "." + args.microbiome]

	workflow.add_task(
		"subset_samples_by_age.R --microbiome [depends[0]] --metadata [depends[1]] --min-a [min_a] --max-a [max_a] --min-b [min_b] --max-b [max_b]",
		depends=[args.microbiome, args.metadata],
		targets=out_files_0,
		min_a=args.min_a,
		max_a=args.max_a,
		min_b=args.min_b,
		max_b=args.max_b
		)

	if args.microbiome_feature == "EC" or args.microbiome_feature == "Pfam":

		out_files_0_a = ["days_" + str(args.min_a) + "-" + str(args.max_a) + "." + args.mpa, "days_" + str(args.min_b) + "-" + str(args.max_b) + "." + args.mpa]

		workflow.add_task(
			"subset_samples_by_age.R --microbiome [depends[0]] --metadata [depends[1]] --min-a [min_a] --max-a [max_a] --min-b [min_b] --max-b [max_b]",
			depends=[args.mpa, args.metadata],
			targets=out_files_0_a,
			min_a=args.min_a,
			max_a=args.max_a,
			min_b=args.min_b,
			max_b=args.max_b
			)

# filter

prefix = []
for element in out_files_0:
        prefix.append(os.path.splitext(element)[0])

out_files_1 = []
for i in prefix:
        out_files_1.append(i + ".major.tsv")

workflow.add_task_group(
   "filter_ab_prev.R -i [depends[0]] -f [feature] -a [abundance] -p [prevalence] -t [transformation]",
    depends=out_files_0,
    targets=out_files_1,
    feature=args.microbiome_feature,
    abundance=args.abundance,
    prevalence=args.prevalence,
    transformation=args.transformation
    )

if args.microbiome_feature == "EC" or args.microbiome_feature == "Pfam":
	
	mpa_list = ["days_" + str(args.min_a) + "-" + str(args.max_a) + ".major_bugs_list.tsv", "days_" + str(args.min_b) + "-" + str(args.max_b) + ".major_bugs_list.tsv"]

	workflow.add_task_group(
		"filter_ab_prev_list.R -i [depends[0]] -a [abundance] -p [prevalence] -t [transformation]",
		depends=out_files_0_a,
		targets=mpa_list,
		abundance=args.abundance,
		prevalence=args.prevalence,
		transformation=args.transformation
		)

# format + run

if args.subset_genetics == "no":

	out_files_2 = []
	for i in prefix:
		out_files_2.append(i + ".major.model_input.tsv")

	workflow.add_task_group(
		"format_input.R -m [depends[0]] -g [genetics] -d [metadata] -o [targets[0]]",
		depends=out_files_1,
		targets=out_files_2,
		genetics=args.genetics,
		metadata=args.metadata
		)

	out_files_3 = []
	for i in prefix:
		out_files_3.append(i + ".major.model_output.tsv")

	workflow.add_task_group(
		"run_models.R -i [depends[0]] -c [cores] -o [targets[0]]",
		depends=out_files_2,
		targets=out_files_3,
		cores=args.cores
		)

else:
    
	count = ( len(open(args.genetics).readlines(  )) - 1 ) / args.size
    
	out_files_3 = []
	for i in range(1, count + 1):
		out_files_3.append("snps." + str(args.size) + ".batch_" + str(i) + ".txt")

	workflow.add_task(
		"subset.sh -i [depends[0]] -p [prefix] -n [chunk_size] -t [threads]",
		depends=args.genetics,
		targets=out_files_3,
		prefix="snps."+str(args.size),
		chunk_size=args.size,
		threads=args.threads
		)

	genetics_list = []
	for i in range(1, count + 1):
		genetics_list.append("snps." + str(args.size) + ".batch_" + str(i))

	suffix = []
	for element in out_files_0:
		suffix.append(element.strip(".tsv"))

	out_files_4 = []
	for i in genetics_list:
		for j in suffix:
			out_files_4.append(i + "." + j + ".major.model_input.tsv")

	workflow.add_task_group(
		"format_input.R -m [microbiome] -g [depends[0]] -d [metadata] -o [targets[0]]",
		depends=out_files_3,
		targets=out_files_4,
		microbiome=out_files_1,
		metadata=args.metadata
		)

	out_files_5 = []
	for i in genetics_list:
		for j in suffix:
			out_files_5.append(i + "." + j + ".major.model_output.tsv")

	if args.system == "server":

		workflow.add_task_group(
			"run_models.R -i [depends[0]] -c [cores] -o [targets[0]]",
			depends=out_files_4,
			targets=out_files_5,
			cores=args.cores
			)
	
	else:

		workflow.do_gridable(
			"run_models.R -i [depends[0]] -o [targets[0]]",
			depends=out_files_4,
			targets=out_files_5,
			cores=args.cores,
			time=args.time,
			mem=args.mem,
			partition=args.partition
			)

	if args.subset_microbiome == "no":

		out_files_6 = args.microbiome_feature + "_" + args.genetics_feature + "." + "days_0-1095.merged_results.tsv"

		workflow.add_task(
			"merge_model_outputs.py -i [extension] -o [targets[0]]",
			depends=out_files_5,
			targets=out_files_6,
			extension="model_output.tsv"
			)

	else:

		out_files_6_a = args.microbiome_feature + "_" + args.genetics_feature + ".days_" + str(args.min_a) + "-" + str(args.max_a) + ".merged_results.tsv"

		extension_6_a = "days_" + str(args.min_a) + "-" + str(args.max_a) + "." + os.path.splitext(args.microbiome)[0] + ".major.model_output.tsv"

		workflow.add_task(
			"merge_model_outputs.py -i [extension] -o [targets[0]]",
			depends=out_files_5,
			targets=out_files_6_a,
			extension=extension_6_a
			)

		out_files_6_b = args.microbiome_feature + "_" + args.genetics_feature + ".days_" + str(args.min_b) + "-" + str(args.max_b) + ".merged_results.tsv"

		extension_6_b = "days_" + str(args.min_b) + "-" + str(args.max_b) + "." + os.path.splitext(args.microbiome)[0] + ".major.model_output.tsv"

		workflow.add_task(
			"merge_model_outputs.py -i [extension] -o [targets[0]]",
			depends=out_files_5,
			targets=out_files_6_b,
			extension=extension_6_b
			)

# visualisation

plots = ["heatmap.png", "model_stats.tsv"]

if args.genetics_feature == "PC":
	
	if args.subset_genetics == "no":

		for i in out_files_3:

			data = i.replace("output", "input")

			out = i.replace(".major.model_output.tsv", ".figures")

			targets = []
			for j in plots:
				targets.append("./" + out + "/" + j)

			workflow.add_task(
				"plot_results_pcs.R -s [depends[0]] -d [model_data] -o [out] -p [padj] -g [group] -t [transformation] --humann [humann] -f [microbiome_feature]",
				depends=i,
				targets=targets,
				model_data=data,
				out=out,
				padj=args.padj,
				group=args.group_by,
				transformation=args.transformation,
				humann=args.path_names,
				microbiome_feature=args.microbiome_feature
				)

		# species contributions to "significant" pathways

		if args.microbiome_feature == "EC" or args.microbiome_feature == "Pfam":

			strat_names = args.path_strat.replace(".tsv", "_names.tsv")
			
			# extract the first column from the stratified HUMAnN2 table

			workflow.add_task(
				"cut -f1 [depends[0]]  | uniq > [targets[0]]",
				depends=args.path_strat,
				targets=strat_names
				)

			# pathways to names

			if args.microbiome_feature == "EC":
				map_names = args.path_names + "/map_level4ec_name.txt"

			if args.microbiome_feature == "Pfam":
				map_names = args.path_names + "/map_pfam_name.txt" 

			# summary of species contributions to "significant" pathways

			for i in out_files_3:

				mpa = i.replace("." + args.microbiome.replace(".tsv", ".major.model_output.tsv"), ".major_bugs_list.tsv")

				sig = i.replace("." + args.microbiome.replace(".tsv", ".major.model_output.tsv"), ".humann2_stratified_names.sig_" + args.microbiome_feature + ".tsv")

				workflow.add_task(
					"plot_path_smry.R --input [depends[0]] --feat [feat] --names [depends[1]] --strat [strat_names] --padj [padj] --mpa [mpa] --out [out] --group [group]",
					depends=[i, map_names],
					out=out,
					feat=args.microbiome_feature,
					strat_names=strat_names,
					targets=sig,
					mpa=mpa,
					padj=args.padj,
					group=args.group_by
					)

			# extract the header from the stratified HUMAnN2 table

			header = re.sub(r".*/", "", "temp_" + args.path_strat.replace(".tsv", ".header.tsv"))

			workflow.add_task(
				"head -1 [depends[0]] > [targets[0]]",
				depends=args.path_strat,
				targets=header
				)

			# extract "significant" pathhways from the stratified HUMAnN2 table

			sig_2 = []
			for i in out_files_3:
				sig_2.append(i.replace("." + args.microbiome.replace(".tsv", ".major.model_output.tsv"), ".humann2_stratified_names.sig_" + args.microbiome_feature + ".tsv"))

			temp_1 = []
			for i in sig_2:
				temp_1.append("temp_" + i.replace(".humann2_stratified_names.sig_" + args.microbiome_feature + ".tsv", ".humann2_stratified_abundances.sig_" + args.microbiome_feature + ".tsv"))

			workflow.add_task_group(
				"grep -f [depends[0]] [path_strat] > [targets[0]]",
				depends=sig_2,
				targets=temp_1,
				path_strat=args.path_strat
				)

			# combine the header + pathways

			res = []
			for i in temp_1:
				res.append(i.replace("temp_", ""))

			workflow.add_task(
				"cat [header] [depends[0]] > [targets[0]]",
				depends=temp_1,
				targets=res,
				header=header
				)

			# make stacked bar charts to show the contributions of species to "significant" pathways

			if args.subset_microbiome == "no":
				
				depends_no_sub = ["days_" + str(args.min_a) + "-" + str(args.max_b) + ".humann2_stratified_abundances.tsv"]

				workflow.add_tasl(
					"plot_path_bars.R --inp [depends[0]] --met [depends[1]] --names [depends[2]] --min [min] --max [max] --feat [feat]",
					depends=[depends_no_sub, args.metadata, map_names],
					min=args.min_a,
					max=args.max_b,
					feat=args.microbiome_feature
					)
			
			else:
				
				depends_a = ["days_" + str(args.min_a) + "-" + str(args.max_a) + ".humann2_stratified_abundances.tsv"]
				
				workflow.add_task(
					"plot_path_bars.R --inp [depends[0]] --met [metadata] --names [map_names] --min [min] --max [max] --feat [feat]",
					depends=res[0],
					metadata=args.metadata,
					map_names=map_names,
					min=args.min_a,
					max=args.max_a,
					feat=args.microbiome_feature
					)

				depends_b = ["days_" + str(args.min_b) + "-" + str(args.max_b) + ".humann2_stratified_abundances.tsv"]
				
				workflow.add_task(
					"plot_path_bars.R --inp [depends[0]] --met [metadata] --names [map_names] --min [min] --max [max] --feat [feat]",
					depends=[res[1]],
					metadata=args.metadata,
					map_names=map_names,
					min=args.min_b,
					max=args.max_b,
					feat=args.microbiome_feature
					)

# else:
# visualisation of "microbiome x SNP" associations

# run the workflow

workflow.go()

print(sig_2)
print(outfiles_3)
print(temp_1)

# END
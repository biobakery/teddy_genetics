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

workflow.add_argument("microbiome", desc="Microbiome data", required=True)
workflow.add_argument("microbiome-feature", desc="Type of microbiome feature", choices=["Bray-Curtis", "EC", "Pfam", "Species"], required=True)
workflow.add_argument("metadata", desc="Metadata for samples", required=True)
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
workflow.add_argument("path-names", desc="Path to directory containing HUMAnN2 pathway names files")
workflow.add_argument("path-strat", desc="HUMAnN2 stratified pathways table")
workflow.add_argument("mpa", desc="MetaPhlAn2 table")

args = workflow.parse_args()

# PART A

if (args.microbiome_feature == "EC" or args.microbiome == "Pfam") and (args.path_names is None or args.path_strat is None or args.mpa is None):
	print("EC / Pfam requires --path-names, --path-strat, and --mpa.")

else:

    # make an output directory

    directory = "./PC_models/" + args.microbiome_feature 
    if not os.path.exists(directory):
        os.makedirs(directory)

	# subset microbiome

    if args.subset_microbiome == "no":

        out_files_0 = [directory + "/days_0-1095." + args.microbiome]

    	workflow.add_task(
    		"samples_0-3_y.R --microbiome [depends[0]] --metadata [depends[1]]",
    		depends=[args.microbiome, args.metadata],
    		targets=out_files_0
    		)

    	if args.microbiome_feature == "EC" or "Pfam":

    		out_files_0_a = [directory + "days_0-1095." + args.mpa]

    		workflow.add_task(
    			"samples_0-3_y.R --microbiome [depends[0]] --metadata [depends[1]] -o [out_dir]",
    			depends=[args.mpa, args.metadata],
    			targets=out_files_0_a,
                out_dir=directory
    			)

    else:

    	out_files_0 = [directory + "/days_" + str(args.min_a) + "-" + str(args.max_a) + "." + args.microbiome, directory + "/days_" + str(args.min_b) + "-" + str(args.max_b) + "." + args.microbiome]

    	workflow.add_task(
    		"subset_samples_by_age.R --microbiome [depends[0]] --metadata [depends[1]] --min-a [min_a] --max-a [max_a] --min-b [min_b] --max-b [max_b] -o [out_dir]",
    		depends=[args.microbiome, args.metadata],
    		targets=out_files_0,
    		min_a=args.min_a,
    		max_a=args.max_a,
    		min_b=args.min_b,
    		max_b=args.max_b,
            out_dir=directory
    		)

        if args.microbiome_feature == "EC" or args.microbiome_feature == "Pfam":

            out_files_0_a = [directory + "/days_" + str(args.min_a) + "-" + str(args.max_a) + "." + args.mpa, directory + "/days_" + str(args.min_b) + "-" + str(args.max_b) + "." + args.mpa]

            workflow.add_task(
                "subset_samples_by_age.R --microbiome [depends[0]] --metadata [depends[1]] --min-a [min_a] --max-a [max_a] --min-b [min_b] --max-b [max_b] -o [out_dir]",
                depends=[args.mpa, args.metadata],
                targets=out_files_0_a,
                min_a=args.min_a,
                max_a=args.max_a,
                min_b=args.min_b,
                max_b=args.max_b,
                out_dir=directory
                )

    # filter

    prefix = []
    for i in out_files_0:
    	prefix.append(os.path.splitext(i)[0])

    out_files_1 = []
    for i in prefix:
    	out_files_1.append(i + ".major.tsv")

    if args.microbiome_feature != "Bray-Curtis":

    	workflow.add_task_group(
    		"filter_ab_prev.R -i [depends[0]] -m [metadata] -f [feature] -a [abundance] -p [prevalence] -t [transformation]",
    		depends=out_files_0,
    		targets=out_files_1,
    		metadata=args.metadata,
    		feature=args.microbiome_feature,
    		abundance=args.abundance,
    		prevalence=args.prevalence,
    		transformation=args.transformation
    		)

    if args.microbiome_feature == "EC" or args.microbiome_feature == "Pfam":

        mpa_list = []
        for i in out_files_0:
            mpa_list.append(i.replace(args.microbiome, "metaphlan2.major_list.tsv"))

    	workflow.add_task_group(
    		"filter_ab_prev_list.R -i [depends[0]] -m [metadata] -a [abundance] -p [prevalence] -o [out_dir]",
    		depends=out_files_0_a,
    		targets=mpa_list,
    		metadata=args.metadata,
    		abundance=args.abundance,
    		prevalence=args.prevalence,
            out_dir=directory
    		)

    # format + run

    out_files_2 = []
    for i in prefix:
    	out_files_2.append(i + ".major.model_input.tsv")

    if args.microbiome_feature == "Bray-Curtis":
    	my_dependency = out_files_0
    else:
    	my_dependency = out_files_1

    workflow.add_task_group(
    	"format_input.R -m [depends[0]] -f [feature] -g [genetics] -d [metadata] -o [targets[0]]",
    	depends=my_dependency,
    	targets=out_files_2,
    	feature=args.microbiome_feature,
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

    # visualisation

    plots = ["heatmap.png", "model_stats.tsv"]

    for i in out_files_3:
	    data = i.replace("output", "input")
	    out = i.replace(".major.model_output.tsv", ".figures")

	    targets = []
	    for j in plots:
		    targets.append(out + "/" + j)

	    workflow.add_task(
		    "plot_results_pcs.R -s [depends[0]] -d [model_data] -o [out] -p [padj] -g [group] -t [transformation] -f [microbiome_feature] --humann [humann]",
		    depends=i,
		    targets=targets,
		    model_data=data,
		    out=out,
		    padj=args.padj,
		    group=args.group_by,
		    transformation=args.transformation,
		    microbiome_feature=args.microbiome_feature,
            humann=args.path_names
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

        sig_list = []

        for i in out_files_3:

            mpa = i.replace(args.microbiome.replace(".tsv", ".major.model_output.tsv"), "metaphlan2.major_list.tsv")

            sig = i.replace(args.microbiome.replace(".tsv", ".major.model_output.tsv"), "humann2_stratified_names.sig_" + args.microbiome_feature + ".tsv")

            sig_list.append(sig)

            workflow.add_task(
			    "plot_path_smry.R --input [depends[0]] --feat [feat] --names [depends[1]] --strat [strat_names] --padj [padj] --mpa [mpa] --out [out] --group [group]",
			    depends=[i, args.path_names],
			    out=out,
			    feat=args.microbiome_feature,
			    strat_names=strat_names,
			    targets=sig,
			    mpa=mpa,
			    padj=args.padj,
			    group=args.group_by
			    )

    if args.microbiome_feature == "EC" or args.microbiome=="Pfam":
        title = "# " + args.microbiome_feature + "-PC models (PART A) #\n"
    else:
        title = "# " + args.microbiome_feature + "-PC models #\n"
    
    print("#" * (len(title) - 1) + "\n" + title + "#" * (len(title) - 1))
    workflow.go()

    # PART B

    if args.microbiome_feature == "EC" or args.microbiome_feature == "Pfam":

        for i in sig_list:

            file_size = os.stat(i).st_size

            if file_size > 0:

                header = directory + "/" + re.sub(r".*/", "", "temp_" + args.path_strat.replace(".tsv", ".header.tsv"))

                workflow.add_task(
                    "head -1 [depends[0]] > [targets[0]]",
                    depends=args.path_strat,
                    targets=header
                    )

                temp_1 = i.replace("days", "temp_days")

                workflow.add_task(
                    "grep -f [depends[0]] [path_strat] > [targets[0]]",
                    depends=i,
                    targets=temp_1,
                    path_strat=args.path_strat
                    )

                res = temp_1.replace("temp_", "")

                workflow.add_task(
                    "cat [header] [depends[0]] > [targets[0]]",
                    depends=temp_1,
                    targets=res,
                    header=header
                    )

                if args.subset_microbiome == "no":

                    depends_no_sub = ["days_" + str(args.min_a) + "-" + str(args.max_b) + ".humann2_stratified_abundances.tsv"]

                    workflow.add_task(
  						"plot_path_bars.R --inp [depends[0]] --met [depends[1]] --names [depends[2]] --min [min] --max [max] --feat [feat] --out [out_dir]",
  						depends=[depends_no_sub, args.metadata, arg.path_names],
  						min=args.min_a,
  						max=args.max_b,
  						feat=args.microbiome_feature,
                        out_dir=directory
  						)

                else:

                    if ("days_" + str(args.min_a) + "-" + str(args.max_a)) in res:

                        targets_bars_a = directory + "/days_" + str(args.min_a) + "-" + str(args.max_a) + "_" + args.microbiome_feature + "_stacked_bar_charts/"
                        
                        workflow.add_task(
                            "plot_path_bars.R --inp [depends[0]] --met [metadata] --names [map_names] --min [min] --max [max] --feat [feat] --out [out_dir]",
                            depends=res,
                            targets=targets_bars_a,
                            metadata=args.metadata,
                            map_names=args.path_names,
							min=args.min_a,
							max=args.max_a,
							feat=args.microbiome_feature,
                            out_dir=directory
							)

                    if ("days_" + str(args.min_b) + "-" + str(args.max_b)) in res:

                        targets_bars_b = directory + "/days_" + str(args.min_b) + "-" + str(args.max_b) + "_" + args.microbiome_feature + "_stacked_bar_charts/"

                        workflow.add_task(
							"plot_path_bars.R --inp [depends[0]] --met [metadata] --names [map_names] --min [min] --max [max] --feat [feat] --out [out_dir]",
							depends=res,
                            targets=targets_bars_b,
							metadata=args.metadata,
							map_names=args.path_names,
							min=args.min_b,
							max=args.max_b,
							feat=args.microbiome_feature,
                            out_dir=directory
							)

        title_b = "# " + args.microbiome_feature + "-PC models (PART B) #\n"
        print("#" * (len(title_b) - 1) + "\n" + title_b + "#" * (len(title_b) - 1))
        workflow.go()

# END

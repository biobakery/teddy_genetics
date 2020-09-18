#!/usr/bin/env python

from anadama2 import Workflow
import glob
import subprocess

class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end

workflow = Workflow(version="0.1", description="A workflow to format input for TEDDY2 models")

workflow.add_argument("microbiome", desc="Microbiome data", required=True)
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
workflow.add_argument("subset-genetics", desc="Should the genetics data be subset into chunks?", choices=["yes", "no"], default="no")
workflow.add_argument("size", desc="Size of the chunks into which the genetics data are subsetted", type=int, default=250)
workflow.add_argument("threads", desc="Number of threads to use when splitting the genetics data into subsets", type=int, default=1)

args = workflow.parse_args()

if args.subset_microbiome == "no":

	out_files_0 = ["days_0-1095." + args.microbiome]

	workflow.add_task(
		"samples_0-3_y.R --microbiome [depends[0]] --metadata [depends[1]]",
        depends=[args.microbiome, args.metadata],
        targets=out_files_0
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

prefix = []
for element in out_files_0:
	prefix.append(element.strip(".tsv"))

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
			out_files_4.append(i + "." + j + ".tsv")

	for i in out_files_1:
		workflow.add_task_group(
			"format_input.R -m [microbiome] -g [depends[0]] -d [metadata] -o [targets[0]]",
			depends=out_files_3,
			targets=out_files_4,
			microbiome=i,
			metadata=args.metadata
			)

# run the workflow

workflow.go()

# remove intermediary files

subprocess.call("rm *.txt", shell=True)

# END
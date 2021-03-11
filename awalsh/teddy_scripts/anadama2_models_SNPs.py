#!/usr/bin/env python

from anadama2 import Workflow
import os
import math
import decimal

class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end

workflow = Workflow(version="0.1", description="A workflow to run SNP models for TEDDY2")

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
workflow.add_argument("subset-genetics", desc="Should the genetics data be subset into chunks?", choices=["yes", "no"], default="no")
workflow.add_argument("size", desc="Size of the chunks into which the genetics data are subsetted", type=int, default=250)
workflow.add_argument("threads", desc="Number of threads to use when splitting the genetics data into subsets", type=int, default=1)
workflow.add_argument("cores", desc="Number of processors to use for modeling on cluster", type=int, default=1)
workflow.add_argument("time", desc="Limit on runtime (minutes) of jobs on cluster", type=int, default=720)
workflow.add_argument("mem", desc="Amount of memory (MB) to allocate to processors on cluster", type=int, default=4000)
workflow.add_argument("partition", desc="Partition to which jobs are submitted on cluster", default="serial_requeue")
workflow.add_argument("mpa", desc="MetaPhlAn2 table")

args = workflow.parse_args()

######################################
# part 1: format the microbiome data #
######################################

# make an output directory

directory = "./microbiome_formatted/" + args.microbiome_feature 
if not os.path.exists(directory):
	os.makedirs(directory)

micro = args.microbiome.split("/")[-1:][0]

# subset microbiome

if args.subset_microbiome == "no":

	out_files_0 = [directory + "/days_0-1095." + micro]

	workflow.add_task(
		"samples_0-3_y.R --microbiome [depends[0]] --metadata [depends[1]]",
		depends=[args.microbiome, args.metadata],
		targets=out_files_0
		)

	if args.microbiome_feature == "EC" or "Pfam":

		mpa = args.mpa.split("/")[-1:][0]

		out_files_0_a = [directory + "days_0-1095." + mpa]

		workflow.add_task(
			"samples_0-3_y.R --microbiome [depends[0]] --metadata [depends[1]] -o [out_dir]",
			depends=[args.mpa, args.metadata],
			targets=out_files_0_a,
			out_dir=directory
			)

else:

	out_files_0 = [directory + "/days_" + str(args.min_a) + "-" + str(args.max_a) + "." + micro, directory + "/days_" + str(args.min_b) + "-" + str(args.max_b) + "." + micro]

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

		mpa = args.mpa.split("/")[-1:][0]

		out_files_0_a = [directory + "/days_" + str(args.min_a) + "-" + str(args.max_a) + "." + mpa, directory + "/days_" + str(args.min_b) + "-" + str(args.max_b) + "." + mpa]

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

#####################################
# part 2: subset the genetics table #
#####################################

snp_chunks = "./SNP_models/snp_chunks/" 
if not os.path.exists(snp_chunks):
	os.makedirs(snp_chunks)

count = 0
for line in open(args.genetics).xreadlines(  ): count += 1

chunk_size = args.size
genetics = args.genetics.split("/")[-1:][0]

if count > chunk_size:

    prefix = snp_chunks + "snp_chunk"

    n = int( math.ceil( decimal.Decimal(count - 1) / decimal.Decimal(chunk_size) ) )

    targets_subsets = []
    for i in range(n):
    	targets_subsets.append(prefix + ".batch_" + str(i + 1) + ".txt")

    workflow.add_task(
        "subset.sh -i [depends[0]] -n [chunk_size] -p [prefix] -t [threads]",
        depends=args.genetics,
        targets=targets_subsets,
        chunk_size=chunk_size,
        prefix=prefix,
        threads=args.threads
        )

else:

    targets_subsets = [snp_chunks + genetics.replace(".tsv", ".txt")]

    for i in targets_subsets:

    	gen = i.replace(".tsv", ".txt")
    	
    	workflow.add_task("cp [input] [targets[0]]",
    		targets=gen,
    		input=args.genetics)
    
##########################
# part 3: run the models #
##########################

results = "./SNP_models/" + args.microbiome_feature 
if not os.path.exists(results):
	os.makedirs(results)

for i in out_files_1:
	for j in targets_subsets:

		part_a = i.split("/")[-1:][0].split(".")[0] + "." + args.microbiome_feature
		part_b = j.split("/")[-1:][0].replace(".txt", ".model_output.tsv.gz")
		targets = results + "/" + part_a + "." + part_b

		workflow.add_task_gridable(
			"run_models_snps.R -m [depends[0]] -f [microbiome_feature] -g [depends[1]] -d [metadata] -c [cpus] -o [results_dir]",
			depends=[i, j],
			targets=targets,
			microbiome_feature=args.microbiome_feature,
			metadata=args.metadata,
			cpus=args.cores,
			cores=args.cores,
			time=args.time,
			mem=args.mem,
			partition=args.partition,
			results_dir=results
			)

############################
# part 4: run the workflow #
############################

workflow.go()

for i in range(n):
	print(i)

#

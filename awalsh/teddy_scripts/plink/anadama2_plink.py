#!/usr/bin/env python

from anadama2 import Workflow
import subprocess
import os

# make tmp directory

directory = "./tmp"

if not os.path.exists(directory):
    os.makedirs(directory)

# make output directory

directory = "./plink_output"

if not os.path.exists(directory):
    os.makedirs(directory)

# add arguments

workflow = Workflow(version="0.1", description="A workflow to process genotype data for TEDDY2 using PLINK")

workflow.add_argument("genetics", desc="Prefix of raw genotype data", required=True)
workflow.add_argument("subjects", desc="Subjects to include")
workflow.add_argument("prefix", desc="Prefix of output", default="plink")
workflow.add_argument("maf", desc="Minor allele frequency threshold", type=float, default=0.05)
workflow.add_argument("hwe", desc="Hardy-Weinberg equilibrium theshold", type=float, default=1E-3)
workflow.add_argument("pihat", desc="PI_HAT threshold", type=float, default=0.2)
workflow.add_argument("ld-prune-pca", desc="Do LD pruning be done prior to PCA", choices=["true", "false"], default="true")
workflow.add_argument("npcs", desc="The number of PCs to use", type=int, default=20)
workflow.add_argument("window", desc="Window size (kb) for LD pruning", type=int, default=50)
workflow.add_argument("step", desc="Step size for LD pruning", type=int, default=5)
workflow.add_argument("r-squared", desc="r^2 threshold for LD pruning", type=float, default=0.5)
workflow.add_argument("metadata", desc="Metadata for samples", required=True)

args = workflow.parse_args()

# 1 - select the samples that are to be included in the analysis

depends_1 = []
for i in [".bed", ".bim", ".fam"]:
	depends_1.append(args.genetics + i)

targets_1 = ["tmp/temp_1.bed", "tmp/temp_1.bim", "tmp/temp_1.fam"]

workflow.add_task(
	"plink --bfile [prefix_in] --keep [subjects] --make-bed --set-hh-missing --out [prefix_out]",
	depends=depends_1,
	targets=targets_1,
	prefix_in=args.genetics,
	subjects=args.subjects,
	prefix_out="tmp/temp_1"
	)

# 2 - run KING autoQC

depends_2 = "tmp/temp_1.bed"

targets_2 = ["tmp/temp_autoQC_snptoberemoved.txt", "tmp/temp_autoQC_sampletoberemoved.txt", "tmp/temp_autoQC_Summary.txt"]

workflow.add_task(
	"king -b [depends[0]] --autoQC --prefix [prefix]",
	depends="tmp/temp_1.bed",
	targets=targets_2,
	prefix="tmp/temp"
	)

# 3 - remove SNPs that failed KING autoQC

depends_3 = ["tmp/temp_1.bed", "tmp/temp_1.bim", "tmp/temp_1.fam", "tmp/temp_autoQC_snptoberemoved.txt"]

targets_3 = ["tmp/temp_2.bed", "tmp/temp_2.bim", "tmp/temp_2.fam"]

workflow.add_task(
	"plink --bfile [prefix_in] --exclude [snps_to_exclude] --make-bed --out [prefix_out]",
	depends=depends_3,
	targets=targets_3,
	prefix_in="tmp/temp_1",
	snps_to_exclude="tmp/temp_autoQC_snptoberemoved.txt",
	prefix_out="tmp/temp_2"
	)

# 4 - remove samples that failed KING autoQC

depends_4 = ["tmp/temp_2.bed", "tmp/temp_2.bim", "tmp/temp_2.fam", "tmp/temp_autoQC_sampletoberemoved.txt"]

targets_4 = ["tmp/temp_3.bed", "tmp/temp_3.bim", "tmp/temp_3.fam"]

workflow.add_task(
	"plink --bfile [prefix_in] --exclude [exclude] --make-bed --out [prefix_out]",
	depends=depends_4,
	targets=targets_4,
	prefix_in="tmp/temp_2",
	exclude="tmp/temp_autoQC_sampletoberemoved.txt",
	prefix_out="tmp/temp_3"
	)

# 5 - reformat temp_autoQC_Summary.txt

workflow.do("sed -e 's/  \\+/\t/g' [d:tmp/temp_autoQC_Summary.txt] | cut -f2,3,4 | tail -n +2 > [t:tmp/temp_autoQC_smry.txt]")

# 6 - summary from KING autoQC

targets_6 = ["tmp/" + args.prefix + "_king_autoQC_Summary.txt", "tmp/" + args.prefix + "_king_autoQC_Breakdown.txt"]

workflow.add_task(
	"plink_king_tidy.R -p tmp/[prefix]",
	depends="tmp/temp_autoQC_smry.txt",
	targets=targets_6,
	prefix=args.prefix
	)

# 7 - caclulate MAF

depends_7 = ["tmp/temp_3.bed", "tmp/temp_3.bim", "tmp/temp_3.fam", "tmp/" + args.prefix + "_king_autoQC_Summary.txt", "tmp/" + args.prefix + "_king_autoQC_Breakdown.txt"]

targets_7 = ["tmp/" + args.prefix + ".frq"]

workflow.add_task(
	"plink --bfile [prefix_in] --freq --out [prefix_out]",
	depends=depends_7,
	targets=targets_7,
	prefix_in="tmp/temp_3",
	prefix_out=["tmp/" + args.prefix]
	)

# 8 - remove SNPs with MAF < args.maf

depends_8 = ["tmp/temp_3.bed", "tmp/temp_3.bim", "tmp/temp_3.fam", "tmp/" + args.prefix + ".frq"]

targets_8 = []
for i in [".bed", ".bim", ".fam"]:
	targets_8.append("tmp/temp_4.maf-" + str(args.maf) + i)

prefix_in_8 = "tmp/temp_3"

prefix_out_8 = "tmp/temp_4.maf-" + str(args.maf)

workflow.add_task(
	"plink --bfile [prefix_in] --maf [maf] --make-bed --out [prefix_out]",
	depends=depends_8,
	targets=targets_8,
	prefix_in=prefix_in_8,
	maf=args.maf,
	prefix_out=prefix_out_8
	)

# 9 - caclulate HWE

depends_9 = targets_8

targets_9 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe"

prefix_in_9 = "tmp/temp_4.maf-" + str(args.maf)

prefix_out_9 = "tmp/" + args.prefix + ".maf-" + str(args.maf)

workflow.add_task(
	"plink --bfile [prefix_in] --hardy --out [prefix_out]",
	depends=depends_9,
	targets=targets_9,
	prefix_in=prefix_in_9,
	prefix_out=prefix_out_9
	)

# 10 - remove samples that deviate from HWE

depends_10 = targets_9

targets_10 = []
for i in [".bed", ".bim", ".fam"]:
	targets_10.append("tmp/temp_5.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + i)

prefix_in_10 = "tmp/temp_4.maf-" + str(args.maf)

prefix_out_10 = "tmp/temp_5.maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

workflow.add_task(
	"plink --bfile [prefix_in] --hwe [hwe] --not-chr x --make-bed --out [prefix_out]",
	depends=depends_10,
	targets=targets_10,
	prefix_in=prefix_in_10,
	hwe=args.hwe,
	prefix_out=prefix_out_10
	)

# 11 - prune LD (prior to heterozygosity calculation)

depends_11 = targets_10

targets_11 = "tmp/temp_5.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".prune.in"

prefix_11 = "tmp/temp_5.maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

workflow.add_task(
	"plink --bfile [prefix] --indep-pairwise [window_size] [step_size] [r2_threshold] --out [prefix]",
	depends=depends_11,
	targets=targets_11,
	prefix=prefix_11,
	window_size=args.window,
	step_size=args.step,
	r2_threshold=args.r_squared
	)

# 12 - caclulate heterozygosity

depends_12 = targets_11

targets_12 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".het"

prefix_in_12 = "tmp/temp_5.maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

prefix_out_12 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

workflow.add_task(
	"plink --bfile [prefix_in] --extract [depends[0]] --het --out [prefix_out]",
	depends=depends_12,
	targets=targets_12,
	prefix_in=prefix_in_12,
	prefix_out=prefix_out_12
	)

# 13 - identify samples that deviate from heterozygosity

depends_13 = targets_12

targets_13 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".qc-fail-het.txt"

workflow.add_task(
	"plink_het.R -i [depends[0]] -o [targets[0]]",
	depends=depends_13,
	targets=targets_13
	)

# 14 - remove samples that deviate from heterozygosity

depends_14 = targets_13

targets_14 = []
for i in [".bed", ".bim", ".fam"]:
	targets_14.append("tmp/temp_6.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + i)

prefix_in_14 = "tmp/temp_5.maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

prefix_out_14 = "tmp/temp_6.maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

workflow.add_task(
	"plink --bfile [prefix_in] --remove [depends[0]] --make-bed --out [prefix_out]",
	depends=depends_14,
	targets=targets_14,
	prefix_in=prefix_in_14,
	prefix_out=prefix_out_14
	)

# 15 - prune LD (prior to IBD calculation)

depends_15 = targets_14

targets_15 = "tmp/temp_6.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".prune.in"

prefix_15 = "tmp/temp_6.maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

workflow.add_task(
	"plink --bfile [prefix] --indep-pairwise [window_size] [step_size] [r2_threshold] --out [prefix]",
	depends=depends_15,
	targets=targets_15,
	prefix=prefix_15,
	window_size=args.window,
	step_size=args.step,
	r2_threshold=args.r_squared
	)

# 16 - calculate IBD

depends_16 = targets_15

targets_16 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".genome"

prefix_in_16 = "tmp/temp_6.maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

prefix_out_16 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

workflow.add_task(
	"plink --bfile [prefix_in] --extract [depends[0]] --genome --out [prefix_out]",
	depends=depends_16,
	targets=targets_16,
	prefix_in=prefix_in_16,
	prefix_out=prefix_out_16
	)

# 17 - missingness per individual

depends_17 = [targets_15, targets_16]

targets_17 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".imiss"

prefix_out_17 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

workflow.add_task(
	"plink --bfile [prefix_in] --extract [depends[0]] --missing --out [prefix_out]",
	depends=depends_17,
	targets=targets_17,
	prefix_in=prefix_in_16,
	prefix_out=prefix_out_17
	)

# 18 - identify related individuals

depends_18 = [targets_16, targets_17]

targets_18 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc-fail-ibd.txt"

workflow.add_task(
	"plink_id_relatives.R -g [depends[0]] -m [depends[1]] -t [pihat] -o [targets[0]]",
	depends=depends_18,
	targets=targets_18,
	pihat=args.pihat
	)

# 19 - remove related individuals

depends_19 = targets_18

targets_19 = []
for i in [".bed", ".bim", ".fam"]:
	targets_19.append("tmp/temp_7.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + i)

prefix_in_19 = "tmp/temp_6.maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

prefix_out_19 = "tmp/temp_7.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat)

workflow.add_task(
	"plink --bfile [prefix_in] --remove [depends[0]] --make-bed --out [prefix_out]",
	depends=depends_19,
	targets=targets_19,
	prefix_in=prefix_in_19,
	prefix_out=prefix_out_19
	)

# 22 - identify outliers

depends_22 = [targets_16, targets_18]

targets_22 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc-fail-ibs.txt"

workflow.add_task(
	"plink_ibs_filter.R -i [depends[0]] -o [targets[0]]",
	depends=depends_22,
	targets=targets_22
	)

# 23 - remove outliers

depends_23 = targets_22

targets_23 = []
for i in [".bed", ".bim", ".fam", ".map", ".ped"]:
	targets_23.append("tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc" + i)

prefix_in_23 = "tmp/temp_7.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat)

prefix_out_23 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc"

workflow.add_task(
	"plink --bfile [prefix_in] --keep [depends[0]] --make-bed --recode tab --out [prefix_out]",
	depends=depends_23,
	targets=targets_23,
	prefix_in=prefix_in_23,
	prefix_out=prefix_out_23
	)

# 40 - make the genetics table

depends_24 = targets_23

targets_24 = "plink_output/genetics.SNPs.tsv"

workflow.add_task(
	"plink_make_genetics_table.R -i [prefix] -o [targets[0]]",
	depends=depends_24,
	prefix="tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc",
	targets=targets_24
	)

# 26/27 - PCA

if args.ld_prune_pca == "true":

	# 26

	depends_26 = targets_23
	
	targets_26 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc" + ".prune.in"
	
	prefix_26 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc"

	workflow.add_task(
		"plink --bfile [prefix] --indep-pairwise [window_size] [step_size] [r2_threshold] --out [prefix]",
		depends=depends_26,
		targets=targets_26,
		prefix=prefix_26,
		window_size=args.window,
		step_size=args.step,
		r2_threshold=args.r_squared
		)
	
	# 27
	
	depends_27 = targets_26

	targets_27 = []
	for i in [".excl_outliers.eigenval", ".excl_outliers.eigenvec", ".excl_outliers.eigenvec.var"]:
		targets_27.append("plink_output/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc" + i)

	prefix_in_27 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc"

	prefix_out_27 = "plink_output/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.excl_outliers"

	workflow.add_task(
		"plink --bfile [prefix_in] --extract [depends[0]] --pca header tabs var-wts --out [prefix_out]",
		depends=depends_27,
		targets=targets_27,
		prefix_in=prefix_in_27,
		prefix_out=prefix_out_27
		)

else:

	depends_27 = targets_23

	targets_27 = []
	for i in [".excl_outliers.eigenval", ".excl_outliers.eigenvec", ".excl_outliers.eigenvec.var"]:
		targets_27.append("plink_output/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc" + i)
	
	prefix_in_27 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc"

	prefix_out_27 = "plink_output/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.excl_outliers"
	
	workflow.add_task(
		"plink --bfile [prefix_in] --pca header tabs var-wts --out [prefix_out]",
		depends=depends_27,
		targets=targets_27,
		prefix_in=prefix_in_27,
		prefix_out=prefix_out_27
		)

# 28 - select the first n PCs (n=args.npcs)

depends_28 = targets_27[1]

targets_28 = "plink_output/" + args.prefix + ".PCA_PC1-PC" + str(args.npcs) + ".tsv"

workflow.add_task(
	"cut -f2-" + str(args.npcs + 2) + " [depends[0]] > [targets[0]]",
	depends=depends_28,
	targets=targets_28
	)

# 29 - identify haplotypes

depends_29 = targets_23

targets_29 = []
for i in [".blocks", ".blocks.det"]:
	targets_29.append("tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc" + i)

prefix_29 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc"

if not os.path.exists(targets_29[0]):
	workflow.add_task_group(
		"plink --bfile [prefix] --blocks 'no-pheno-req' --out [prefix]",
		depends=depends_29,
		targets=targets_29,
		prefix=prefix_29
		)

# 30 - plots

depends_30 = targets_27

targets_30 = []
for i in [".PLINK_PCA_fig.png", ".QC_bar_summary.png", ".QC_histograms.png"]:
	targets_30.append("plink_output/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + i)

king_file = "tmp/" + args.prefix + "_king_autoQC_Summary.txt"
frq_file = "tmp/" + args.prefix + ".frq"
hwe_file = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe"
genome_file = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".genome"
het_file = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".het"
pca_files = "plink_output/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat)

workflow.add_task(
	"plink_plots.R --ki [king_file] --mi [frq_file] --mt [maf_thresh] --ei [hwe_file] --et [hwe_thresh] --gi [genome_file] --gt [ibd_thresh] --hi [het_file] --pc [pca_files] --md [metadata]",
	depends=depends_30,
	targets=targets_30,
	king_file=king_file,
	frq_file=frq_file,
	maf_thresh=args.maf,
	hwe_file=hwe_file,
	hwe_thresh=args.hwe,
	genome_file=genome_file,
	ibd_thresh=args.pihat,
	het_file=het_file,
	pca_files=pca_files,
	metadata=args.metadata
	)

# 31 - download + unzip + sort the reference gff for humans

targets_31_a = "tmp/GCF_000001405.39_GRCh38.p13_genomic.gff.gz"

if not os.path.exists("tmp/GCF_000001405.39_GRCh38.p13_genomic.gff"):
	workflow.add_task(
		"wget --directory-prefix tmp https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
		targets=targets_31_a,
		)

	depends_31_b = targets_31_a
	targets_31_b = "tmp/temp_GCF_000001405.39_GRCh38.p13_genomic.gff.gz"

	workflow.add_task("mv [depends[0]] [targets[0]]",
		depends=depends_31_b,
		targets=targets_31_b)

	depends_31_c = targets_31_b
	targets_31_c = "tmp/temp_GCF_000001405.39_GRCh38.p13_genomic.gff"

	workflow.add_task(
		"gunzip [depends[0]]",
		depends=targets_31_b,
		targets=targets_31_c,
		)

	workflow.do("grep 'HGNC' [d:tmp/temp_GCF_000001405.39_GRCh38.p13_genomic.gff] > [t:tmp/GCF_000001405.39_GRCh38.p13_genomic.gff]")

depends_31_d = "tmp/GCF_000001405.39_GRCh38.p13_genomic.gff"
targets_31_d = "tmp/GCF_000001405.39_GRCh38.p13_genomic.sorted.gff"

workflow.add_task("bedtools sort -i [depends[0]] > [targets[0]]",
	depends=depends_31_d,
	targets=targets_31_d
	)

# 34 - convert our .map file to a .bed file

targets_34_a = "tmp/GCF_000001405.39_GRCh38.p13_assembly_report.txt"

if not os.path.exists(targets_34_a):
	workflow.add_task(
		"wget --directory-prefix tmp https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt",
		targets=targets_34_a,
		)

	depends_34_b = targets_34_a
	targets_34_b = "tmp/temp_GCF_000001405.39_GRCh38.p13_assembly_report.txt"

	workflow.add_task("mv [depends[0]] [targets[0]]",
		depends=depends_34_b,
		targets=targets_34_b)

	depends_34_c = targets_34_b
	targets_34_c = "tmp/GCF_000001405.39_GRCh38.p13_assembly_report.txt"

	workflow.add_task("tail -n +63 [depends[0]] > [targets[0]]",
		depends=depends_34_c,
		targets=targets_34_c)

depends_34 = ["tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.map", "tmp/GCF_000001405.39_GRCh38.p13_assembly_report.txt"]

rsid_positions_34 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.map"

targets_34 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.map_bed"

workflow.add_task(
	"plink_map.R -i [rsid_positions] -o [targets[0]]",
	depends=depends_34,
	targets=targets_34,
	rsid_positions=rsid_positions_34
	)

# 35 - sort our .bed file

depends_35 = targets_34

targets_35 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.map_sorted.bed"

rsid_positions_35 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.map_bed"

workflow.add_task(
	"bedtools sort -i [rsid_positions] > [targets[0]]",
	depends=depends_35,
	targets=targets_35,
	rsid_positions=rsid_positions_35
	)

# 36 - identify closest genes

depends_36 = [targets_35, targets_31_d]

targets_36 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.closest_gene.bed"

workflow.add_task(
	"bedtools closest -a [depends[0]] -b [depends[1]] -d > [targets[0]]",
	depends=depends_36,
	targets=targets_36
	)

# 37 - download HGNC data

if not os.path.exists("tmp/temp_hgnc_complete_set.txt"):

	workflow.add_task(
		"wget --directory-prefix tmp ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt",
		targets="tmp/hgnc_complete_set.txt",
		)

	workflow.do("mv [d:tmp/hgnc_complete_set.txt] [t:tmp/temp_hgnc_complete_set.txt]")

	workflow.do("cut -f1-3 [d:tmp/temp_hgnc_complete_set.txt] > [t:tmp/hgnc_complete_set.txt]")

# 38 - annotate SNPs

depends_38 = [targets_36, "tmp/hgnc_complete_set.txt"]

targets_38 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.rsid_to_gene.map.txt"

workflow.add_task(
	"plink_rsid_to_gene.R -i [depends[0]] -o [targets[0]]",
	depends=depends_38,
	targets=targets_38
	)

# 39 - visualise the PC loadings of SNPs

depends_39 = [targets_27[1], targets_38]

rsid_loadings_39 = "plink_output/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.excl_outliers.eigenvec.var"
rsid_names_39 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.rsid_to_gene.map.txt"
rsid_positions_39 = "tmp/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.map"

targets_39 = "plink_output/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.pc_loadings_manhattan.png"

prefix_out_39 = "plink_output/" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc"

workflow.add_task(
	"plink_pc_loadings_viz.R -i [loadings] -a [names] -m [positions] -o [prefix]",
	depends=depends_39,
	targets=targets_39,
	loadings=rsid_loadings_39,
	names=rsid_names_39,
	positions=rsid_positions_39,
	prefix=prefix_out_39
	)

# 40 - rsid

depends_40 = [rsid_names_39, rsid_positions_39]

targets_40 = "./plink_output/rsid_annotations_coordinates.tsv"

workflow.add_task(
	"make_rsid_map.R --map [rsid_positions] --rsid [rsid_position]",
	depends=depends_40,
	targets=targets_40,
	rsid_positions=rsid_names_39,
	rsid_position=rsid_positions_39
	)

# run the workflow

print("##################\n# PLINK WORKFLOW #\n##################")
workflow.go()

# END

#!/usr/bin/env python

from anadama2 import Workflow
import subprocess
import os

# make tmp directory OR move files from the tmp directory to the working directory 

directory = "./tmp"

if not os.path.exists(directory):
    os.makedirs(directory)

subprocess.call("mv tmp/* .", shell=True)

# add arguments

workflow = Workflow(version="0.1", description="A workflow to process genotype data for TEDDY2 using PLINK")

workflow.add_argument("genetics", desc="Prefix of raw genotype data", required=True)
workflow.add_argument("subjects", desc="Subjects to include")
workflow.add_argument("prefix", desc="Prefix of output", default="plink")
workflow.add_argument("maf", desc="Minor allele frequency threshold", type=float, default=0.05)
workflow.add_argument("hwe", desc="Hardy-Weinberg equilibrium theshold", type=float, default=1E-3)
workflow.add_argument("pihat", desc="PI_HAT threshold", type=float, default=0.2)
workflow.add_argument("ld-prune-pca", desc="Do LD pruning be done prior to PCA", choices=["true", "false"], default="false")
workflow.add_argument("npcs", desc="The number of PCs to use", type=int, default=5)

args = workflow.parse_args()

# 1 - select the samples that are to be included in the analysis

depends_1 = []
for i in [".bed", ".bim", ".fam"]:
	depends_1.append(args.genetics + i)

targets_1 = ["temp_1.bed", "temp_1.bim", "temp_1.fam"]

workflow.add_task(
	"plink --bfile [prefix_in] --keep [subjects] --make-bed --set-hh-missing --out [prefix_out]",
	depends=depends_1,
	targets=targets_1,
	prefix_in=args.genetics,
	subjects=args.subjects,
	prefix_out="temp_1"
	)

# 2 - run KING autoQC

depends_2 = "temp_1.bed"

targets_2 = ["temp_autoQC_snptoberemoved.txt", "temp_autoQC_sampletoberemoved.txt", "temp_autoQC_Summary.txt"]

workflow.add_task(
	"king -b [depends[0]] --autoQC --prefix [prefix]",
	depends="temp_1.bed",
	targets=targets_2,
	prefix="temp"
	)

# 3 - remove SNPs that failed KING autoQC

depends_3 = ["temp_1.bed", "temp_1.bim", "temp_1.fam", "temp_autoQC_snptoberemoved.txt"]

targets_3 = ["temp_2.bed", "temp_2.bim", "temp_2.fam"]

workflow.add_task(
	"plink --bfile [prefix_in] --exclude [snps_to_exclude] --make-bed --out [prefix_out]",
	depends=depends_3,
	targets=targets_3,
	prefix_in="temp_1",
	snps_to_exclude="temp_autoQC_snptoberemoved.txt",
	prefix_out="temp_2"
	)

# 4 - remove samples that failed KING autoQC

depends_4 = ["temp_2.bed", "temp_2.bim", "temp_2.fam", "temp_autoQC_sampletoberemoved.txt"]

targets_4 = ["temp_3.bed", "temp_3.bim", "temp_3.fam"]

workflow.add_task(
	"plink --bfile [prefix_in] --exclude [exclude] --make-bed --out [prefix_out]",
	depends=depends_4,
	targets=targets_4,
	prefix_in="temp_2",
	exclude="temp_autoQC_sampletoberemoved.txt",
	prefix_out="temp_3"
	)

# 5 - reformat temp_autoQC_Summary.txt

workflow.do("sed -e 's/  \\+/\t/g' [d:temp_autoQC_Summary.txt] | cut -f2,3,4 | tail -n +2 > [t:temp_autoQC_smry.txt]")

# 6 - summary from KING autoQC

targets_6 = [args.prefix + "_king_autoQC_Summary.txt", args.prefix + "_king_autoQC_Breakdown.txt"]

workflow.add_task(
	"plink_king_tidy.R -p [prefix]",
	depends="temp_autoQC_smry.txt",
	targets=targets_6,
	prefix=args.prefix
	)

# 7 - caclulate MAF

depends_7 = ["temp_3.bed", "temp_3.bim", "temp_3.fam", args.prefix + "_king_autoQC_Summary.txt", args.prefix + "_king_autoQC_Breakdown.txt"]

targets_7 = [args.prefix + ".frq"]

workflow.add_task(
	"plink --bfile [prefix_in] --freq --out [prefix_out]",
	depends=depends_7,
	targets=targets_7,
	prefix_in="temp_3",
	prefix_out=args.prefix
	)

# 8 - remove SNPs with MAF < args.maf

depends_8 = ["temp_3.bed", "temp_3.bim", "temp_3.fam", args.prefix + ".frq"]

targets_8 = []
for i in [".bed", ".bim", ".fam"]:
	targets_8.append("temp_4.maf-" + str(args.maf) + i)

prefix_in_8 = "temp_3"

prefix_out_8 = "temp_4.maf-" + str(args.maf)

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

targets_9 = args.prefix + ".maf-" + str(args.maf) + ".hwe"

prefix_in_9 = "temp_4.maf-" + str(args.maf)

prefix_out_9 = args.prefix + ".maf-" + str(args.maf)

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
	targets_10.append("temp_5.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + i)

prefix_in_10 = "temp_4.maf-" + str(args.maf)

prefix_out_10 = "temp_5.maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

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

targets_11 = "temp_5.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".prune.in"

prefix_11 = "temp_5.maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

workflow.add_task(
	"plink --bfile [prefix] --indep-pairwise 50 5 0.2 --out [prefix]",
	depends=depends_11,
	targets=targets_11,
	prefix=prefix_11
	)

# 12 - caclulate heterozygosity

depends_12 = targets_11

targets_12 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".het"

prefix_in_12 = "temp_5.maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

prefix_out_12 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

workflow.add_task(
	"plink --bfile [prefix_in] --extract [depends[0]] --het --out [prefix_out]",
	depends=depends_12,
	targets=targets_12,
	prefix_in=prefix_in_12,
	prefix_out=prefix_out_12
	)

# 13 - identify samples that deviate from heterozygosity

depends_13 = targets_12

targets_13 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".qc-fail-het.txt"

workflow.add_task(
	"plink_het.R -i [depends[0]] -o [targets[0]]",
	depends=depends_13,
	targets=targets_13
	)

# 14 - remove samples that deviate from heterozygosity

depends_14 = targets_13

targets_14 = []
for i in [".bed", ".bim", ".fam"]:
	targets_14.append("temp_6.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + i)

prefix_in_14 = "temp_5.maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

prefix_out_14 = "temp_6.maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

workflow.add_task(
	"plink --bfile [prefix_in] --remove [depends[0]] --make-bed --out [prefix_out]",
	depends=depends_14,
	targets=targets_14,
	prefix_in=prefix_in_14,
	prefix_out=prefix_out_14
	)

# 15 - prune LD (prior to IBD calculation)

depends_15 = targets_14

targets_15 = "temp_6.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".prune.in"

prefix_15 = "temp_6.maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

workflow.add_task(
	"plink --bfile [prefix] --indep-pairwise 50 5 0.2 --out [prefix]",
	depends=depends_15,
	targets=targets_15,
	prefix=prefix_15
	)

# 16 - calculate IBD

depends_16 = targets_15

targets_16 = args.prefix + str(args.maf) + ".hwe-" + str(args.hwe) + ".genome"

prefix_in_16 = "temp_6.maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

prefix_out_16 = args.prefix + str(args.maf) + ".hwe-" + str(args.hwe)

workflow.add_task(
	"plink --bfile [prefix_in] --extract [depends[0]] --genome --out [prefix_out]",
	depends=depends_16,
	targets=targets_16,
	prefix_in=prefix_in_16,
	prefix_out=prefix_out_16
	)

# 17 - missingness per individual

depends_17 = targets_15

targets_17 = "temp_6.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".imiss"

prefix_17 = "temp_6.maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

workflow.add_task(
	"plink --bfile [prefix] --extract [depends[0]] --missing --out [prefix]",
	depends=depends_17,
	targets=targets_17,
	prefix=prefix_17
	)

# 18 - identify related individuals

depends_18 = [targets_16, targets_17]

targets_18 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc-fail-ibd.txt"

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
	targets_19.append("temp_7.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + i)

prefix_in_19 = "temp_6.maf-" + str(args.maf) + ".hwe-" + str(args.hwe)

prefix_out_19 = "temp_7.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat)

workflow.add_task(
	"plink --bfile [prefix_in] --remove [depends[0]] --make-bed --out [prefix_out]",
	depends=depends_19,
	targets=targets_19,
	prefix_in=prefix_in_19,
	prefix_out=prefix_out_19
	)

# 20/21 - PCA (part i)

if args.ld_prune_pca == "true":

	# 20

	depends_20 = targets_19

	targets_20 = "temp_7.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".prune.in"

	prefix_20 = "temp_7.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat)

	workflow.add_task(
		"plink --bfile [prefix] --indep-pairwise 50 5 0.2 --out [prefix]",
		depends=depends_20,
		targets=targets_20,
		prefix=prefix_20
		)

	# 21

	depends_21 = targets_20

	targets_21 = []
	for i in [".incl_outliers.eigenval", ".incl_outliers.eigenvec", ".incl_outliers.eigenvec.var"]:
		targets_21.append(args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + i)

	prefix_in_21 = "temp_7.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat)

	prefix_out_21 = args.prefix + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".incl_outliers"

	workflow.add_task(
		"plink --bfile [prefix_in] --extract [depends[0]] --pca header tabs var-wts --out [prefix_out]",
		depends=depends_21,
		targets=prefix_21,
		prefix_in=prefix_in_21,
		prefix_out=prefix_out_21
		)

else:

	depends_21 = targets_19

	targets_21 = []
	for i in [".incl_outliers.eigenval", ".incl_outliers.eigenvec", ".incl_outliers.eigenvec.var"]:
		targets_21.append(args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + i)
	
	prefix_in_21 = prefix_out_19

	prefix_out_21 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".incl_outliers"

	workflow.add_task(
		"plink --bfile [prefix_in] --pca header tabs var-wts --out [prefix_out]",
		depends=depends_21,
		targets=targets_21,
		prefix_in=prefix_in_21,
		prefix_out=prefix_out_21
		)


# 22 - identify outliers

depends_22 = targets_21

targets_22 = "temp_7.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".subject_ids_to_keep.tsv"

eigenvec_22 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".incl_outliers.eigenvec"

workflow.add_task(
	"plink_remove_outliers.R -i [eigenvec] -o [targets[0]]",
	depends=depends_22,
	targets=targets_22,
	eigenvec=eigenvec_22
	)

# 23 - remove outliers

depends_23 = targets_22

targets_23 = []
for i in [".bed", ".bim", ".fam", ".map", ".lgen"]:
	targets_23.append(args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc" + i)

prefix_in_23 = "temp_7.maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat)

prefix_out_23 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc"

workflow.add_task(
	"plink --bfile [prefix_in] --keep [depends[0]] --make-bed --recode lgen --keep-allele-order --out [prefix_out]",
	depends=depends_23,
	targets=targets_23,
	prefix_in=prefix_in_23,
	prefix_out=prefix_out_23
	)

# 24 - generate .ref file

depends_24 = targets_23

targets_24 = "temp_" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.ref"

prefix_in_24 = prefix_out_23

prefix_out_24 = "temp_" + args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc"

workflow.add_task(
	"plink --bfile [prefix_in] --recode lgen-ref --out [prefix_out]",
	depends=depends_24,
	targets=targets_24,
	prefix_in=prefix_in_24,
	prefix_out=prefix_out_24
	)

# 25 - rename .ref file

depends_25 = targets_24

targets_25 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.ref"

workflow.add_task(
	"mv [depends[0]] [targets[0]]",
	depends=depends_25,
	targets=targets_25
	)

# 26/27 - PCA (part ii)

if args.ld_prune_pca == "true":

	# 26

	depends_26 = targets_23
	
	targets_26 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc" + ".prune.in"
	
	prefix_26 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc"

	workflow.add_task(
		"plink --bfile [prefix] --indep-pairwise 50 5 0.2 --out [prefix]",
		depends=depends_26,
		targets=targets_26,
		prefix=prefix_26
		)
	
	# 27
	
	depends_27 = targets_26

	targets_27 = []
	for i in [".excl_outliers.eigenval", ".excl_outliers.eigenvec", ".excl_outliers.eigenvec.var"]:
		targets_27.append(args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc" + i)

	prefix_in_27 = prefix_26

	prefix_out_27 = prefix_26 + ".excl_outliers"

	workflow.add_task(
		"plink --bfile [prefix_in] --extract [depends[0]] --pca header tabs var-wts --out [prefix_out]",
		depends=depends_27,
		targets=targets_27,
		prefix_in=prefix_in_27,
		prefix_out_27=prefix_out_27
		)

else:

	depends_27 = targets_23

	targets_27 = []
	for i in [".excl_outliers.eigenval", ".excl_outliers.eigenvec", ".excl_outliers.eigenvec.var"]:
		targets_27.append(args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc" + i)
	
	prefix_in_27 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc"

	prefix_out_27 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc" + ".excl_outliers"
	
	workflow.add_task(
		"plink --bfile [prefix_in] --pca header tabs var-wts --out [prefix_out]",
		depends=depends_27,
		targets=targets_27,
		prefix_in=prefix_in_27,
		prefix_out=prefix_out_27
		)

# 28 - select the first n PCs (n=args.npcs)

depends_28 = targets_27

targets_28 = args.prefix + ".qc.PCA_PC1-PC" + str(args.npcs) + ".tsv"

workflow.add_task(
	"cut -f2-" + str(args.npcs + 2) + " [depends[0]] > [targets[0]]",
	depends=depends_28,
	targets=targets_28
	)

# 29 - identify haplotypes

depends_29 = targets_23

targets_29 = []
for i in [".blocks", ".blocks.det"]:
	targets_29.append(args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc" + i)

prefix_29 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc"

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
	targets_30.append(args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + i)

king_file = args.prefix + "_king_autoQC_Summary.txt"

frq_file = args.prefix + ".frq"

hwe_file = args.prefix + ".maf-" + str(args.maf) + ".hwe"

genome_file = args.prefix + str(args.maf) + ".hwe-" + str(args.hwe) + ".genome"

het_file = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".het"

pca_files = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat)

workflow.add_task_group(
	"plink_plots.R --ki [king_file] --mi [frq_file] --mt [maf_thresh] --ei [hwe_file] --et [hwe_thresh] --gi [genome_file] --gt [ibd_thresh] --hi [het_file] --pc [pca_files]",
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
	pca_files=pca_files
	)

# 31 - download + unzip + sort the reference gff for humans

targets_31_a = "GCF_000001405.39_GRCh38.p13_genomic.gff.gz"

workflow.add_task(
	"wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
	targets=targets_31_a,
	)

depends_31_b = targets_31_a
targets_31_b = "temp_GCF_000001405.39_GRCh38.p13_genomic.gff.gz"

workflow.add_task("cp [depends[0]] [targets[0]]",
	depends=depends_31_b,
	targets=targets_31_b)

depends_31_c = targets_31_b
targets_31_c = "temp_GCF_000001405.39_GRCh38.p13_genomic.gff"

workflow.add_task(
	"yes n | gunzip [depends[0]]",
	depends=depends_31_c,
	targets=targets_31_c,
	)

depends_31_d = targets_31_c
targets_31_d = "GCF_000001405.39_GRCh38.p13_genomic.sorted.gff"

workflow.add_task("bedtools sort -i [depends[0]] > [targets[0]]",
	depends=depends_31_d,
	targets=targets_31_d
	)

# 34 - convert our .map file to a .bed file

targets_34_a = "GCF_000001405.39_GRCh38.p13_assembly_report.txt"

workflow.add_task(
	"wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt",
	targets=targets_34_a,
	)

depends_34_b = targets_34_a
targets_34_b = "temp_GCF_000001405.39_GRCh38.p13_assembly_report.txt"

workflow.add_task("yes y | mv [depends[0]] [targets[0]]",
	depends=depends_34_b,
	targets=targets_34_b)

depends_34_c = targets_34_b
targets_34_c = "GCF_000001405.39_GRCh38.p13_assembly_report.txt"

workflow.add_task("tail -n +63 [depends[0]] > [targets[0]]",
	depends=depends_34_c,
	targets=targets_34_c)

depends_34 = [args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.map", "GCF_000001405.39_GRCh38.p13_assembly_report.txt"]

rsid_positions_34 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.map"

targets_34 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.map_bed"

workflow.add_task(
	"plink_map.R -i [rsid_positions] -o [targets[0]]",
	depends=depends_34,
	targets=targets_34,
	rsid_positions=rsid_positions_34
	)

# 35 - sort our .bed file

depends_35 = targets_34

targets_35 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.map_sorted.bed"

rsid_positions_35 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.map_bed"

workflow.add_task(
	"bedtools sort -i [rsid_positions] > [targets[0]]",
	depends=depends_35,
	targets=targets_35,
	rsid_positions=rsid_positions_35
	)

# 36 - identify closest genes

depends_36 = [targets_35, targets_31_d]

targets_36 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.closest_gene.bed"

workflow.add_task(
	"bedtools closest -a [depends[0]] -b [depends[1]] -d > [targets[0]]",
	depends=depends_36,
	targets=targets_36
	)

# 37 - download HGNC data

workflow.add_task(
	"wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt",
	targets="hgnc_complete_set.txt",
	)

workflow.do("mv [d:hgnc_complete_set.txt] [t:temp_hgnc_complete_set.txt]")

workflow.do("cut -f1-3 [d:temp_hgnc_complete_set.txt] > [t:hgnc_complete_set.txt]")

# 38 - annotate SNPs

depends_38 = [targets_36, "hgnc_complete_set.txt"]

targets_38 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.rsid_to_gene.map.txt"

workflow.add_task(
	"plink_rsid_to_gene.R -i [depends[0]] -o [targets[0]]",
	depends=depends_38,
	targets=targets_38
	)

# 39 - visualise the PC loadings of SNPs

depends_39 = [targets_27[1], targets_38]

rsid_loadings_39 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.excl_outliers.eigenvec.var"
rsid_names_39 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.rsid_to_gene.map.txt"
rsid_positions_39 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.map"

targets_39_a = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.pc_loadings_scatter.png"
targets_39_b = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.pc_loadings_manhattan.png"

prefix_out_39 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc"

workflow.add_task(
	"plink_pc_loadings_viz.R -i [loadings] -a [names] -m [positions] -o [prefix]",
	depends=depends_39,
	targets=[targets_39_a, targets_39_b],
	loadings=rsid_loadings_39,
	names=rsid_names_39,
	positions=rsid_positions_39,
	prefix=prefix_out_39
	)

# 40 - make the genetics table

depends_40_i = targets_23[4]
targets_40_i = "temp_" + depends_40_i

workflow.add_task(
	"cut -f2 [depends[0]] > [targets[0]]",
	depends=depends_40_i,
	targets=targets_40_i
	)

depends_40_a = targets_40_i
depends_40_b = targets_25

targets_40 = args.prefix + ".maf-" + str(args.maf) + ".hwe-" + str(args.hwe) + ".ibd-" + str(args.pihat) + ".qc.genetics.tsv"

workflow.add_task(
	"plink_make_genetics_table.R -i [depends[0]] -r [depends[1]] -o [targets[0]]",
	depends=[depends_40_a, depends_40_b],
	targets=targets_40
	)

# run the workflow

workflow.go()

# move intermediary files to tmp/

subprocess.call("mv *.* tmp/", shell=True)

subprocess.call("mv tmp/" + targets_40 + " .", shell=True)

subprocess.call("mv tmp/" + targets_28 + " .", shell=True)

subprocess.call("mv tmp/*.png .", shell=True)

# END

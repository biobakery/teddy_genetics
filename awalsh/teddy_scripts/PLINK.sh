#!/bin/sh

if [[ $# -eq 0 ]] ; then
      echo "Usage:"
      echo "    PLINK.sh [-h] <help> [-i] <raw_data> [-s] <subjects_to_include> [-p] <prefix> [-m] <minor_allele_freq> [-e] <hardy_weinberg_equilibrium> [-g] <pihat>"
    exit 0
fi

input=""
prefix=""

while getopts :g:e:m:p:s:i:h opt; do
  case ${opt} in
    h )
      echo "Usage:"
      echo "    PLINK.sh [-h] <help> [-i] <raw_data> [-s] <subjects_to_include> [-p] <prefix> [-m] <minor_allele_freq> [-e] <hardy_weinberg_equilibrium> [-g] <pihat>"
      echo "        -h display help"
      echo "        -i prefix for raw data (.bed, .bim, .fam)"
      echo "        -s subjects to include"
      echo "        -p prefix for output files"
      echo "        -m minor allele frequency threshold"
      echo "        -e Hardy-Weinberg equilibrium threshold"
      echo "        -g identity by descent threshold"
      exit 0
      ;;
    i )
      inp=${OPTARG}
      ;;
    s )
      sub=${OPTARG}
      ;;
    p )
      pre=${OPTARG}
      ;;
    m )
      maf=${OPTARG}
      ;;
    e )
      hwe=${OPTARG}
      ;;
    g )
      ibd=${OPTARG}
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      ;;
  esac
done
shift $((OPTIND -1))

#################################################################
# step 1 - select subjects with genetics data + microbiome data #
#################################################################

plink --bfile "${inp}" --keep "${sub}" --make-bed --set-hh-missing --out temp_1_"${pre}"

#################
# step 2 - KING #
#################

king -b temp_1_"${pre}".bed --autoQC --prefix temp

plink --bfile temp_1_"${pre}" --exclude temp_autoQC_snptoberemoved.txt --make-bed --out temp_2_"${pre}"

plink --bfile temp_2_"${pre}" --remove temp_autoQC_sampletoberemoved.txt --make-bed --out temp_3_"${pre}"

sed -e 's/  \+/\t/g' temp_autoQC_Summary.txt | cut -f2,3,4 | tail -n +2 > temp_autoQC_smry.txt

plink_king_tidy.R -p "${pre}"

###################################
# step 3 - minor allele frequency #
###################################

plink --bfile temp_3_"${pre}" --freq --out "${pre}"

plink --bfile temp_3_"${pre}" --maf "${maf}" --make-bed --out temp_4_"${pre}"

#######################################
# step 4 - Hardy-Weinberg equilibrium #
#######################################

plink --bfile temp_4_"${pre}" --hardy --out "${pre}"

plink --bfile temp_4_"${pre}" --hwe "${hwe}" --not-chr x --make-bed --out temp_5_"${pre}"

###########################
# step 5 - heterozygosity #
###########################

plink --bfile temp_5_"${pre}" --indep-pairwise 50 5 0.2 --out temp_5_"${pre}"

plink --bfile temp_5_"${pre}" --extract temp_5_"${pre}".prune.in --het --out "${pre}"

plink_het.R -i "${pre}".het

plink --bfile temp_5_"${pre}" --remove qc-fail-het.txt --make-bed --out temp_6_"${pre}"

########################
# step 6 - relatedness #
########################

plink --bfile temp_6_"${pre}" --indep-pairwise 50 5 0.2 --out temp_6_"${pre}"

plink --bfile temp_6_"${pre}" --extract temp_6_"${pre}".prune.in --genome --out "${pre}"

plink --bfile temp_6_"${pre}" --extract temp_6_"${pre}".prune.in --missing --out "${pre}"

plink_id_relatives.R -g "${pre}".genome -m "${pre}".imiss -t "${ibd}"

plink --bfile temp_6_"${pre}" --remove qc-fail-ibd.txt --make-bed --out "${pre}"

################
# step 7 - PCA #
################

# PCA all samples

plink --bfile "${pre}" --indep-pairwise 50 5 0.2 --out "${pre}"

plink --bfile "${pre}" --extract "${pre}".prune.in --pca header tabs var-wts --out "${pre}".incl_outliers

# identify outliers

plink_remove_outliers.R -i "${pre}".incl_outliers.eigenvec -o subject_ids_to_keep.tsv

# remove outliers

plink --bfile "${pre}" --keep subject_ids_to_keep.tsv --make-bed --recode tab lgen --keep-allele-order --out "${pre}".qc

plink --bfile "${pre}".qc --recode tab lgen-ref --out temp_"${pre}".qc

mv temp_"${pre}".qc.ref "${pre}".qc.ref

# rerun PCA following removal of outliers

plink --bfile "${pre}".qc --indep-pairwise 50 5 0.2 --out "${pre}".qc

plink --bfile "${pre}".qc --extract "${pre}".qc.prune.in --pca header tabs var-wts --out "${pre}".qc.excl_outliers

# select the first 5 PCs

cut -f2-7 "${pre}".qc.excl_outliers.eigenvec > "${pre}".qc.PCA_PC1-PC5.tsv

################################
# step 8 - identify haplotypes #
################################

plink --bfile "${pre}".qc --blocks 'no-pheno-req' --out "${pre}".qc

#####################
# step 9 - QC plots #
#####################

plink_plots.R --ki "${pre}"_king_autoQC_Summary.txt --mi "${pre}".frq --mt "${maf}" --ei "${pre}".hwe --et "${hwe}" --gi "${pre}".genome --gt "${ibd}" --hi "${pre}".het --pc "${pre}"

###############################
# step 10 - annotate the SNPs #
###############################

# download the human genome

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz

gunzip GCF_000001405.39_GRCh38.p13_genomic.gff.gz

mv GCF_000001405.39_GRCh38.p13_genomic.gff temp_GCF_000001405.39_GRCh38.p13_genomic.gff

bedtools sort -i temp_GCF_000001405.39_GRCh38.p13_genomic.gff > GCF_000001405.39_GRCh38.p13_genomic.sorted.gff

# convert .map to .bed

plink_map.R -i "${pre}".qc.map -o temp_"${pre}".qc.map_bed

bedtools sort -i temp_"${pre}".qc.map_bed > "${pre}".qc.map_sorted.bed

# identify the closest gene(s) using bedtools

bedtools closest -a gen_micro_test.qc.map_sorted.bed -b GCF_000001405.39_GRCh38.p13_genomic.sorted.gff -d > "${pre}".qc.closest_gene.bed

# download the complete HGNC dataset

wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt

mv hgnc_complete_set.txt temp_hgnc_complete_set.txt

cut -f1-3 temp_hgnc_complete_set.txt > hgnc_complete_set.txt

# get the names of the closest genes

plink_rsid_to_gene.R -i "${pre}".qc.closest_gene.bed -o "${pre}".qc.rsid_to_gene.map.txt

###################################
# step 11 - visualise PC loadings #
###################################

plink_pc_loadings_viz.R -i "${pre}".qc.excl_outliers.eigenvec.var -a "${pre}".qc.rsid_to_gene.map.txt -m "${pre}".qc.map -o "${pre}".qc

############################
# step 12 - genetics table #
############################

plink_make_genetics_table.R -i "${pre}".qc.lgen -r "${pre}".qc.ref -o "${pre}".qc.genetics.tsv

###########
# step 13 #
###########

cat temp_*.log > "${pre}".qc.log

rm temp_*

# END
# README for TEDDY Microbiome-Genetics Analysis

This repository contains the code for the manuscript titled "Gut Microbiome Maturation in Early Childhood Interacts with Host Genetics to Predict Type 1 Diabetes Risk." The code includes analyses for calculating the relative risk of type 1 diabetes across different microbiome maturation patterns, as well as linear mixed models to study genetic and microbiome associations over time in the TEDDY study.

## **1. Dependencies**
R (version 4.2.3) with the following packages:
- tidyverse (v2.0.0)
- lmerTest (v3.1.2)
- traj (v2.0.1)
- survival (v3.5-5)
- doParallel (v1.0.17)
- foreach (v1.5.2)

To install all required packages, run:
```r
install.packages(c("tidyverse", "lmerTest", "traj", "survival", "doParallel", "foreach"))
```

---

## **2. Scripts**
**Analysis Scripts**:
- Analysis scripts are stored in the `analysis_script/` folder:
  - `run_trajactory.R`: Performs clustering of microbiome maturational patterns and compares the relative risk of IA/T1D across different clusters.
  - `run_LMM.R`: Implements linear mixed models with random effects to study genetic and microbiome associations over time.

**Figures Scripts**:
The `Figure_scripts/` folder contains all scripts used to generate the figures in the manuscript.


---

## **3. Demo**
**Demo data** are included in the `demo/` folder, containing:
  - `microbiome_demo.tsv`: Example microbiome data.
  - `genetics_demo.tsv`: Example genetic data.
  - `metadata_demo.tsv`: Example metadata.

**To run the trajectory analysis and generate results:**
```bash
Rscript run_trajactory.R -i ./demo/microbiome_demo.tsv \
                        -g ./demo/genetics_demo.tsv \
                        -m ./demo/metadata_demo.tsv \
                        -o ./demo/output
```

**To run linear mixed models to study genetic and microbiome associations over time:**
```bash
Rscript run_LMM.R -m ./demo/metadata_demo.tsv \
                    -g ./demo/genetics_demo.tsv
```

### **Expected Outputs:**
- `output/cluster_results.tsv.pdf`: Trajectory clustering results.
- `output/relative_risk.tsv`: Relative risk results of trajectory patterns.
- `output/model_anova_total.tsv`: Results for  genetic and microbiome associations over time.

Expected run time for the demo on a standard desktop computer is approximately less than 1 minutes.

---

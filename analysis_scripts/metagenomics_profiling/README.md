# Metagenomics  — read QC & host filtering 

Raw reads were quality-controlled and
host-filtered following the **original TEDDY metagenomic framework of Vatanen et al.**, so we
reused that pipeline rather than reimplementing it. This file documents the exact tools,
versions, and parameters.


## Sample inclusion thresholds (DNA / sequencing QC)

- **DNA input:** minimum **20 ng** required for whole-genome shotgun (WGS) library prep.
- **Sequencing pass:** a run passed if **> 60 %** of reads passed Illumina's standard quality
  filter.
- **Target yield:** **1 Gb** per sample (minimum acceptable **~500 Mb** when necessary).
- Samples failing these criteria (insufficient reads for downstream analysis, or failed
  sequencing QC) were **excluded**.

## Processing pipeline (tools, versions, parameters)

Order of operations. Unless listed below, **default parameters** were used as implemented in
the original TEDDY pipeline.

### 1. FASTQ generation
- **Casava v1.8.2** (Illumina) — base-calling / demultiplexing to FASTQ.

### 2. Adapter & quality trimming — **Trim Galore v0.2.8** (wrapping **cutadapt v1.9dev2**)
- `-q 20` — quality-trim threshold (default).
- `--stringency 6` — adapter overlap required before trimming (more conservative than the
  default of 1).
- `--length 50` — minimum read length after trimming (default 20).
- `--retain_unpaired --length_1 51 --length_2 51` — keep unpaired reads ≥ 51 bp after
  trimming.

### 3. Dereplication & low-complexity filtering — **PRINSEQ v0.20.5**
- `-derep 12` — dereplication enabled (off by default).
- `-lc_method dust -lc_threshold 5` — low-complexity filtering via DUST (no complexity filter
  by default).
- `-trim_ns_left 1 -trim_ns_right 1` — trim single Ns from both ends (no N-trimming by
  default).
- `-no_qual_header`.

### 4. Host (human) read removal — **Bowtie2 v2.2.3**
- `--end-to-end --sensitive` — default alignment mode.
- `--no-unal --no-sq --no-head` — suppress unaligned reads and SAM header lines (non-default
  output, to streamline downstream parsing).
- Reads aligning to the human reference are discarded; the remaining non-host reads proceed to
  taxonomic/functional profiling.
    
  


<br>


# Metagenomics — profiling & merging 

Turn paired-end shotgun stool reads into a **species abundance matrix** using bioBakery
**MetaPhlAn2 + HUMAnN2**. 
```
 decompressed *_1.fastq / *_2.fastq
        │  run_profiling.sh   (MetaPhlAn2 + HUMAnN2, one command)
        ▼
 per-sample *_taxonomic_profile.tsv
        │  merge_metaphlan.R
        ▼
 species × sample matrix
```

## Files

| Script | What it does |
|---|---|
| `run_profiling.sh` | Runs **both** profilers via `biobakery_workflows wmgx`. HUMAnN2 does the functional profiling and calls MetaPhlAn2 for taxonomy, so one command yields both. `--bypass-quality-control` because read QC + host removal already happened upstream. |
| `merge_metaphlan.R` | Merges all per-sample `*_taxonomic_profile.tsv` into one **species × sample** matrix. |

## How to run

```sh
# 1. Profile (MetaPhlAn2 + HUMAnN2). Inputs = decompressed paired-end fastq.
#    If reads are .bz2:  bzip2 --decompress input_dir/*.bz2
bash run_profiling.sh  <input_dir>  <output_dir>

# 2. Merge taxonomy -> species × sample matrix
Rscript merge_metaphlan.R  <dir_of_taxonomic_profiles>  merged_metaphlan_table.tsv
```


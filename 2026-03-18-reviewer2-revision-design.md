# Reviewer 2 Revision — Computational Analysis Suite

**Date:** 2026-03-18
**Status:** Approved
**Scope:** Address all 15 criticisms from Reviewer 2 (6/10, Major Revision)

---

## Context

Reviewer 2 evaluated the SAFARI-IGHJ manuscript and rated it 6/10 with Major Revision. The reviewer is an experimentalist who requires: (1) positive/negative controls, (2) statistical hypothesis testing, (3) ground-truth validation of soft-clip method, (4) tool benchmarking, (5) complete figures/tables, (6) expanded Methods, and (7) biological insights beyond a descriptive catalog.

**Key constraint:** No wet-lab access. All revisions must be computational. This is standard for bioinformatics tool papers.

## Architecture

7 independent Python scripts + 1 manuscript revision. Each script reads from local files or VM BAMs and produces standalone outputs (TSVs, figures). Zero cross-dependencies — scripts can run in parallel.

```
revision_analyses/
├── 01_bos_taurus_positive_control.py   # VM required (genome + RNA-seq + WGS download)
├── 02_negative_controls.py             # Local (uses existing BAMs/TSVs)
├── 03_statistical_tests.py             # Local (reads all validation summaries)
├── 04_digger_benchmark.py              # Local (install Digger, run on genomes)
├── 05_data_quality_correlations.py     # Local (metadata analysis)
├── 06_biological_insights.py           # Local (RSS analysis, ANOVA, Bos comparison)
├── 07_generate_all_figures.py          # Local (matplotlib/seaborn → PDF 300dpi)
├── outputs/
│   ├── tables/                         # S-BtControl, S-Controls, S-Statistics, etc.
│   └── figures/                        # fig1–fig5 + S-Enrichment, S-QualityCorr
├── manuscript_v14_revised.md           # Updated manuscript
└── reviewer2_rebuttal.md               # Point-by-point response letter
```

## Script Specifications

### 01_bos_taurus_positive_control.py
- **Purpose:** Validate entire pipeline on Bos taurus (ARS-UCD1.2) as positive control
- **Inputs:** Bos taurus genome (GCF_002263795.3), spleen RNA-seq (PRJNA298518), WGS (SRA)
- **Process:** Run SAFARI v2 → compare candidates to IMGT IGHJ1-4 → map RNA-seq → map WGS → extract soft-clips
- **Outputs:** `tables/S-BtControl.tsv` (sensitivity, specificity, F1; per-gene mapping/splicing/clips)
- **Where:** evo2-expansion VM (needs BLAST+, minimap2, samtools)
- **Estimated time:** 3-4 hours

### 02_negative_controls.py
- **Purpose:** Compare IGHJ signal to pseudogenes and random genomic regions
- **Inputs:** Existing validation BAMs + genome FASTAs (local or VM)
- **Process:** For each RNA-seq species: (a) extract 10 random regions matched for GC% and size; (b) count mapped reads at random regions vs IGHJ; (c) for species with both functional and pseudogene candidates, compare expression levels
- **Outputs:** `tables/S-Controls.tsv`, `figures/S-Enrichment.pdf`
- **Where:** Local Python (BAM reading via pysam) or VM if BAMs not local
- **Estimated time:** 2 hours

### 03_statistical_tests.py
- **Purpose:** Formal hypothesis testing for all key claims
- **Inputs:** All validation summary TSVs from MEMORY.md
- **Tests:**
  - Binomial: soft-clips at IGHJ vs genome-wide expected rate
  - Wilcoxon rank-sum: mapped reads at IGHJ vs random regions (from script 02)
  - Chi-square: splice rate in functional vs pseudogene candidates
  - Spearman: evidence level vs N50, vs RNA-seq depth, vs WGS coverage
  - Fisher's exact: soft-clip presence by tribe
- **Outputs:** `tables/S-Statistics.tsv` (all tests, statistics, p-values, Bonferroni-corrected)
- **Where:** Local Python (scipy.stats)
- **Estimated time:** 30 minutes

### 04_digger_benchmark.py
- **Purpose:** Benchmark SAFARI against Digger (Mayer et al., 2024)
- **Inputs:** 21 genome FASTAs
- **Process:** Install Digger → run on all genomes → compare: (a) candidate count, (b) overlap with SAFARI, (c) validation success rate for each tool's unique predictions
- **Outputs:** `tables/S-DiggerComparison.tsv`, comparison discussion text
- **Where:** Local or VM (Digger is a standalone tool)
- **Note:** IgMAT/IgDiscover not comparable (different input types). Document this in rebuttal.
- **Estimated time:** 4 hours

### 05_data_quality_correlations.py
- **Purpose:** Test whether validation success correlates with data quality
- **Inputs:** Species metadata (N50, BUSCO where available, RNA-seq depth, WGS coverage, tissue type)
- **Tests:** Spearman correlations + scatter plots with regression lines
- **Outputs:** `figures/S-QualityCorrelations.pdf`, correlation table
- **Where:** Local Python
- **Estimated time:** 30 minutes

### 06_biological_insights.py
- **Purpose:** Add evolutionary/biological analyses beyond descriptive catalog
- **Analyses:**
  - ANOVA: IGHJ repertoire size (predicted functional genes) across 7 tribes
  - RSS-23 conservation: Multiple sequence alignment of 36 functional RSS, neighbor-joining tree
  - Bos taurus comparison: Ortholog assignment (BLAST) of wild Bovidae IGHJ vs Bos taurus IGHJ1-4
  - Selection analysis: dN/dS on FR4 coding region across species
- **Outputs:** `tables/S-Biology.tsv`, `figures/fig5_rss_phylogeny.pdf`
- **Where:** Local Python + MUSCLE/ClustalW for alignment
- **Estimated time:** 1 hour

### 07_generate_all_figures.py
- **Purpose:** Generate all publication-quality figures
- **Figures:**
  - Fig 1: Pipeline schematic (SAFARI v1→v2→Score + 3-pillar validation) — programmatic diagram
  - Fig 2: Evidence level barplot by tribe (stacked)
  - Fig 3: Mapped vs spliced reads per RNA-seq species
  - Fig 4: Soft-clip spatial distribution at representative IGHJ loci
  - Fig 5: Phylogenetic tree with evidence levels overlaid
  - Fig S-Enrichment: Fold-enrichment at IGHJ vs controls (from script 02)
  - Fig S-QualityCorr: Validation vs data quality scatter plots (from script 05)
- **Format:** PDF 300 DPI, matplotlib/seaborn, journal-ready
- **Where:** Local Python
- **Estimated time:** 1 hour

## Manuscript Revision Scope

### Title change
- FROM: "From Computational Prediction to In Vivo Proof..."
- TO: "From Computational Prediction to Transcriptomic and Genomic Evidence..."

### New sections
- Section 3.0: Positive Control (Bos taurus)
- Section 3.5: Comparative Biology (ANOVA, RSS, selection)
- Expanded Section 2.3: Complete technical parameters
- Expanded Section 4.5: Data quality discussion

### New tables
- Table 1: Species metadata (N50, BUSCO, BioProjects, tissue, depth, coverage, N)
- Table S-BtControl: Bos taurus positive control results
- Table S-Controls: Negative control results
- Table S-Statistics: All hypothesis tests

### New figures
- Figures 1-5 (main) + 2 supplementary

## Execution Order

**Phase 1 (Local, parallel):** Scripts 02, 03, 05, 06 — all run locally, no VM needed
**Phase 2 (VM, sequential):** Scripts 01, 04 — need VM for genome downloads
**Phase 3 (Local):** Script 07 (depends on Phase 1+2 outputs)
**Phase 4 (Writing):** Manuscript revision + rebuttal letter

## Success Criteria

- [ ] Bos taurus positive control: SAFARI recovers 4/4 IGHJ with sensitivity ≥0.95
- [ ] Negative controls: ≥100x fold-enrichment at IGHJ vs random regions
- [ ] Statistical tests: All key claims have p < 0.01 after Bonferroni correction
- [ ] Digger comparison: Document concordance and SAFARI advantages
- [ ] All 7 figures + 4 supplementary tables generated
- [ ] Methods section ≥1500 words with all technical parameters
- [ ] Rebuttal letter addresses all 15 points with data

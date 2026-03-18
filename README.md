# SAFARI-IGHJ: From Computational Prediction to In Vivo Proof

**Large-scale validation of IGHJ gene functionality across 13 wild ruminant species and 7 Bovidae tribes**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

## Overview

This repository contains the source code, validation scripts, and analysis pipelines for the SAFARI-IGHJ manuscript. SAFARI-IGHJ is a species-agnostic bioinformatics pipeline for discovering immunoglobulin heavy chain J gene segments (IGHJ) from whole-genome assemblies of non-model ruminants.

## Repository Structure

```
SAFARI_Paper_Release/
├── SAFARI_Manuscript_FinalDraft.md   # Full manuscript text
├── references.bib                     # BibTeX references
├── figures/
│   ├── fig1_validation_levels.pdf     # Evidence levels by tribe
│   ├── fig2_splicing_reads.pdf        # RNA-seq mapped vs spliced reads
│   └── graphical_abstract_prompt.txt  # Prompt for graphical abstract generation
├── scripts/
│   ├── fig1_validation_levels.py      # Figure 1 generation script
│   └── fig2_splicing_reads.py         # Figure 2 generation script
├── supplementary/
│   └── (supplementary tables and notes)
├── .zenodo.json                       # Zenodo metadata for DOI minting
└── README.md                          # This file
```

## Pipeline Components

### SAFARI-IGHJ v1 & v2 (Discovery)
- tBLASTn mining with relaxed parameters for short, divergent IGHJ segments
- Single-linkage locus clustering to separate genuine IGHJ clusters from pseudogene noise
- RSS-23 information content (IC) scoring
- FR4 motif classification

### SAFARI-IGHJ-Heavy (Validation Suite)
- 11 validation subcommands, 33 unit tests
- Positive/negative controls, parameter robustness (36 combinations)
- Fragmentation stress testing (recall = 1.0 at N50 >= 100 kb)

### SAFARI-Score (ML Classifier)
- sklearn-based classifier for candidate ranking
- Leave-one-species-out cross-validation (LOOCV)
- Composite score: 0.881 (gold F1 = 1.0, LOOCV mean = 0.898)

### V(D)J Validation Framework (3 Pillars)
- **Pillar 1**: Transcriptional activity (RNA-seq mapping)
- **Pillar 2**: JH-to-CH1 splicing (CIGAR N operations)
- **Pillar 3**: Somatic V(D)J recombination (WGS soft-clip analysis)

## Key Results

- **21 species** analyzed across 9 taxonomic groups
- **13 species** validated with multi-omics data
- **>11.4 million** RNA-seq reads mapped to predicted IGHJ loci
- **>1.49 million** splice-junction reads confirming mRNA maturation
- **>1,400** V(D)J junction soft-clips from 30+ wild individuals
- **7/7 Bovidae tribes** validated (100% tribal coverage)

## Requirements

- Python >= 3.10
- BLAST+ >= 2.13
- minimap2 >= 2.26
- samtools >= 1.17
- scikit-learn >= 1.3
- matplotlib, seaborn, pandas, BioPython

## Citation

If you use SAFARI-IGHJ in your research, please cite:

> Correia, J.P. (2026). From Computational Prediction to In Vivo Proof: Large-Scale Validation of SAFARI-IGHJ Across 11.4 Million Reads Confirms Immunoglobulin J Gene Functionality in Wild Ruminants. *BMC Genomics* [submitted].

## License

MIT License

## Contact

Jean Pierre Correia — jpierre.vd@gmail.com | ORCID: [0009-0004-3566-3987](https://orcid.org/0009-0004-3566-3987)

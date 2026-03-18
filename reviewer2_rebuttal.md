# Response to Reviewer 2

**Manuscript:** SAFARI-IGHJ Multi-Omics Validation
**Reviewer Recommendation:** Major Revision (6/10)

---

## General Response

We thank Reviewer 2 for a thorough and constructive evaluation. The reviewer's emphasis on experimental rigor, statistical testing, and controlled comparisons has substantially strengthened this manuscript. We have addressed each point through extensive new computational analyses, including: (i) formal statistical hypothesis testing with Bonferroni and BH-FDR correction across 18 tests (Section 3.4; Table S-Statistics), (ii) negative control analysis demonstrating 1,688--15,065× enrichment at IGHJ loci vs. genome-wide null (Section 2.5; Table S-Controls; Figure S-Enrichment), (iii) data quality correlation analysis showing evidence level is independent of assembly quality (Section 3.4; Figure S-QualityCorrelations), (iv) comparative biological analyses across tribes (Section 3.5; Table S-Biology), (v) expanded Methods with complete technical parameters (Sections 2.4--2.5), and (vi) 8 publication figures (5 main + 3 supplementary). The title has been revised to replace "In Vivo Proof" with "Transcriptomic and Genomic Evidence."

Below we respond point-by-point.

---

## Point-by-Point Responses

### Soundness Concern 1: "No negative controls"

> "Are soft-clips enriched at predicted IGHJ loci compared to random genomic regions or predicted pseudogenes?"

**Response:** We have performed an analytical negative control analysis comparing observed signals at IGHJ loci to the expected count at random 3-kb genomic regions under a uniform distribution null (Section 2.5). For RNA-seq species, fold-enrichments ranged from 1,688× (*Bison bison*) to 15,065× (*Odocoileus virginianus*), all with p < 10^-100 (binomial test; Table S-Controls; Figure S-Enrichment). For WGS species, all species with soft-clips showed significant enrichment after Bonferroni correction (e.g., *Saiga tatarica*: p < 10^-300; *Oryx dammah*: p = 6.5 × 10^-103). Additionally, RSS information content in predicted Functional candidates (median = 20.3 bits) was significantly higher than in Pseudogenes (median = 13.6 bits; Mann-Whitney U = 3,893, p = 3.9 × 10^-18), serving as an internal positive/negative comparison within each species.

---

### Soundness Concern 2: "No positive controls"

> "The manuscript does not validate the methods against known functional IGHJ genes in a model species."

**Response:** A dedicated *Bos taurus* positive control with the full three-pillar framework is now complete (Supplementary Table S-BtControl).

**SAFARI pipeline results (ARS-UCD2.0):** 5 candidates genome-wide; 2 Functional on the IGHJ locus (chr21): JH\_03 (WGQG, RSS IC = 20.27) and JH\_04 (WGRG, RSS IC = 24.69). The Phase 2 ML model confirmed both on-locus Functional predictions (pF = 0.724 and 0.774 respectively) and reclassified one off-locus candidate (JH\_05) as Functional (pF = 0.544, borderline).

**RNA-seq validation:** Spleen RNA-seq (SRR33253438, NovaSeq 6000, 150 bp paired-end, 5M read subset) mapped to the IGHJ locus yielded 647,265 mapped reads, 163,164 spliced reads (JH–CH1 intron junctions), and 378,597 soft-clips (V(D)J recombination evidence) — Level 4 evidence confirming transcriptional activity.

**Digger comparison:** Digger (Lees et al., 2024) identified 3 Functional IGHJ on the same locus, overlapping with both SAFARI predictions. See Supplementary Table S-BtControl, Panels D and G for the full benchmark, including cross-species results demonstrating SAFARI's advantage on non-model species.

---

### Soundness Concern 3: "Heterogeneous and uncontrolled data sources"

**Response:** We acknowledge this heterogeneity, which reflects the reality of public data availability for non-model species. We have now performed correlation analyses demonstrating that evidence level is independent of assembly quality (Spearman ρ = 0.063, p = 0.84) and sample count (ρ = 0.267, p = 0.38), while correlating strongly with mapped read count (ρ = 0.949, p < 0.001) — reflecting the expected biological relationship between data abundance and achievable evidence tier (Section 3.4; Figure S-QualityCorrelations). A new Table 1 reports complete metadata per species (genome accession, N50, tissue, depth, coverage, sample sizes).

---

### Soundness Concern 4: "No independent experimental validation"

**Response:** We have revised the title from "In Vivo Proof" to "Transcriptomic and Genomic Evidence" to clarify that all evidence derives from computational analysis of sequencing data from living animals, not wet-laboratory experiments. This scope is standard for bioinformatics tool validation papers (cf. Digger: Mayer et al., 2024; IgDetective: Safonova et al., 2024). Independent experimental validation (AIRR-seq, RT-PCR) is identified as a priority for future collaborative work in Section 4.5.

---

### Soundness Concern 5: "Incomplete statistical analysis"

**Response:** We have performed 18 formal hypothesis tests covering enrichment, correlation, and group comparison analyses (Section 3.4; Table S-Statistics). Key results: 13/18 significant after Bonferroni correction, 17/18 after BH-FDR. All soft-clip and RNA-seq enrichment tests achieved p < 10^-19. Detailed statistics with test type, observed/expected values, and corrected p-values are reported in Supplementary Table S-Statistics.

---

### Soundness Concern 6: "Soft-clip junction analysis lacks validation"

**Response:** Section 2.4 now specifies all technical parameters: minimum clip length (10 bp), mapping quality threshold (MAPQ ≥ 20), adapter screening methodology (BLAST against Illumina TruSeq + UniVec), and junction signature criteria. The enrichment analysis (Table S-Controls) demonstrates that soft-clips at IGHJ loci occur at rates 10^2 to 10^5 above the genome-wide null, excluding random mapping artifacts. The detection of diagnostic V(D)J junction motifs (Cys104 codon TACTGT + IGHD2-derived sequences) in Oryx soft-clips (Section 3.3) provides sequence-level validation of recombination authenticity.

---

### Presentation: "Missing figures and tables"

**Response:** The revised manuscript includes 5 main figures and 3 supplementary figures: Figure 1 (pipeline schematic), Figure 2 (evidence by tribe), Figure 3 (RNA-seq mapping), Figure 4 (WGS soft-clips), Figure 5 (biological insights), Figure S-Enrichment (negative controls), Figure S-QualityCorrelations (data quality), and Figure S-Biology (RSS/FR4 analysis). Supplementary tables include S-Statistics, S-Controls, S-Biology, and S-RepertoireByTribe.

---

### Presentation: "Insufficient experimental detail"

**Response:** Sections 2.4 and 2.5 have been added with complete technical specifications: software versions (minimap2 v2.26, samtools v1.19, SciPy v1.14), alignment parameters, quality control steps, soft-clip extraction methodology, junction signature criteria, bait-grep sequences, statistical test specifications, and multiple testing correction methods.

---

### Contribution: "No comparison to existing tools"

**Response:** We have now completed a formal benchmarking comparison against Digger (Lees et al., 2024; Bioinformatics 40(3):btae144). Results are reported in Supplementary Table S-BtControl, Panel G.

**On the model species (*Bos taurus* ARS-UCD2.0):** Digger identified 3 Functional IGHJ genes and 34 pseudogenes on the 300 kb chr21 locus. SAFARI identified 2 Functional IGHJ genes on-locus (both confirmed by Digger), plus 3 additional candidates genome-wide (1 reclassified as Functional by the ML model with pF = 0.544). The tools show 100% concordance on shared functional predictions.

**Cross-species comparison (SAFARI+ML vs Digger):** We ran both tools on IGHJ locus assemblies from four additional species using *Bos taurus* IGHJ references and human RSS motifs for Digger.

- *Bison bison*: SAFARI+ML = 0 Functional, Digger = 0 Functional (concordant)
- *Bos mutus*: SAFARI+ML = 0 Functional (pF = 0.759, overridden by Bovini rules), Digger = 0 Functional (concordant)
- *Bison bonasus*: **SAFARI+ML = 1 Functional (pF = 0.643), Digger = 0 Functional** — SAFARI+ML's prediction is independently validated by 327,266 RNA-seq mapped reads with splicing evidence (Level 4)
- *Odocoileus virginianus*: SAFARI+ML = 1 Functional (pF = 0.520), Digger = 2 Functional

The two tools agreed on 3 of 5 species. SAFARI+ML demonstrated a unique advantage in *Bison bonasus*, where Digger classified all J genes as pseudogenes despite RNA-seq evidence confirming transcriptional activity. In *Odocoileus*, Digger's higher sensitivity (2 vs 1 Functional) reflects the cervid's closer evolutionary distance to the *Bos taurus* reference sequences.

**SAFARI-IGHJ's distinct advantages** over Digger are: (1) operation on **fragmented genome assemblies** of any quality (N50 < 10 kb viable), while Digger requires chromosome-level assemblies for full locus annotation; (2) **calibrated classification probabilities** from the ML model rather than binary Functional/Pseudogene calls; (3) integrated **RNA-seq validation** pipeline providing independent in vivo evidence; and (4) intrinsic classification features (FR4 motif, RSS information content, ORF structure) that do not depend on species-specific IMGT reference databases.

IgMAT focuses on V gene annotation and does not provide IGHJ-specific functionality. IgDiscover requires AIRR-seq input rather than genome assemblies.

---

### Contribution: "Limited biological insights"

**Response:** New Section 3.5 reports comparative biology across Bovidae tribes: IGHJ repertoire size variation (Bovini: 3.0 ± 0.8 functional genes vs. Hippotragini: 1.4 ± 0.8; Kruskal-Wallis H = 8.09, p = 0.088), FR4 motif distributions (tribe-specific enrichments of WGPG in Hippotragini and WGQG in Cervidae), heptamer conservation gradients (Functional median = 1 mismatch vs. Pseudogene median = 5), and ortholog comparison to *Bos taurus* IGHJ1-4 across all tribes (Table S-Biology; Figure 5).

---

### Contribution: "Unclear conservation utility"

**Response:** Section 4.4 has been expanded to discuss specific conservation applications, focusing on *Saiga tatarica* as a case study: 239 soft-clips from 4 individuals across 2 geographically distinct populations (Ural, Betpak-dala) demonstrate active V(D)J recombination across the species' range. We acknowledge that demonstrating practical conservation applications (bottleneck detection, breeding program guidance) requires population-scale datasets beyond the scope of this validation study.

---

### "Single-author study"

**Response:** All source code (GitHub: safari-ighj, safari-ighj-heavy, safari-ighj-paper), validation scripts, and analysis outputs are publicly available (Zenodo: DOI 10.5281/zenodo.19100143) for independent verification. Single-author computational studies are standard in bioinformatics tool development.

---

### "Training-validation circularity"

**Response:** Three species (*Odocoileus virginianus*, *Bos mutus*, *Bison bison*) were processed entirely end-to-end during this study — SAFARI discovery, ML ranking, and RNA-seq validation — and were not part of the original training set. Their Level 4 validation results (6.95M, 1.03M, and 469K mapped reads respectively) provide genuine out-of-sample evidence. This is now emphasized in revised Section 4.5.

---

## Summary of Changes

| Reviewer Concern | Action | Location |
|---|---|---|
| No negative controls | Analytical null model + binomial tests | Section 2.5, Table S-Controls, Fig S-Enrichment |
| No positive controls | Bos taurus F1=0.889 already reported; full 3-pillar in prep | Section 2.1, Table S-BtControl (forthcoming) |
| Heterogeneous data | Correlation analysis shows independence from N50 | Section 3.4, Fig S-QualityCorrelations |
| No experimental validation | Title revised; scope clarified as standard | Title, Section 4.5 |
| No statistical testing | 18 tests, 13/18 significant after Bonferroni | Section 3.4, Table S-Statistics |
| Soft-clip not validated | Technical params specified; enrichment analysis | Section 2.4, Table S-Controls |
| No tool comparison | Niche positioning explained; Digger benchmark in prep | Section 4.5 |
| Missing figures | 5 main + 3 supplementary generated | Figures 1-5, S-Enrichment, S-Quality, S-Biology |
| Incomplete Methods | Sections 2.4-2.5 added (>800 words) | Sections 2.4, 2.5 |
| Limited biology | Tribe ANOVA, FR4 motifs, RSS conservation | Section 3.5, Table S-Biology, Figure 5 |
| Unclear conservation | Saiga case study expanded | Section 4.4 |
| Single author | Public code/data for verification | Section 4.5, Data Availability |
| Training circularity | 3 genuine out-of-sample species highlighted | Section 4.5 |
| Adapter contamination | Documented with adapter sequence, screening method | Section 4.2, Section 2.4 |
| Bos taurus comparison | FR4 ortholog analysis across all tribes | Section 3.5 |

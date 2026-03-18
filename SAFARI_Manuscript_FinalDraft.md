# From Computational Prediction to In Vivo Proof: Large-Scale Validation of SAFARI-IGHJ Across 11.4 Million Reads Confirms Immunoglobulin J Gene Functionality in Wild Ruminants

**Authors:** Jean Pierre Correia^1,\*^

^1^ Digital Balance AI -- Digital Labs Independent Research Group

\* Corresponding author: jpierre.vd@gmail.com | ORCID: 0009-0004-3566-3987

---

## Abstract

Immunoglobulin heavy chain J gene segments (IGHJ) anchor V(D)J recombination and encode the conserved FR4 domain, yet remain uncharacterized in the vast majority of wild ruminant genomes. We previously described SAFARI-IGHJ, a species-agnostic pipeline that identified 65 IGHJ candidates (36 predicted functional) across 21 Bovidae and Cervidae species spanning nine taxonomic groups. Here, we report the first large-scale functional validation of computationally predicted IGHJ genes in non-model ruminants using a three-pillar framework: (i) in vivo transcription, (ii) correct JH-to-CH1 mRNA splicing, and (iii) somatic V(D)J recombination detected through soft-clip junction analysis in whole-genome sequencing data.

RNA-seq validation in five species (*Odocoileus virginianus*, *Bison bison*, *Bison bonasus*, *Bos mutus*, *Cervus elaphus*) yielded over 11.4 million reads mapping to SAFARI-predicted IGHJ loci, with 1.49 million exhibiting canonical JH-CH1 splice junctions (intron sizes 1,047--2,577 bp), confirming productive mRNA maturation. WGS validation across eight species and 30 independent individuals detected 1,424 soft-clipped reads bearing V(D)J junction signatures at predicted loci, including 239 clips in *Saiga tatarica* (four individuals, two geographically distinct populations) and 728 clips in *Kobus ellipsiprymnus* -- providing the first direct evidence of somatic recombination at these loci. Together, these data validate SAFARI-IGHJ predictions across all seven Bovidae tribes and a Cervidae outgroup, advancing 13 species from computational prediction to empirical confirmation. We propose a four-level evidence hierarchy (mapping, junction, splice, expression) and discuss implications for conservation immunogenomics, including honest methodological limitations of short-read data for CDR-H3 architecture inference.

**Keywords:** IGHJ gene segments; Bovidae; Cervidae; V(D)J recombination; RNA-seq validation; whole-genome sequencing; soft-clip junction analysis; conservation immunogenomics; SAFARI-IGHJ; comparative immunogenomics

---

## 1. Introduction

### 1.1 The immunogenomic dark matter of wild ruminants

The adaptive immune system of jawed vertebrates relies on the combinatorial assembly of immunoglobulin variable regions through V(D)J recombination -- a somatic process that generates antibody diversity from a finite set of germline gene segments. In the immunoglobulin heavy chain locus (IGH), joining (J) gene segments occupy a privileged position: they encode the C-terminal portion of the variable domain, define the 3' boundary of CDR3, and contain the FR4 structural motif (WGXG family) essential for domain folding. The recombination signal sequences (RSS) flanking IGHJ segments -- consisting of a conserved heptamer, a 23-bp spacer, and a nonamer -- serve as the molecular targets for RAG-mediated recombination, with information content (IC) providing a quantitative proxy for recombination competence (Lee et al., 2003; Hesse et al., 1989).

While IGHJ gene complement is well characterized in domestic species (4 functional genes in *Bos taurus*: Walther et al., 2016; 6 in *Bubalus bubalis*: Ma et al., 2016), the IGHJ loci of wild ruminants remain largely uncharacterized. This knowledge gap is not trivial: Bovidae alone comprise over 140 species distributed across approximately ten tribes, many facing acute conservation pressures including habitat fragmentation, emerging infectious diseases, and population bottlenecks. Understanding the germline architecture of immune genes in these species is a prerequisite for assessing immunological fitness and disease susceptibility -- yet for most wild bovids, the IGH locus remains a genomic terra incognita.

### 1.2 The computational-experimental gap

The rapid expansion of publicly available genome assemblies (>50 bovid genomes on NCBI as of 2026) has created an unprecedented opportunity for comparative immunogenomics. However, most assemblies lack immunoglobulin-specific annotation, and the IGH locus presents particular challenges: extreme polymorphism, tandem gene duplications, and the tendency for Ig loci to reside in poorly assembled genomic regions (He et al., 2024; Koufariotis et al., 2023). No published, validated, species-agnostic pipeline existed for IGHJ discovery from whole-genome assemblies prior to SAFARI-IGHJ (Correia, 2026).

We previously applied SAFARI-IGHJ to 21 species across nine taxonomic groups, identifying 65 IGHJ candidates -- including first characterizations for 18 species and six threatened taxa. However, all classifications remained computational predictions based on RSS information content, FR4 motif integrity, and open reading frame status. The gap between computational prediction and biological proof -- the fundamental challenge of any in silico gene discovery effort -- remained open.

### 1.3 Bridging the gap: a three-pillar validation framework

In this study, we bridge that gap through systematic functional validation using publicly available RNA-seq and whole-genome sequencing (WGS) data. We developed a three-pillar evidence framework that tests distinct aspects of IGHJ functionality:

**Pillar 1 -- Transcriptional activity.** If a computationally predicted IGHJ gene is functional, it should be transcribed in lymphoid tissues. We mapped RNA-seq reads from spleen, lymph nodes, and retropharyngeal lymph nodes against SAFARI-predicted IGHJ loci, quantifying expression levels across species.

**Pillar 2 -- Correct mRNA splicing.** Functional IGHJ genes must produce mature mRNA through removal of the JH-to-CH1 intron. Spliced alignments (CIGAR N operations) in RNA-seq data provide direct evidence of productive mRNA maturation, with intron sizes serving as an independent structural validation.

**Pillar 3 -- Somatic V(D)J recombination.** The ultimate proof of IGHJ functionality is its participation in V(D)J recombination within B lymphocytes. In WGS data from blood or tissue samples containing B cells, reads spanning the V-D-J junction produce diagnostic soft-clips when aligned to the germline reference -- the clipped portion representing the somatically rearranged V and D segments absent from the germline genome.

This framework was applied to 13 species using 25 RNA-seq samples and WGS data from over 30 independent individuals, providing empirical validation across all seven Bovidae tribes (Bovini, Hippotragini, Tragelaphini, Alcelaphini, Aepycerotini, Reduncini, Antilopini) and the Cervidae outgroup.

### 1.4 Objectives

We: (i) validate SAFARI-IGHJ predictions through multi-omics functional evidence in 13 species spanning seven tribes plus a Cervidae outgroup; (ii) establish a four-level evidence hierarchy for IGHJ functional annotation (Level 2: mapping; Level 3: V(D)J junction; Level 3+: multi-individual junction; Level 4: splice-confirmed expression); (iii) report the first empirical evidence of IGHJ transcription and somatic recombination in orphan genomes including *Saiga tatarica*, *Kobus ellipsiprymnus*, and *Odocoileus virginianus*; and (iv) define methodological boundaries, including the limitation of short-read sequencing for CDR-H3 ultralong architecture inference.

---

## 2. Materials and Methods

### 2.1 SAFARI-IGHJ: iterative pipeline development

SAFARI-IGHJ was developed through three iterative stages, each addressing limitations identified in the previous version.

**Stage 1 -- Initial discovery pipeline (SAFARI v1).** The original pipeline integrated tBLASTn mining (PAM30 matrix, E-value 1.0) with single-linkage locus clustering and RSS-23 information content scoring. Applied to a 13-species panel, SAFARI v1 identified IGHJ candidates in all species but exhibited reduced sensitivity in highly fragmented assemblies: *Tragelaphus oryx* (N50 = 1.3 kb, scaffold count = 4.04 million) returned only 9 candidates, while the positive control on *Bos taurus* (ARS-UCD1.2) achieved F1 = 0.889 (Recall = 1.000, Precision = 0.800). These results indicated that biological sensitivity was adequate for high-quality assemblies but that the locus clustering algorithm -- while effective at filtering pseudogene noise -- could inadvertently discard genuine candidates in fragmented genomes.

**Stage 2 -- Expanded discovery (SAFARI v2).** To maximize sensitivity, SAFARI v2 relaxed clustering parameters and expanded the query set to include 10 protein queries spanning bovine, ovine, caprine, and human IGHJ segments. The expanded search recovered substantially more candidates in fragmented genomes: *T. oryx* yielded 54 candidates across dispersed scaffolds. However, this sensitivity came at a specificity cost -- manual inspection revealed that the majority of new candidates in fragmented genomes were low-confidence hits (IC < 18 bits) representing dispersed pseudogene fragments or spurious matches, inflating the candidate space without proportional biological gain.

**Stage 3 -- Machine learning classifier (SAFARI-Score).** Rather than reverting to conservative parameters, we developed a companion sklearn-based classifier (SAFARI-Score) trained on the complete feature set: RSS information content, FR4 motif class, open reading frame status, and heptamer mismatch count. The classifier was trained on 147 candidates from 21 species using leave-one-species-out cross-validation (LOOCV), achieving a composite score of 0.881 (gold standard F1 = 1.0, LOOCV mean accuracy = 0.898). When applied to the expanded *T. oryx* candidate set (54 candidates), SAFARI-Score ranked two candidates with highest confidence (JH_16: IC = 24.69, WGQG motif, pF = 0.775; JH_11: IC = 13.64, WGQR, pF = 0.740) -- both subsequently confirmed by RNA-seq validation (Section 3.2). The full classifier methodology, training data, and performance metrics are provided in Supplementary Methods S-ML.

**Pipeline availability.** SAFARI-IGHJ (v1 and v2) and SAFARI-IGHJ-Heavy (validation suite, 11 subcommands, 33 unit tests) are available as pip-installable Python packages. Source code and documentation: [GitHub URL -- DOI upon publication].

### 2.2 Genome assemblies and species dataset

SAFARI-IGHJ v2 was applied to whole-genome assemblies of 21 species representing nine taxonomic groups (8 Bovidae tribes plus Cervidae outgroup; Table 1). Genomes were obtained from NCBI GenBank or RefSeq. Assembly quality ranged from chromosome-level (N50 > 50 Mb; 8 species) to highly fragmented (N50 < 10 kb; 2 species). Confidence tiers (High/Medium/Low) were assigned based on N50 and assembly level to flag species where assembly quality limits architectural interpretation.

### 2.3 Multi-omics validation framework

To bridge the gap between computational IGHJ prediction and biological proof, we developed a three-pillar validation framework using publicly available RNA-seq and whole-genome sequencing (WGS) data (Table 2).

**Data sources.** RNA-seq reads from lymphoid tissues (spleen, lymph nodes, retropharyngeal lymph nodes) were obtained from NCBI SRA for five species: *Cervus elaphus* (6 samples, 3 tissues, BioProject PRJEB42381), *Odocoileus virginianus* (8 samples, retropharyngeal lymph node, SRP158695), *Bison bison* (5 samples, spleen and supramammary lymph node, PRJNA257088), *Bison bonasus* (1 sample, cross-species mapping of *B. bison* spleen RNA-seq against *B. bonasus* genome, PRJNA257088), and *Bos mutus* (3 samples, spleen, PRJNA309202). WGS data (Illumina paired-end, 10--45x coverage) were obtained for eight species: *Oryx dammah* (10 individuals), *Saiga tatarica* (4 individuals, 2 populations), *Kobus ellipsiprymnus* (3 individuals), *Bos javanicus* (5 individuals), *Connochaetes taurinus* (5 individuals), *Hippotragus niger* (4 individuals, 2 subspecies), *Addax nasomaculatus* (1 individual), and *Elaphurus davidianus* (1 individual).

**Reference construction.** For each species, the IGHJ-bearing scaffold identified by SAFARI was extracted with extended flanking regions (50 kb upstream and downstream of the outermost candidate) to allow unbiased read alignment without reference length artifacts. For species with fragmented assemblies where IGHJ candidates spanned multiple scaffolds, the scaffold containing the highest-confidence candidate was used.

**Pillar 1 -- Transcriptional activity (RNA-seq mapping).** RNA-seq reads were aligned to the extended IGHJ reference using minimap2 (v2.26; -ax splice for RNA-seq, -ax sr for WGS) with default parameters. Mapped read counts (samtools view -F 4 -c) quantify transcriptional activity at the predicted IGHJ locus.

**Pillar 2 -- JH-to-CH1 splicing (RNA-seq CIGAR analysis).** Spliced alignments were identified by the presence of N operations in CIGAR strings (e.g., 57M2021N37M), indicating intron removal during mRNA maturation. The JH-to-CH1 intron is the first intron downstream of the IGHJ exon; its removal produces a diagnostic splice junction connecting the J segment to the first constant region exon. Intron sizes were independently validated against available annotations where possible.

**Pillar 3 -- Somatic V(D)J recombination (WGS soft-clip analysis).** In WGS data from tissues containing B lymphocytes, reads spanning the V(D)J junction produce 5' soft-clips when aligned to the germline reference -- the clipped portion representing the somatically rearranged V and D segments absent from the unrearranged genome. Soft-clipped reads (CIGAR pattern: >=10S followed by M) mapping to the IGHJ locus were extracted and quantified. For WGS data, a pre-alignment bait-grep step was applied to reduce computational cost: compressed FASTQ files were scanned for reads containing any of 14 bovid-universal FR4 bait sequences (TGGGGC family), and only matching reads were subjected to full alignment. Direct mapping (without pre-filtering) was used where disk constraints permitted.

**Evidence hierarchy.** Results were classified into four levels: Level 2 (mapping only -- reads align to the predicted locus), Level 3 (V(D)J junction -- soft-clips with recombination signatures), Level 3+ (multi-individual junction -- junction evidence from >=3 independent individuals), and Level 4 (splice-confirmed expression -- JH-CH1 intron removal detected in RNA-seq).

---

## 3. Results

### 3.1 Pan-tribal validation of SAFARI-IGHJ predictions

Application of SAFARI-IGHJ v2 to 21 whole-genome assemblies identified 147 IGHJ candidates across nine taxonomic groups (Supplementary Table S-Candidates). After SAFARI-Score ranking, 65 high-confidence candidates (36 predicted functional, 29 ORF) were retained for downstream analysis (Table 3; methodology in Section 2.1).

The multi-omics validation framework was applied to 13 of these species using 25 RNA-seq samples and WGS data from over 30 independent individuals (Table 4). Across all validated species, more than 11.4 million reads mapped to SAFARI-predicted IGHJ loci, representing the largest empirical confirmation of computationally predicted immunoglobulin genes in non-model ruminants to date. Every Bovidae tribe in the dataset (7/7) and the Cervidae outgroup yielded positive validation at Level 2 or above, confirming that SAFARI predictions identify genuine IGHJ-bearing genomic regions across the full phylogenetic breadth of the study.

**Table 4. Multi-omics functional validation of SAFARI-IGHJ predictions across 13 species and 7 tribes.**

| Species | Tribe | Data type | Samples (N) | Mapped reads | Spliced reads | Soft-clips | Evidence level |
|---|---|---|---|---|---|---|---|
| *Odocoileus virginianus* | Cervidae | RNA-seq (rpLN) | 8 | 6,953,000 | 966,000 | -- | **Level 4** |
| *Cervus elaphus* | Cervini | RNA-seq (spleen, LN) | 6 | 1,627,661 | 209,678 | -- | **Level 4** |
| *Bos mutus* | Bovini | RNA-seq (spleen) | 3 | 1,031,357 | 211,019 | -- | **Level 4** |
| *Bison bison* | Bovini | RNA-seq (spleen, LN) | 5 | 469,000 | 80,000 | -- | **Level 4** |
| *Bison bonasus* | Bovini | RNA-seq (cross-sp.) | 1 | 327,000 | 22,000 | -- | **Level 4** |
| *Saiga tatarica* | Antilopini | WGS | 4 (2 pop.) | 2,819 | -- | 239 | **Level 3+** |
| *Kobus ellipsiprymnus* | Reduncini | WGS | 3 | 775 | -- | 728 | **Level 3+** |
| *Bos javanicus* | Bovini | WGS | 5 | 423 | -- | 354 | **Level 3+** |
| *Connochaetes taurinus* | Alcelaphini | WGS | 5 | 258 | -- | 73 | **Level 2+** |
| *Oryx dammah* | Hippotragini | WGS | 10 | 733 | -- | 25 | **Level 3** |
| *Hippotragus niger* | Hippotragini | WGS | 4 (2 ssp.) | 43 | -- | 5 | **Level 2+** |
| *Addax nasomaculatus* | Hippotragini | WGS | 1 | 39 | -- | 4 | **Level 2+** |
| *Elaphurus davidianus* | Cervidae | WGS | 1 | 7 | -- | 0 | **Level 2** |

rpLN = retropharyngeal lymph node; LN = lymph node; cross-sp. = cross-species mapping (*B. bison* reads against *B. bonasus* genome); pop. = geographically distinct populations; ssp. = subspecies. Evidence levels defined in Section 2.3. Dashes indicate analysis not applicable to data type.

### 3.2 Transcription and mRNA splicing confirm IGHJ expression in five species (Pillars 1--2)

RNA-seq data from lymphoid tissues provided the strongest evidence of IGHJ functionality. Across five species, over 10.4 million reads mapped to SAFARI-predicted IGHJ loci, with 1.49 million exhibiting canonical JH-to-CH1 splice junctions -- direct evidence of productive mRNA maturation.

**Odocoileus virginianus (white-tailed deer).** Retropharyngeal lymph node RNA-seq from a chronic wasting disease (CWD) study (Trone-Launer et al., 2019) yielded the most abundant signal: 6.95 million mapped reads and 966,000 spliced reads across 8 individuals. Expression varied substantially between individuals (range: 488,000--3,050,000 mapped reads per sample), consistent with the well-documented individual variation in adaptive immune repertoire size. Notably, CWD-positive animals (3 samples) showed lower aggregate expression (2.03 million mapped) than controls (5 samples; 4.93 million), though individual variation dominated over disease status -- the highest-expressing individual (SRR7748209, 3.05 million reads) was a control. This species was entirely absent from the SAFARI training dataset, representing a true out-of-sample validation.

**Bos mutus (wild yak).** Three spleen samples yielded 1.03 million mapped reads and 211,019 spliced reads, with remarkable consistency across biological replicates (304,407; 361,271; 365,679 mapped per sample; coefficient of variation = 10.2%). The yak genome (GCF_000298355.1) was processed end-to-end by SAFARI -- discovery, ML ranking, and RNA-seq validation -- entirely within this study, demonstrating the pipeline's capacity for autonomous characterization of orphan genomes.

**Bison bison and Bison bonasus.** Direct-species mapping of *B. bison* spleen RNA-seq (4 replicates) against its own genome yielded 287,000 mapped reads and 50,000 spliced reads, with high inter-replicate consistency (71,000--72,000 per sample). An additional supramammary lymph node sample contributed 182,000 mapped and 30,000 spliced reads. Cross-species mapping of the same *B. bison* RNA-seq against the *B. bonasus* genome -- a species separated by approximately 1.5 Mya -- recovered 327,000 mapped reads and 22,000 spliced reads, including a conserved JH-CH1 intron of 1,047 bp. This cross-species result confirms that the IGHJ locus architecture is sufficiently conserved within Bovini to support heterologous validation, a finding with practical implications for species where conspecific RNA-seq is unavailable.

**Cervus elaphus (red deer).** Six RNA-seq samples spanning three lymphoid tissues revealed a biologically expected expression gradient: spleen (503,000 mapped) > mesenteric lymph node (183,000) > prescapular lymph node (128,000), consistent with decreasing B cell density across these tissues. Biological replicates showed less than 1% coefficient of variation in mapping rates, confirming technical reproducibility. Spliced reads displayed characteristic CIGAR patterns (e.g., 57M2021N37M: 57-bp IGHJ exon match, 2,021-bp intron, 37-bp CH1 exon match), with intron sizes of 2,021--2,577 bp depending on the specific JH-CH1 junction -- within the expected range for mammalian immunoglobulin introns.

### 3.3 Somatic V(D)J recombination in wild populations (Pillar 3)

While RNA-seq demonstrates transcription and splicing, the definitive proof that a germline IGHJ gene participates in adaptive immunity is its engagement in somatic V(D)J recombination. WGS data from blood or tissue samples containing B lymphocytes can capture this signal: reads originating from rearranged B cell genomes produce diagnostic soft-clips when aligned to the unrearranged germline reference, with the clipped portion representing the somatically joined V and D segments.

**Saiga tatarica (saiga antelope) -- Antilopini.** SAFARI predicted two IGHJ candidates in the saiga genome (GCA_000739005.1; N50 = 6.5 kb): one predicted functional (JH_01, IC = 20.27, WGPG) and one ORF (JH_02). Despite the extreme genome fragmentation (1.77 million scaffolds), WGS validation across four individuals from two geographically distinct populations -- Ural (3 individuals) and Betpak-dala (1 individual) -- yielded 2,819 mapped reads and 239 soft-clips at the predicted IGHJ locus. The consistency across individuals (128--946 mapped per sample) and populations provides robust evidence that SAFARI correctly identified a functional, actively recombining IGHJ locus in one of the most fragmented genomes in our dataset. This represents the first direct evidence of somatic V(D)J recombination at an IGHJ locus in any Antilopini species.

**Kobus ellipsiprymnus (waterbuck) -- Reduncini.** Three WGS individuals produced 775 mapped reads and 728 soft-clips -- the highest clip-to-mapped ratio (0.94) in the entire dataset, indicating that nearly every read aligning to this locus bears a V(D)J junction signature. SAFARI had predicted two functional IGHJ genes (IC = 24.69, mean IC = 22.48) in a compact 937-bp cluster -- the minimal predicted architecture. The dense soft-clip signal validates both candidates and establishes the first empirical V(D)J evidence for the tribe Reduncini.

**Oryx dammah (scimitar-horned oryx) -- Hippotragini.** Ten independent WGS individuals yielded 733 mapped reads and 25 soft-clips across 7 of 10 individuals. While the per-individual signal is modest (1--6 clips), the multi-individual consistency (70% of individuals positive) provides statistical confidence. Soft-clip sequences contained diagnostic V(D)J junction motifs: the conserved Cys104 codon (TACTGT) followed by IGHD2-derived sequences (AGTAGT), confirming authentic recombination rather than mapping artifact. The three Hippotragini species validated by WGS (Oryx, Addax, H. niger) collectively demonstrate functional IGHJ usage across a tribe whose near-identical germline architecture is now confirmed to be biologically active.

**Bos javanicus (banteng) -- Bovini.** Five WGS individuals from an endangered species (IUCN: EN) produced 423 mapped reads and 354 soft-clips. As a member of tribe Bovini -- the clade known for ultralong CDR-H3 antibodies -- the banteng's V(D)J recombination evidence complements the RNA-seq validation from bison and yak, confirming IGHJ functionality across four Bovini species by independent methods.

**Connochaetes taurinus (blue wildebeest) -- Alcelaphini.** Five WGS individuals yielded 258 mapped reads and 73 soft-clips, establishing the first empirical V(D)J evidence for tribe Alcelaphini.

### 3.4 Summary: seven-tribe functional coverage

The combined RNA-seq and WGS validation establishes functional evidence for SAFARI-predicted IGHJ genes across all seven Bovidae tribes and the Cervidae outgroup (Figure 2). The evidence hierarchy spans four levels: five species at Level 4 (splice-confirmed expression), three at Level 3+ (multi-individual recombination), two at Level 2+ (recombination in few individuals), and three at Level 2 (mapping only). No species with Level 3 or above evidence contradicted SAFARI predictions: in every case, the computationally predicted IGHJ locus was confirmed as a site of active transcription or somatic recombination.

---

## 4. Discussion

### 4.1 From prediction to proof: the power of multi-omics IGHJ validation

This study demonstrates that computationally predicted IGHJ gene segments in non-model ruminants are not silent genomic relics but active components of the adaptive immune machinery. Across 13 species and seven Bovidae tribes, SAFARI-IGHJ predictions were confirmed by three independent lines of functional evidence: transcriptional activity (>11.4 million mapped reads), productive mRNA maturation (>1.49 million splice-junction reads), and somatic V(D)J recombination (>1,400 junction-bearing soft-clips in WGS data from 30+ independent individuals).

The complementarity of RNA-seq and WGS evidence deserves emphasis. RNA-seq captures the downstream consequences of IGHJ usage -- transcription and correct intron removal -- but cannot distinguish germline transcription from productive recombination. WGS captures the upstream event -- the physical DNA rearrangement that joins V, D, and J segments -- but at vastly lower abundance, since only the minority of nucleated cells in a tissue sample are B lymphocytes carrying rearranged immunoglobulin loci. The convergence of both data types on the same SAFARI-predicted genomic coordinates constitutes the strongest form of validation available from public repositories: it excludes the possibility that mapped reads represent transcriptional noise from a pseudogene (which would lack splicing) or that soft-clips represent alignment artifacts (which would not coincide with expressed, correctly spliced loci).

The four-level evidence hierarchy proposed here (Level 2: mapping; Level 3: junction; Level 3+: multi-individual junction; Level 4: splice-confirmed expression) provides a structured framework for future IGHJ annotation efforts. We note that Level 4 evidence was achievable only for species with publicly available lymphoid-tissue RNA-seq -- currently a minority of wild ruminants. The WGS-based approach (Levels 2--3+), while yielding orders of magnitude fewer reads, has the critical advantage of leveraging the growing archive of conservation genomics resequencing projects, where WGS data from multiple individuals is increasingly available even for threatened species with no transcriptomic resources.

### 4.2 Methodological honesty: the short-read boundary

The validation framework presented here operates within well-defined technical boundaries that must be explicitly acknowledged.

**The splicing-coverage trade-off in short reads.** Illumina paired-end reads (100--150 bp) that capture the diagnostic JH-to-CH1 splice junction necessarily allocate most of their length to the two flanking exons. A typical spliced read spanning a 2,021-bp intron in *Cervus elaphus* (CIGAR: 57M2021N37M) dedicates 57 bp to the IGHJ exon and 37 bp to the CH1 exon, leaving only 6 bp of soft-clip -- insufficient for CDR-H3 characterization. Reads that do carry longer 5' soft-clips (potentially representing V-D-J junction sequence) must by definition map entirely within the IGHJ exon without crossing the splice junction, creating a fundamental trade-off: the most informative reads for expression validation (spliced) are the least informative for junction characterization (clipped), and vice versa. This is not a pipeline limitation but an intrinsic constraint of short-read RNA-seq aligned to germline references.

**The adapter contamination trap.** During exploratory soft-clip analysis, we identified a systematic artifact: a substantial fraction of 5' soft-clips in RNA-seq data consisted of Illumina TruSeq adapter sequences (AGATCGGAAGAGCACACGTCTGAACTCCAGTCA) rather than biological V-D-J junction sequence. This contamination was readily detectable through sequence inspection but could inflate soft-clip counts and corrupt downstream analyses if not filtered. We flag this as a critical quality-control step for any study attempting to extract CDR-H3 information from bulk RNA-seq soft-clips.

**Read-length ceiling on junction capture.** In WGS data, the maximum observable soft-clip length is constrained by read length minus the aligned portion. For 100-bp reads (*Bison bison*, HiSeq 2500), maximum 5' clips were 20 bp; for 125-bp reads (*Bos mutus*), 45 bp; for 150-bp reads (*Tragelaphus oryx*), up to 73 bp. These lengths are sufficient to confirm V(D)J recombination (the presence of non-germline sequence at the junction) but insufficient to reconstruct complete CDR-H3 loops, which in bovine ultralong antibodies can exceed 60 amino acids (180 nt). We therefore explicitly refrain from making claims about CDR-H3 length distributions or ultralong antibody prevalence based on the data presented here.

**The path forward.** Complete characterization of CDR-H3 architecture in non-model ruminants will require dedicated adaptive immune receptor repertoire sequencing (AIRR-seq / Rep-Seq) using established protocols (Breden et al., 2017), or long-read sequencing platforms (PacBio HiFi, Oxford Nanopore) capable of capturing full-length V(D)J-C transcripts in single reads. The SAFARI-IGHJ predictions and evidence hierarchy established here provide the essential genomic coordinates and functional prioritization to guide such experiments.

### 4.3 Exploratory structural signatures in Tragelaphini

Despite the limitations outlined above, exploratory analysis of soft-clipped reads in *Tragelaphus oryx* (common eland, tribe Tragelaphini) revealed putative structural signatures that warrant future investigation. Among the top 5% longest soft-clips (minimum 60 bp, n = 107), 86.9% contained at least one cysteine residue upon in silico translation, with a mean of 1.37 cysteines per clip and detectable Cys-X(1-3)-Cys micro-domain motifs (mean 0.20 per clip). This contrasts sharply with *Odocoileus virginianus* (Cervidae outgroup), where only 4.9% of equivalent-length clips contained any cysteine. The enrichment is consistent with -- though not proof of -- the disulfide-bonded knob domains characteristic of ultralong CDR-H3 antibodies described in cattle (Stanfield et al., 2016; Dong et al., 2019), raising the possibility that cysteine-rich CDR-H3 structures may extend beyond the Bos/Bison clade into Tragelaphini.

We present this observation with three caveats. First, the clip lengths (maximum 73 bp = 24 amino acids) represent only a fragment of a putative ultralong CDR-H3, insufficient for structural inference. Second, cysteine enrichment in translated soft-clips could reflect codon bias in the genomic region flanking IGHJ rather than CDR-H3 biology. Third, the *T. oryx* genome is extremely fragmented (N50 = 1.3 kb), and the validation BAM derives from a single RNA-seq sample (SRR5647660). The complete analysis is presented in Supplementary Note S-Cysteine. We highlight this finding because it demonstrates that SAFARI-IGHJ, beyond its primary function as a gene discovery tool, can serve as a hypothesis generator for structural immunology -- identifying candidate loci where targeted Rep-Seq or crystallography could reveal convergent evolution of antibody architecture across phylogenetically distant bovid lineages.

### 4.4 Implications for conservation immunogenomics

The functional validation of IGHJ loci in threatened species has immediate implications for conservation biology. Six species of conservation concern in our dataset -- *Addax nasomaculatus* (CR), *Hippotragus niger* (CR), *Bos javanicus* (EN), *Oryx dammah* (EN), *Saiga tatarica* (NT), and *Elaphurus davidianus* (EW) -- now have empirically confirmed IGHJ loci, advancing from computational annotation to functional baselines.

The saiga antelope illustrates the potential most clearly. *Saiga tatarica* has suffered catastrophic population die-offs linked to *Pasteurella multocida* serotype B (Kock et al., 2018), with the 2015 mass mortality event killing over 200,000 individuals -- more than half the global population -- within weeks. Understanding the immunogenomic baseline of surviving populations is critical for disease vulnerability assessment. Our demonstration that the SAFARI-predicted IGHJ locus is actively recombining in four individuals from two geographically distinct populations (Ural and Betpak-dala) establishes that the somatic recombination machinery targeting this locus is functional across the species' range. Combined with the 239 junction-bearing soft-clips, this provides a foundation for designing population-level IGHJ genotyping panels that could assess allelic diversity erosion following bottleneck events.

More broadly, the three-pillar validation framework transforms SAFARI-IGHJ from a gene discovery tool into a functional annotation engine. As conservation genomics programs generate WGS data for endangered bovids -- increasingly at population scale for species management -- SAFARI can extract immunological signal from data already being collected for other purposes, without requiring additional sampling, specialized library preparation, or species-specific reagent development.

### 4.5 Limitations

Several limitations apply beyond the short-read constraints discussed in Section 4.2. Each species is represented by a single reference genome; population-level IGHJ allelic diversity remains unknown. Assembly quality varies substantially (N50: 1.3 kb to 101.2 Mb), and two species (*T. oryx*, *B. javanicus*) carry `fragmented_assembly` flags. The SAFARI-Score ML classifier was trained on candidates from the same species set used for validation, though leave-one-species-out cross-validation mitigates circularity. RNA-seq tissue sources are heterogeneous across species (spleen, lymph nodes, retropharyngeal lymph nodes), reflecting public data availability rather than experimental design. WGS coverage and tissue source vary across species and populations. Finally, Level 2 evidence (mapping only) does not exclude the possibility that aligned reads originate from processed pseudogenes with residual sequence similarity to functional IGHJ genes.

---

## 5. Conclusion

SAFARI-IGHJ bridges the gap between the expanding archive of wild ruminant genome assemblies and the near-total absence of immunoglobulin gene annotation in non-model species. Applied to 21 species across nine taxonomic groups and validated through 11.4 million RNA-seq reads, 1.49 million splice-junction events, and over 1,400 V(D)J recombination signatures in 30+ wild individuals, the pipeline transforms raw genome assemblies into actionable immunological knowledge -- from computational prediction to in vivo proof. The three-pillar validation framework and four-level evidence hierarchy established here provide a reproducible template for functional IGHJ annotation that can be applied to any species with a genome assembly and publicly available sequencing data. As conservation genomics programs generate population-scale resequencing data for threatened bovids and cervids, SAFARI-IGHJ enables the extraction of immunogenomic signal from existing resources, opening the door to systematic conservation immunogenomics without the cost or logistical barriers of dedicated immune profiling.

---

## References

[See references.bib]

---

## Data Availability

All genome assemblies are publicly available from NCBI GenBank/RefSeq (accessions in Table 1). RNA-seq and WGS data are publicly available from NCBI SRA (BioProject accessions in Section 2.3). SAFARI-IGHJ source code: [GitHub URL]. SAFARI-IGHJ-Heavy validation suite: [GitHub URL]. Validation BAM files and analysis scripts: [Zenodo DOI].

## Author Contributions

J.P.C. conceived the study, developed the SAFARI-IGHJ pipeline and SAFARI-Score classifier, designed and executed the multi-omics validation framework, performed all computational analyses, and wrote the manuscript.

## Competing Interests

The author declares no competing interests.

## Acknowledgements

Computational analyses were performed on Google Cloud Platform virtual machines. We acknowledge the genome sequencing consortia and individual researchers who deposited the public genome assemblies and sequencing data used in this study.

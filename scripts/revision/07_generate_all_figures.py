#!/usr/bin/env python3
"""
07_generate_all_figures.py — Publication-quality figures for revised manuscript.
Addresses Reviewer 2 criticism #8: "Missing figures and tables"
"""
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import seaborn as sns

BASE = Path(__file__).parent.parent
DATA = BASE / "data"
OUT_F = BASE / "outputs" / "figures"
OUT_F.mkdir(parents=True, exist_ok=True)

CANDIDATES = Path("/Users/jpierrevd/Documents/projectos  AI Science/safari_autoresearch_v2/phase2/candidates_all_with_embeddings.tsv")

# Color palette (colorblind-safe)
COLORS = {
    "Level 4": "#1B5E20",   # Dark green
    "Level 3+": "#4CAF50",  # Green
    "Level 3": "#8BC34A",   # Light green
    "Level 2+": "#FFC107",  # Amber
    "Level 2": "#FF9800",   # Orange
    "RNA-seq": "#1565C0",   # Blue
    "WGS": "#E65100",       # Deep orange
}

TRIBE_ORDER = ["Bovini", "Hippotragini", "Tragelaphini", "Alcelaphini",
               "Aepycerotini", "Reduncini", "Antilopini", "Cervidae"]


def fig2_evidence_by_tribe():
    """Fig 2: Evidence levels by tribe (stacked horizontal bar)."""
    val = pd.read_csv(DATA / "validation_summary.tsv", sep="\t")

    # Map evidence levels to labels
    def level_label(lev):
        if lev >= 4: return "Level 4"
        if lev >= 3.5: return "Level 3+"
        if lev >= 3: return "Level 3"
        if lev >= 2.5: return "Level 2+"
        return "Level 2"

    val["level_label"] = val["evidence_level"].apply(level_label)

    fig, ax = plt.subplots(figsize=(10, 5))
    sns.set_style("whitegrid")

    tribes = []
    for tribe in TRIBE_ORDER:
        sub = val[val.tribe == tribe]
        if len(sub) > 0:
            tribes.append(tribe)

    level_order = ["Level 4", "Level 3+", "Level 3", "Level 2+", "Level 2"]
    bottom = np.zeros(len(tribes))

    for level in level_order:
        counts = []
        for tribe in tribes:
            sub = val[(val.tribe == tribe) & (val.level_label == level)]
            counts.append(len(sub))
        counts = np.array(counts)
        if counts.sum() > 0:
            ax.barh(range(len(tribes)), counts, left=bottom,
                    color=COLORS[level], label=level, edgecolor="white", linewidth=0.5)
            # Add count labels
            for i, c in enumerate(counts):
                if c > 0:
                    ax.text(bottom[i] + c/2, i, str(c), ha="center", va="center",
                            fontsize=9, fontweight="bold", color="white")
            bottom += counts

    ax.set_yticks(range(len(tribes)))
    ax.set_yticklabels(tribes, fontsize=11)
    ax.set_xlabel("Number of Validated Species", fontsize=12)
    ax.set_title("Functional Evidence Levels Across Bovidae Tribes and Cervidae Outgroup",
                 fontsize=13, fontweight="bold")
    ax.legend(loc="lower right", fontsize=9, framealpha=0.9)
    ax.set_xlim(0, max(bottom) + 0.5)

    plt.tight_layout()
    fig.savefig(OUT_F / "fig2_evidence_by_tribe.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_F / "fig2_evidence_by_tribe.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  Fig 2: evidence_by_tribe OK")


def fig3_rnaseq_mapping():
    """Fig 3: Mapped vs Spliced reads for RNA-seq species."""
    val = pd.read_csv(DATA / "validation_summary.tsv", sep="\t")
    rna = val[val.data_type == "RNA-seq"].sort_values("mapped", ascending=True)

    fig, ax = plt.subplots(figsize=(10, 5))
    sns.set_style("whitegrid")

    species_labels = [s.replace("_", " ").replace("bison bonasus", "B. bonasus\n(cross-species)")
                      for s in rna["species"]]
    species_labels = [f"*{s.split()[0][0]}. {' '.join(s.split()[1:])}*" if "cross" not in s
                      else s for s in species_labels]
    # Simplified labels
    labels = []
    for _, r in rna.iterrows():
        sp = r["species"].replace("_", " ")
        parts = sp.split()
        label = f"{parts[0][0]}. {parts[1]}"
        if r["species"] == "bison_bonasus":
            label += "\n(cross-sp.)"
        labels.append(label)

    y = np.arange(len(rna))
    height = 0.35

    bars1 = ax.barh(y - height/2, rna["mapped"] / 1e6, height,
                     color="#1565C0", label="Mapped reads", edgecolor="white")
    bars2 = ax.barh(y + height/2, rna["spliced"] / 1e6, height,
                     color="#E91E63", label="Spliced reads (JH-CH1)", edgecolor="white")

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=10, fontstyle="italic")
    ax.set_xlabel("Reads (millions)", fontsize=12)
    ax.set_title("RNA-seq Validation: Mapped and Spliced Reads at SAFARI-Predicted IGHJ Loci",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=10, loc="lower right")

    # Add value labels
    for bar in bars1:
        w = bar.get_width()
        if w > 0.1:
            ax.text(w + 0.05, bar.get_y() + bar.get_height()/2,
                    f"{w:.1f}M", va="center", fontsize=8, color="#1565C0")

    for bar in bars2:
        w = bar.get_width()
        if w > 0.01:
            ax.text(w + 0.05, bar.get_y() + bar.get_height()/2,
                    f"{w*1000:.0f}K", va="center", fontsize=8, color="#E91E63")

    plt.tight_layout()
    fig.savefig(OUT_F / "fig3_rnaseq_mapping.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_F / "fig3_rnaseq_mapping.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  Fig 3: rnaseq_mapping OK")


def fig4_wgs_softclips():
    """Fig 4: WGS soft-clip counts across species and individuals."""
    val = pd.read_csv(DATA / "validation_summary.tsv", sep="\t")
    wgs = val[(val.data_type == "WGS") & (val.softclips > 0)].sort_values("softclips", ascending=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [2, 1]})
    sns.set_style("whitegrid")

    # Panel A: Soft-clip counts
    labels = []
    for _, r in wgs.iterrows():
        sp = r["species"].replace("_", " ").split()
        label = f"{sp[0][0]}. {sp[1]}\n(n={int(r['samples_n'])})"
        labels.append(label)

    colors = [COLORS.get({"3.5": "Level 3+", "3.0": "Level 3", "3": "Level 3",
                          "2.5": "Level 2+"}.get(str(r["evidence_level"]), "Level 2"), "#999")
              for _, r in wgs.iterrows()]

    ax1.barh(range(len(wgs)), wgs["softclips"], color=colors, edgecolor="white", linewidth=0.5)
    ax1.set_yticks(range(len(wgs)))
    ax1.set_yticklabels(labels, fontsize=9, fontstyle="italic")
    ax1.set_xlabel("V(D)J Junction Soft-Clips", fontsize=11)
    ax1.set_title("A. Somatic V(D)J Recombination Evidence\n(WGS Soft-Clips at IGHJ Loci)", fontsize=11, fontweight="bold")

    # Add values
    for i, (_, r) in enumerate(wgs.iterrows()):
        ax1.text(r["softclips"] + 5, i, str(int(r["softclips"])), va="center", fontsize=9)

    # Panel B: Clip-to-mapped ratio
    wgs_plot = wgs.copy()
    wgs_plot["clip_ratio"] = wgs_plot["softclips"] / wgs_plot["mapped"]
    ax2.barh(range(len(wgs_plot)), wgs_plot["clip_ratio"], color=colors, edgecolor="white", linewidth=0.5)
    ax2.set_yticks(range(len(wgs_plot)))
    ax2.set_yticklabels(["" for _ in range(len(wgs_plot))])
    ax2.set_xlabel("Soft-Clip / Mapped Ratio", fontsize=11)
    ax2.set_title("B. Junction Fraction\n(Higher = More Active Recombination)", fontsize=11, fontweight="bold")

    for i, (_, r) in enumerate(wgs_plot.iterrows()):
        ax2.text(r["clip_ratio"] + 0.02, i, f"{r['clip_ratio']:.2f}", va="center", fontsize=9)

    plt.tight_layout()
    fig.savefig(OUT_F / "fig4_wgs_softclips.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_F / "fig4_wgs_softclips.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  Fig 4: wgs_softclips OK")


def fig1_pipeline_schematic():
    """Fig 1: Pipeline schematic (SAFARI v1 -> v2 -> Score -> 3 pillars)."""
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 7)
    ax.axis("off")

    # Title
    ax.text(7, 6.7, "SAFARI-IGHJ: Pipeline Architecture and Three-Pillar Validation Framework",
            ha="center", fontsize=14, fontweight="bold")

    # --- ACT 1: Pipeline Evolution ---
    ax.text(3.5, 6.2, "ACT 1: Computational Discovery", ha="center", fontsize=12,
            fontweight="bold", color="#1565C0")

    boxes = [
        (0.5, 5.0, 2.0, 0.8, "#E3F2FD", "SAFARI v1\ntBLASTn + RSS\n(13 species)", "#1565C0"),
        (3.5, 5.0, 2.0, 0.8, "#E8F5E9", "SAFARI v2\nExpanded queries\n(21 species)", "#2E7D32"),
        (6.5, 5.0, 2.0, 0.8, "#FFF3E0", "SAFARI-Score\nsklearn ML\n(LOOCV=0.898)", "#E65100"),
    ]
    for x, y, w, h, fc, text, ec in boxes:
        box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1", facecolor=fc, edgecolor=ec, linewidth=2)
        ax.add_patch(box)
        ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=8, fontweight="bold")

    # Arrows between boxes
    for x1, x2 in [(2.5, 3.5), (5.5, 6.5)]:
        ax.annotate("", xy=(x2, 5.4), xytext=(x1, 5.4),
                    arrowprops=dict(arrowstyle="->", color="#333", lw=2))

    # Output box
    box = FancyBboxPatch((9.5, 5.0), 2.5, 0.8, boxstyle="round,pad=0.1",
                          facecolor="#F3E5F5", edgecolor="#7B1FA2", linewidth=2)
    ax.add_patch(box)
    ax.text(10.75, 5.4, "65 IGHJ Candidates\n36 Functional\n21 Species, 9 Groups",
            ha="center", va="center", fontsize=8, fontweight="bold", color="#7B1FA2")
    ax.annotate("", xy=(9.5, 5.4), xytext=(8.5, 5.4),
                arrowprops=dict(arrowstyle="->", color="#333", lw=2))

    # --- ACT 2: Validation ---
    ax.text(7, 3.8, "ACT 2: Multi-Omics Functional Validation", ha="center", fontsize=12,
            fontweight="bold", color="#C62828")

    # Three pillars
    pillars = [
        (1.0, 1.5, 3.5, 1.8, "#E3F2FD", "#1565C0",
         "PILLAR 1: Transcription\n(RNA-seq Mapping)\n\n11.4M reads mapped\n5 species, 25 samples\nLevel 4 evidence"),
        (5.25, 1.5, 3.5, 1.8, "#FCE4EC", "#C62828",
         "PILLAR 2: Splicing\n(JH-CH1 Intron Removal)\n\n1.49M splice junctions\nIntrons 1,047-2,577 bp\nCIGAR N operations"),
        (9.5, 1.5, 3.5, 1.8, "#FFF3E0", "#E65100",
         "PILLAR 3: Recombination\n(WGS Soft-Clips)\n\n1,424 V(D)J junctions\n30+ individuals, 8 species\nCIGAR S operations"),
    ]
    for x, y, w, h, fc, ec, text in pillars:
        box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.1", facecolor=fc, edgecolor=ec, linewidth=2)
        ax.add_patch(box)
        ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=7.5, fontweight="bold")

    # Arrow from output to pillars
    for px in [2.75, 7.0, 11.25]:
        ax.annotate("", xy=(px, 3.3), xytext=(10.75, 5.0),
                    arrowprops=dict(arrowstyle="->", color="#333", lw=1.5, connectionstyle="arc3,rad=0.1"))

    # Bottom: result
    box = FancyBboxPatch((3.0, 0.2), 8.0, 0.9, boxstyle="round,pad=0.1",
                          facecolor="#E8F5E9", edgecolor="#1B5E20", linewidth=3)
    ax.add_patch(box)
    ax.text(7, 0.65, "RESULT: 13 Species Validated | 7/7 Bovidae Tribes + Cervidae Outgroup\n"
            "Level 4 (5 spp.) | Level 3+ (3 spp.) | Level 3 (1 sp.) | Level 2+ (3 spp.) | Level 2 (1 sp.)",
            ha="center", va="center", fontsize=9, fontweight="bold", color="#1B5E20")

    for px in [2.75, 7.0, 11.25]:
        ax.annotate("", xy=(7, 1.1), xytext=(px, 1.5),
                    arrowprops=dict(arrowstyle="->", color="#1B5E20", lw=1.5))

    plt.tight_layout()
    fig.savefig(OUT_F / "fig1_pipeline_schematic.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_F / "fig1_pipeline_schematic.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  Fig 1: pipeline_schematic OK")


def fig5_phylogeny_evidence():
    """Fig 5: Simplified phylogenetic layout with evidence levels."""
    val = pd.read_csv(DATA / "validation_summary.tsv", sep="\t")

    fig, ax = plt.subplots(figsize=(12, 8))
    sns.set_style("white")
    ax.set_xlim(0, 12)
    ax.set_ylim(-0.5, len(val) + 0.5)
    ax.axis("off")

    # Sort by tribe then evidence level
    tribe_order_map = {t: i for i, t in enumerate(TRIBE_ORDER)}
    val["tribe_rank"] = val["tribe"].map(tribe_order_map)
    val = val.sort_values(["tribe_rank", "evidence_level"], ascending=[True, False]).reset_index(drop=True)

    level_colors = {4: "#1B5E20", 3.5: "#4CAF50", 3: "#8BC34A", 2.5: "#FFC107", 2: "#FF9800"}

    current_tribe = None
    for i, (_, row) in enumerate(val.iterrows()):
        y = len(val) - 1 - i

        # Tribe separator
        if row["tribe"] != current_tribe:
            current_tribe = row["tribe"]
            ax.text(0.5, y + 0.35, current_tribe, fontsize=10, fontweight="bold", color="#333")

        # Species name (italic)
        sp = row["species"].replace("_", " ")
        parts = sp.split()
        label = f"{parts[0][0]}. {parts[1]}"
        ax.text(1.0, y, label, fontsize=9, fontstyle="italic", va="center")

        # Evidence level colored circle
        color = level_colors.get(row["evidence_level"], "#999")
        ax.scatter(4.5, y, s=200, c=color, edgecolors="white", linewidth=1, zorder=5)

        # Level label
        lev = row["evidence_level"]
        lev_str = {4: "L4", 3.5: "L3+", 3: "L3", 2.5: "L2+", 2: "L2"}.get(lev, "?")
        ax.text(4.5, y, lev_str, ha="center", va="center", fontsize=7, fontweight="bold", color="white")

        # Data type
        ax.text(5.5, y, row["data_type"], fontsize=8, va="center", color="#666")

        # Key metric
        if row["data_type"] == "RNA-seq":
            metric = f"{row['mapped']/1e6:.1f}M mapped, {row['spliced']/1e3:.0f}K spliced"
        else:
            metric = f"{int(row['mapped'])} mapped, {int(row['softclips'])} clips (n={int(row['samples_n'])})"
        ax.text(6.5, y, metric, fontsize=8, va="center")

    # Legend
    for i, (lev, color) in enumerate(level_colors.items()):
        lev_str = {4: "Level 4: Splice-confirmed", 3.5: "Level 3+: Multi-individual junction",
                   3: "Level 3: V(D)J junction", 2.5: "Level 2+: Junction (few ind.)",
                   2: "Level 2: Mapping only"}.get(lev, "?")
        ax.scatter(9.5, len(val) - 1 - i * 0.8, s=120, c=color, edgecolors="white")
        ax.text(10.0, len(val) - 1 - i * 0.8, lev_str, fontsize=8, va="center")

    ax.set_title("Phylogenetic Coverage: SAFARI-IGHJ Validation Across 7 Bovidae Tribes + Cervidae",
                 fontsize=13, fontweight="bold", pad=20)

    plt.tight_layout()
    fig.savefig(OUT_F / "fig5_phylogeny_evidence.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_F / "fig5_phylogeny_evidence.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("  Fig 5: phylogeny_evidence OK")


def main():
    print("=" * 60)
    print("GENERATING PUBLICATION FIGURES")
    print("=" * 60)
    fig1_pipeline_schematic()
    fig2_evidence_by_tribe()
    fig3_rnaseq_mapping()
    fig4_wgs_softclips()
    fig5_phylogeny_evidence()
    print(f"\nAll figures saved to: {OUT_F}")
    print(f"Files: {list(OUT_F.glob('*.pdf'))}")

if __name__ == "__main__":
    main()

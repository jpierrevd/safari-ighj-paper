#!/usr/bin/env python3
"""
06_biological_insights.py — Evolutionary and biological analyses.
Addresses Reviewer 2 criticism #13: "Limited biological insights"
"""
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

BASE = Path(__file__).parent.parent
OUT_T = BASE / "outputs" / "tables"
OUT_F = BASE / "outputs" / "figures"

CANDIDATES = Path("/Users/jpierrevd/Documents/projectos  AI Science/safari_autoresearch_v2/phase2/candidates_all_with_embeddings.tsv")
PREDICTIONS = Path("/Users/jpierrevd/Documents/projectos  AI Science/safari_autoresearch_v2/phase2/sklearn_predictions_all147.tsv")

# Known Bos taurus IGHJ (IMGT reference)
BOS_TAURUS_IGHJ = {
    "IGHJ1": {"fr4": "WGQG", "ic": 24.69, "status": "Functional"},
    "IGHJ2": {"fr4": "WGQG", "ic": 20.27, "status": "Functional"},
    "IGHJ3": {"fr4": "WGPG", "ic": 20.27, "status": "Functional"},
    "IGHJ4": {"fr4": "WGRG", "ic": 15.85, "status": "Functional"},
}

def main():
    cand = pd.read_csv(CANDIDATES, sep="\t")
    pred = pd.read_csv(PREDICTIONS, sep="\t")
    merged = cand.merge(pred[["candidate_id", "classification_pred", "prob_functional"]], on="candidate_id", how="left")

    results = []

    # 1. IGHJ repertoire size across tribes
    print("=" * 60)
    print("BIOLOGICAL INSIGHTS — REVIEWER 2 RESPONSE")
    print("=" * 60)

    # Count functional candidates per species
    func = merged[merged.classification_safari == "Functional"].copy()
    species_counts = func.groupby(["species", "tribe"]).size().reset_index(name="n_functional")
    tribe_counts = species_counts.groupby("tribe")["n_functional"].apply(list).to_dict()

    print("\n1. IGHJ FUNCTIONAL GENE COUNT BY TRIBE")
    print("-" * 40)
    tribe_data = []
    for tribe, counts in sorted(tribe_counts.items()):
        mean_c = np.mean(counts)
        std_c = np.std(counts) if len(counts) > 1 else 0
        print(f"  {tribe:20s}: {counts} (mean={mean_c:.1f} ± {std_c:.1f}, n_species={len(counts)})")
        tribe_data.append({"tribe": tribe, "n_species": len(counts), "mean_functional": mean_c,
                          "sd_functional": std_c, "counts": str(counts)})

    # Kruskal-Wallis test (non-parametric ANOVA)
    groups = [np.array(counts) for counts in tribe_counts.values() if len(counts) >= 1]
    if len(groups) >= 3:
        # Need at least 2 observations per group for meaningful test
        groups_valid = [g for g in groups if len(g) >= 2]
        if len(groups_valid) >= 3:
            h_stat, h_p = stats.kruskal(*groups_valid)
            print(f"\n  Kruskal-Wallis H={h_stat:.2f}, p={h_p:.4f}")
        else:
            h_stat, h_p = np.nan, np.nan
            print(f"\n  Insufficient data for Kruskal-Wallis (need ≥2 species per tribe)")
    results.append({"analysis": "IGHJ repertoire ANOVA", "statistic": f"H={h_stat:.2f}" if not np.isnan(h_stat) else "NA",
                    "p_value": h_p if not np.isnan(h_p) else "NA", "note": "Kruskal-Wallis across tribes"})

    # 2. FR4 motif distribution
    print("\n2. FR4 MOTIF DISTRIBUTION")
    print("-" * 40)
    fr4_by_tribe = merged.groupby(["tribe", "fr4_motif"]).size().reset_index(name="count")
    fr4_pivot = fr4_by_tribe.pivot_table(index="tribe", columns="fr4_motif", values="count", fill_value=0)
    print(fr4_pivot.to_string())

    # 3. RSS IC distribution: Functional vs Pseudogene
    print("\n3. RSS INFORMATION CONTENT")
    print("-" * 40)
    for cls in ["Functional", "ORF", "Pseudogene"]:
        sub = merged[merged.classification_safari == cls]["rss_ic"]
        if len(sub) > 0:
            print(f"  {cls:15s}: median={sub.median():.1f}, mean={sub.mean():.1f}, "
                  f"range=[{sub.min():.1f}-{sub.max():.1f}], n={len(sub)}")

    # Mann-Whitney: Functional IC > Pseudogene IC
    func_ic = merged[merged.classification_safari == "Functional"]["rss_ic"]
    pseudo_ic = merged[merged.classification_safari == "Pseudogene"]["rss_ic"]
    if len(func_ic) >= 3 and len(pseudo_ic) >= 3:
        u_stat, u_p = stats.mannwhitneyu(func_ic, pseudo_ic, alternative="greater")
        print(f"\n  Mann-Whitney U (Functional > Pseudogene): U={u_stat:.0f}, p={u_p:.2e}")
        results.append({"analysis": "IC Functional > Pseudogene", "statistic": f"U={u_stat:.0f}",
                        "p_value": f"{u_p:.2e}", "note": "Mann-Whitney U, one-sided"})

    # 4. Bos taurus comparison: which wild species share which IGHJ orthologs?
    print("\n4. BOS TAURUS IGHJ ORTHOLOG COMPARISON")
    print("-" * 40)
    print("  Known Bos taurus functional IGHJ: 4 genes (IGHJ1-4)")
    print("  FR4 motifs: WGQG (×2), WGPG (×1), WGRG (×1)")

    # Count FR4 motif matches to Bos taurus
    bt_motifs = Counter([v["fr4"] for v in BOS_TAURUS_IGHJ.values()])
    for tribe in sorted(merged["tribe"].unique()):
        tribe_func = merged[(merged.tribe == tribe) & (merged.classification_safari == "Functional")]
        tribe_motifs = Counter(tribe_func["fr4_motif"].values)
        shared = sum(min(bt_motifs[m], tribe_motifs.get(m, 0)) for m in bt_motifs)
        print(f"  {tribe:20s}: {len(tribe_func)} functional, FR4={dict(tribe_motifs)}, "
              f"shared_with_Bt={shared}/4")

    # 5. Heptamer conservation across candidates
    print("\n5. HEPTAMER MISMATCH DISTRIBUTION")
    print("-" * 40)
    for cls in ["Functional", "ORF", "Pseudogene"]:
        sub = merged[merged.classification_safari == cls]["heptamer_mm"]
        if len(sub) > 0:
            print(f"  {cls:15s}: median={sub.median():.0f}, mean={sub.mean():.1f}, "
                  f"range=[{sub.min():.0f}-{sub.max():.0f}]")

    # Generate figure: IC distribution by classification
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    sns.set_style("whitegrid")

    # Panel A: IC by classification
    ax = axes[0]
    class_order = ["Functional", "ORF", "Pseudogene"]
    class_colors = {"Functional": "#4CAF50", "ORF": "#FFC107", "Pseudogene": "#F44336"}
    for cls in class_order:
        sub = merged[merged.classification_safari == cls]["rss_ic"]
        ax.hist(sub, bins=15, alpha=0.6, label=f"{cls} (n={len(sub)})", color=class_colors[cls])
    ax.set_xlabel("RSS Information Content (bits)", fontsize=11)
    ax.set_ylabel("Count", fontsize=11)
    ax.set_title("A. RSS IC Distribution by Classification", fontsize=11)
    ax.legend(fontsize=9)

    # Panel B: Functional gene count by tribe
    ax = axes[1]
    tribe_df = pd.DataFrame(tribe_data).sort_values("mean_functional", ascending=True)
    ax.barh(range(len(tribe_df)), tribe_df["mean_functional"],
            xerr=tribe_df["sd_functional"], color="#2196F3", alpha=0.8, capsize=3)
    ax.set_yticks(range(len(tribe_df)))
    ax.set_yticklabels(tribe_df["tribe"], fontsize=9)
    ax.set_xlabel("Mean Functional IGHJ Genes per Species", fontsize=11)
    ax.set_title("B. IGHJ Repertoire Size by Tribe", fontsize=11)

    # Panel C: FR4 motif distribution heatmap
    ax = axes[2]
    func_only = merged[merged.classification_safari == "Functional"]
    motif_tribe = func_only.groupby(["tribe", "fr4_motif"]).size().unstack(fill_value=0)
    sns.heatmap(motif_tribe, annot=True, fmt="d", cmap="YlOrRd", ax=ax, cbar_kws={"label": "Count"})
    ax.set_title("C. FR4 Motif Distribution (Functional Only)", fontsize=11)
    ax.set_ylabel("")

    plt.suptitle("Biological Insights: IGHJ Gene Architecture Across Bovidae Tribes",
                 fontsize=13, fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(OUT_F / "fig5_biological_insights.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_F / "fig5_biological_insights.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Save results table
    results_df = pd.DataFrame(results)
    results_df.to_csv(OUT_T / "S-Biology.tsv", sep="\t", index=False)
    tribe_df_out = pd.DataFrame(tribe_data)
    tribe_df_out.to_csv(OUT_T / "S-RepertoireByTribe.tsv", sep="\t", index=False)

    print(f"\nFigure: {OUT_F / 'fig5_biological_insights.pdf'}")
    print(f"Tables: {OUT_T / 'S-Biology.tsv'}, {OUT_T / 'S-RepertoireByTribe.tsv'}")

if __name__ == "__main__":
    main()

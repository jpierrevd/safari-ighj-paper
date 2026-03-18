#!/usr/bin/env python3
"""
05_data_quality_correlations.py — Data quality vs validation success.
Addresses Reviewer 2 criticism #3: "Heterogeneous and uncontrolled data"
"""
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

BASE = Path(__file__).parent.parent
DATA = BASE / "data"
OUT_T = BASE / "outputs" / "tables"
OUT_F = BASE / "outputs" / "figures"

def main():
    val = pd.read_csv(DATA / "validation_summary.tsv", sep="\t")
    val["n50_kb"] = pd.to_numeric(val["n50_kb"], errors="coerce")
    val["log_mapped"] = np.log10(val["mapped"].clip(lower=1))
    val["log_n50"] = np.log10(val["n50_kb"].clip(lower=0.1))

    sns.set_style("whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Panel A: Evidence level vs log(N50)
    ax = axes[0, 0]
    colors = {"RNA-seq": "#2196F3", "WGS": "#FF9800"}
    for dt in val["data_type"].unique():
        sub = val[val.data_type == dt]
        ax.scatter(sub["log_n50"], sub["evidence_level"], c=colors[dt], s=80, alpha=0.8, label=dt, edgecolors="white")
    rho, p = stats.spearmanr(val["log_n50"].dropna(), val.loc[val["log_n50"].notna(), "evidence_level"])
    ax.set_xlabel("log₁₀(N50, kb)", fontsize=11)
    ax.set_ylabel("Evidence Level", fontsize=11)
    ax.set_title(f"A. Evidence Level vs Assembly Quality\n(Spearman ρ={rho:.3f}, p={p:.3f})", fontsize=11)
    ax.legend(fontsize=9)
    ax.set_yticks([2, 2.5, 3, 3.5, 4])

    # Panel B: Evidence level vs sample count
    ax = axes[0, 1]
    for dt in val["data_type"].unique():
        sub = val[val.data_type == dt]
        ax.scatter(sub["samples_n"], sub["evidence_level"], c=colors[dt], s=80, alpha=0.8, label=dt, edgecolors="white")
    rho2, p2 = stats.spearmanr(val["samples_n"], val["evidence_level"])
    ax.set_xlabel("Number of Samples", fontsize=11)
    ax.set_ylabel("Evidence Level", fontsize=11)
    ax.set_title(f"B. Evidence Level vs Sample Count\n(Spearman ρ={rho2:.3f}, p={p2:.3f})", fontsize=11)
    ax.legend(fontsize=9)

    # Panel C: log(mapped) vs evidence level
    ax = axes[1, 0]
    for dt in val["data_type"].unique():
        sub = val[val.data_type == dt]
        ax.scatter(sub["evidence_level"], sub["log_mapped"], c=colors[dt], s=80, alpha=0.8, label=dt, edgecolors="white")
    rho3, p3 = stats.spearmanr(val["evidence_level"], val["log_mapped"])
    ax.set_xlabel("Evidence Level", fontsize=11)
    ax.set_ylabel("log₁₀(Mapped Reads)", fontsize=11)
    ax.set_title(f"C. Mapped Reads vs Evidence Level\n(Spearman ρ={rho3:.3f}, p={p3:.3f})", fontsize=11)
    ax.legend(fontsize=9)

    # Panel D: Tissue type effect (RNA-seq only)
    ax = axes[1, 1]
    rna = val[val.data_type == "RNA-seq"].copy()
    tissue_order = rna.sort_values("mapped", ascending=False)["tissue"].values
    bar_colors = ["#1565C0", "#42A5F5", "#90CAF9", "#BBDEFB", "#E3F2FD"]
    ax.barh(range(len(rna)), rna.sort_values("mapped", ascending=True)["mapped"] / 1e6,
            color=bar_colors[:len(rna)], edgecolor="white")
    ax.set_yticks(range(len(rna)))
    ax.set_yticklabels([f"{s.replace('_',' ').title()}\n({t})" for s, t in
                        zip(rna.sort_values("mapped", ascending=True)["species"],
                            rna.sort_values("mapped", ascending=True)["tissue"])], fontsize=9)
    ax.set_xlabel("Mapped Reads (millions)", fontsize=11)
    ax.set_title("D. RNA-seq Mapping by Species and Tissue", fontsize=11)

    plt.suptitle("Data Quality and Validation Success Correlations", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(OUT_F / "S-QualityCorrelations.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_F / "S-QualityCorrelations.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Save correlation table
    corr_df = pd.DataFrame([
        {"comparison": "evidence_level vs log10(N50)", "rho": rho, "p_value": p, "n": len(val.dropna(subset=["n50_kb"]))},
        {"comparison": "evidence_level vs sample_count", "rho": rho2, "p_value": p2, "n": len(val)},
        {"comparison": "evidence_level vs log10(mapped)", "rho": rho3, "p_value": p3, "n": len(val)},
    ])
    corr_df.to_csv(OUT_T / "S-QualityCorr.tsv", sep="\t", index=False)

    print("=" * 60)
    print("DATA QUALITY CORRELATIONS")
    print("=" * 60)
    for _, r in corr_df.iterrows():
        sig = "***" if r["p_value"] < 0.001 else ("*" if r["p_value"] < 0.05 else "ns")
        print(f"  [{sig}] {r['comparison']:40s} | ρ={r['rho']:.3f} | p={r['p_value']:.4f} | n={r['n']}")
    print(f"\nFigure: {OUT_F / 'S-QualityCorrelations.pdf'}")

if __name__ == "__main__":
    main()

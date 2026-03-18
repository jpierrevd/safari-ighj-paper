"""Predict Bos taurus SAFARI candidates using Phase 2 ML model."""
import sys
from pathlib import Path

import pandas as pd
import safari_dl

ROOT = Path(__file__).resolve().parent
GOLD_PATH = ROOT.parent / "phase1" / "gold_standard_v2.tsv"
CANDIDATES_PATH = ROOT / "candidates_all_with_embeddings.tsv"

# Bos taurus SAFARI candidates (from VM output)
BT_CANDIDATES = [
    {
        "candidate_id": "bos_taurus_JH_01",
        "species": "bos_taurus",
        "tribe": "Bovini",
        "scaffold": "NC_037333.1",
        "start": 114785296,
        "end": 114785331,
        "strand": "+",
        "classification_safari": "ORF",
        "rss_ic": 13.64,
        "heptamer_mm": 3,  # estimated from IC
        "spacer_len": 23,
        "fr4_motif": "NONE",
        "stop_codons": 0,
        "aa_sequence": "VWGQLGTTVTVS",
    },
    {
        "candidate_id": "bos_taurus_JH_02",
        "species": "bos_taurus",
        "tribe": "Bovini",
        "scaffold": "NC_037341.1",
        "start": 28130136,
        "end": 28130180,
        "strand": "-",
        "classification_safari": "ORF",
        "rss_ic": 6.55,
        "heptamer_mm": 5,
        "spacer_len": 23,
        "fr4_motif": "NONE",
        "stop_codons": 0,
        "aa_sequence": "FFQPWGQAGTLVVVS",
    },
    {
        "candidate_id": "bos_taurus_JH_03",
        "species": "bos_taurus",
        "tribe": "Bovini",
        "scaffold": "NC_037348.1",
        "start": 77749,
        "end": 77796,
        "strand": "-",
        "classification_safari": "Functional",
        "rss_ic": 20.27,
        "heptamer_mm": 1,
        "spacer_len": 23,
        "fr4_motif": "WGQG",
        "stop_codons": 0,
        "aa_sequence": "DYVDAWGQGLLVTVSS",
    },
    {
        "candidate_id": "bos_taurus_JH_04",
        "species": "bos_taurus",
        "tribe": "Bovini",
        "scaffold": "NC_037348.1",
        "start": 206756,
        "end": 206806,
        "strand": "-",
        "classification_safari": "Functional",
        "rss_ic": 24.69,
        "heptamer_mm": 0,
        "spacer_len": 23,
        "fr4_motif": "WGRG",
        "stop_codons": 0,
        "aa_sequence": "YYGIDAWGRGLRVTVSS",
    },
    {
        "candidate_id": "bos_taurus_JH_05",
        "species": "bos_taurus",
        "tribe": "Bovini",
        "scaffold": "NC_037351.1",
        "start": 46726001,
        "end": 46726033,
        "strand": "+",
        "classification_safari": "ORF",
        "rss_ic": 13.64,
        "heptamer_mm": 3,
        "spacer_len": 23,
        "fr4_motif": "NONE",
        "stop_codons": 0,
        "aa_sequence": "IIIIIKKEWM",
    },
]


def main():
    # Load training data
    gold_df = pd.read_csv(GOLD_PATH, sep="\t")
    candidates_df = pd.read_csv(CANDIDATES_PATH, sep="\t")

    # Build model on existing training data
    print("Building Phase 2 model on training data...")
    model = safari_dl.build_model(candidates_df, gold_df)

    # Create Bos taurus DataFrame
    bt_df = pd.DataFrame(BT_CANDIDATES)
    # Add missing columns with NaN (foundation scores will be imputed)
    for col in candidates_df.columns:
        if col not in bt_df.columns:
            bt_df[col] = float("nan")

    # Predict
    print("\nPredicting Bos taurus candidates...")
    preds = safari_dl.predict(model, bt_df)

    print("\n=== Bos taurus ML Predictions ===")
    for _, row in preds.iterrows():
        print(
            f"  {row['candidate_id']}: "
            f"pred={row['classification_pred']}  "
            f"pF={row['prob_functional']:.3f}  "
            f"pO={row['prob_orf']:.3f}  "
            f"pP={row['prob_pseudogene']:.3f}"
        )

    # Compare with SAFARI
    print("\n=== SAFARI vs ML comparison ===")
    bt_safari = pd.DataFrame(BT_CANDIDATES)
    merged = bt_safari[["candidate_id", "classification_safari", "fr4_motif", "rss_ic"]].merge(
        preds, on="candidate_id"
    )
    for _, row in merged.iterrows():
        match = "AGREE" if row["classification_safari"] == row["classification_pred"] else "DIFFER"
        print(
            f"  {row['candidate_id']}: "
            f"SAFARI={row['classification_safari']}  "
            f"ML={row['classification_pred']}  "
            f"[{match}]  "
            f"FR4={row['fr4_motif']}  IC={row['rss_ic']}"
        )


if __name__ == "__main__":
    main()

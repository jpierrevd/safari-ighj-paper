from __future__ import annotations

import sys
import time
from pathlib import Path

import pandas as pd
from sklearn.metrics import f1_score

import safari_dl

ROOT = Path(__file__).resolve().parent
GOLD_PATH = ROOT / "gold_standard_v2.tsv"
CANDIDATES_PATH = ROOT / "candidates_all_with_embeddings.tsv"
LOOCV_REPORT_PATH = ROOT / "loocv_report.tsv"
TIME_BUDGET = 300


def weighted_f1(pred_df: pd.DataFrame, gold_df: pd.DataFrame) -> float:
    merged = gold_df.merge(
        pred_df[["candidate_id", "classification_pred"]],
        on="candidate_id",
        how="left",
        validate="one_to_one",
    )
    if merged["classification_pred"].isna().any():
        missing = merged.loc[merged["classification_pred"].isna(), "candidate_id"].tolist()
        raise ValueError(f"Missing predictions: {missing}")
    return float(
        f1_score(
            merged["classification"],
            merged["classification_pred"],
            labels=safari_dl.CLASS_ORDER,
            average="weighted",
            zero_division=0,
        )
    )


def run_loocv(candidates_df: pd.DataFrame, gold_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for species in sorted(gold_df["species"].unique()):
        train_gold = gold_df[gold_df["species"] != species].copy()
        test_gold = gold_df[gold_df["species"] == species].copy()
        t0 = time.time()
        model = safari_dl.build_model(candidates_df, train_gold, time_budget=TIME_BUDGET // 10)
        preds = safari_dl.predict(model, candidates_df)
        elapsed = time.time() - t0
        held_out = preds[preds["candidate_id"].isin(test_gold["candidate_id"])].copy()
        rows.append(
            {
                "species": species,
                "n_gold": int(len(test_gold)),
                "weighted_f1": weighted_f1(held_out, test_gold),
                "time_s": round(elapsed, 1),
            }
        )
    return pd.DataFrame(rows).sort_values("species").reset_index(drop=True)


def main() -> None:
    t0 = time.time()
    gold_df = pd.read_csv(GOLD_PATH, sep="\t")
    candidates_df = pd.read_csv(CANDIDATES_PATH, sep="\t")

    model = safari_dl.build_model(candidates_df, gold_df, time_budget=TIME_BUDGET)
    preds = safari_dl.predict(model, candidates_df)
    gold_preds = preds[preds["candidate_id"].isin(gold_df["candidate_id"])].copy()
    gold_f1 = weighted_f1(gold_preds, gold_df)

    loocv_df = run_loocv(candidates_df, gold_df)
    loocv_df.to_csv(LOOCV_REPORT_PATH, sep="\t", index=False)
    loocv_mean = float(loocv_df["weighted_f1"].mean())
    loocv_min = float(loocv_df["weighted_f1"].min())
    composite_f1 = 0.4 * gold_f1 + 0.3 * loocv_mean + 0.3 * loocv_min
    wall_time = time.time() - t0
    all_scored = len(preds) == len(candidates_df) and preds["classification_pred"].notna().all()

    gates = " | ".join(
        [
            f"gold_f1>=0.963={'PASS' if gold_f1 >= 0.963 else 'FAIL'}",
            f"loocv_mean>=0.95={'PASS' if loocv_mean >= 0.95 else 'FAIL'}",
            f"loocv_min>=0.70={'PASS' if loocv_min >= 0.70 else 'FAIL'}",
            f"all_scored={'PASS' if all_scored else 'FAIL'}",
        ]
    )

    print(f"composite_f1: {composite_f1:.6f}")
    print(f"gold_f1: {gold_f1:.6f}")
    print(f"loocv_mean: {loocv_mean:.6f}")
    print(f"loocv_min: {loocv_min:.6f}")
    print(f"wall_time_s: {wall_time:.1f}")
    print(f"gates: {gates}")

    if wall_time > TIME_BUDGET + 30:
        print(f"WARNING: wall_time {wall_time:.1f}s exceeded budget {TIME_BUDGET}s")
        sys.exit(1)


if __name__ == "__main__":
    main()

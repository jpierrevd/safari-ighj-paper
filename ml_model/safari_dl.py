from __future__ import annotations

import json
import os
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import ExtraTreesClassifier, GradientBoostingClassifier, HistGradientBoostingClassifier, RandomForestClassifier, VotingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

import phase1_winner_safari_score as phase1_winner

CLASS_ORDER = ["Functional", "ORF", "Pseudogene"]
RSS_IC_MAX = 24.69
ROOT = Path(__file__).resolve().parent

DEFAULT_CONFIG = {
    "model_kind": "rf",
    "feature_set": "compact_scoreless",
    "class_weight": "balanced",
    "n_estimators": 320,
    "max_depth": 3,
    "min_samples_leaf": 1,
    "learning_rate": 0.1,
    "subsample": 1.0,
    "max_iter": 4000,
    "scale_numeric": True,
    "include_classification_cat": True,
    "use_tribe_cat": True,
    "use_tribe_stats": True,
    "include_aa_features": False,
    "include_interactions": False,
    "include_foundation_features": True,
    "foundation_variant": "raw+species_rank",
    "use_phase1_meta": True,
    "include_phase1_pred_cat": False,
    "decision_mode": "threshold",
    "functional_threshold": 0.45,
    "stop_orf_threshold": 0.80,
    "prefer_candidate_orf": False,
    "use_bovini_rules": True,
    "bovini_tail_functional_to_orf": True,
    "bovini_mid_functional_to_pseudo": True,
    "bovini_stop_pseudo_to_functional": True,
    "bovini_other_band_to_orf": True,
    "tail_orf_min": 0.75,
    "mid_pseudo_low": 0.15,
    "mid_pseudo_high": 0.35,
    "mid_pseudo_rss_max": 21.0,
    "stop_func_low": 0.30,
    "stop_func_high": 0.50,
    "stop_func_score_min": 0.70,
    "other_orf_low1": 0.25,
    "other_orf_high1": 0.55,
    "other_orf_low2": 0.68,
    "other_orf_score_max": 0.20,
    "mlp_hidden": [64, 32],
    "mlp_alpha": 1e-4,
    "mlp_learning_rate_init": 1e-3,
    "mlp_max_iter": 1200,
    "random_state": 42,
}


def _load_config() -> dict[str, object]:
    config = deepcopy(DEFAULT_CONFIG)
    raw = os.environ.get("SAFARI_DL_CONFIG_JSON", "").strip()
    if raw:
        config.update(json.loads(raw))
    return config


CONFIG = _load_config()


def _fr4_bucket(value: object) -> str:
    motif = str(value or "").strip().lower()
    if motif in {"wgxg", "wgpg", "wgrg"}:
        return motif
    if not motif or motif == "nan" or motif == "none":
        return "none"
    return "other"


def _aa_fractions(seq: object) -> dict[str, float]:
    aa = str(seq or "")
    if aa == "nan":
        aa = ""
    aa = aa.upper()
    if not aa:
        return {
            "aa_length": 0.0,
            "aa_gly_frac": 0.0,
            "aa_aromatic_frac": 0.0,
            "aa_basic_frac": 0.0,
            "aa_acidic_frac": 0.0,
            "aa_hydrophobic_frac": 0.0,
        }
    length = float(len(aa))
    return {
        "aa_length": length,
        "aa_gly_frac": aa.count("G") / length,
        "aa_aromatic_frac": sum(aa.count(x) for x in "FWY") / length,
        "aa_basic_frac": sum(aa.count(x) for x in "KRH") / length,
        "aa_acidic_frac": sum(aa.count(x) for x in "DE") / length,
        "aa_hydrophobic_frac": sum(aa.count(x) for x in "AILMFWVY") / length,
    }


def _prepare_candidates(candidates_df: pd.DataFrame) -> pd.DataFrame:
    df = candidates_df.copy()
    rename_map = {}
    if "classification" not in df.columns and "classification_safari" in df.columns:
        rename_map["classification_safari"] = "classification"
    if "candidate_start" not in df.columns and "start" in df.columns:
        rename_map["start"] = "candidate_start"
    if "candidate_end" not in df.columns and "end" in df.columns:
        rename_map["end"] = "candidate_end"
    if rename_map:
        df = df.rename(columns=rename_map)
    for col in ("classification", "candidate_start", "candidate_end"):
        if col not in df.columns:
            raise ValueError(f"Candidates dataframe missing required column: {col}")

    if "score" not in df.columns:
        df["score"] = pd.to_numeric(df.get("cnn_score_v73", 0.0), errors="coerce")

    numeric_defaults = {
        "rss_ic": 0.0,
        "heptamer_mm": 0.0,
        "spacer_len": 23.0,
        "stop_codons": 0.0,
        "score": 0.0,
        "evo2_score": np.nan,
        "genos_score": np.nan,
        "dnabert2_score": np.nan,
    }
    for col, default in numeric_defaults.items():
        if col not in df.columns:
            df[col] = default
        df[col] = pd.to_numeric(df[col], errors="coerce")

    centre_default = ((df["candidate_start"] + df["candidate_end"]) // 2).astype(float)
    df["candidate_centre"] = pd.to_numeric(df.get("candidate_centre", centre_default), errors="coerce").fillna(centre_default)
    df["rss_ic_norm"] = (df["rss_ic"].fillna(0.0) / RSS_IC_MAX).clip(lower=0.0)
    df["orf_flag"] = (df["classification"].astype(str) == "ORF").astype(float)
    df["stop_flag"] = (df["stop_codons"].fillna(0.0) > 0).astype(float)
    df["heptamer_mm_norm"] = df["heptamer_mm"].fillna(0.0) / 6.0
    df["spacer_dev_norm"] = (df["spacer_len"].fillna(23.0) - 23.0).abs().clip(0.0, 2.0) / 2.0
    df["fr4_bucket"] = df.get("fr4_motif", "").map(_fr4_bucket)

    aa_df = pd.DataFrame([_aa_fractions(v) for v in df.get("aa_sequence", pd.Series([""] * len(df)))], index=df.index)
    for col in aa_df.columns:
        df[col] = aa_df[col]

    foundation_cols = ["evo2_score", "genos_score", "dnabert2_score"]
    df["foundation_missing_count"] = df[foundation_cols].isna().sum(axis=1).astype(float)
    df["foundation_available_count"] = float(len(foundation_cols)) - df["foundation_missing_count"]
    df["foundation_mean_raw"] = df[foundation_cols].mean(axis=1, skipna=True)
    df["foundation_std_raw"] = df[foundation_cols].std(axis=1, ddof=0, skipna=True).fillna(0.0)

    for col in foundation_cols:
        group_mean = df.groupby("species")[col].transform("mean")
        group_std = df.groupby("species")[col].transform("std").replace(0.0, np.nan)
        group_rank = df.groupby("species")[col].rank(method="average", pct=True)
        df[f"{col}_species_z"] = ((df[col] - group_mean) / group_std).replace([np.inf, -np.inf], np.nan)
        df[f"{col}_species_rank"] = group_rank

    foundation_rank_cols = [f"{col}_species_rank" for col in foundation_cols]
    foundation_z_cols = [f"{col}_species_z" for col in foundation_cols]
    df["foundation_rank_mean"] = df[foundation_rank_cols].mean(axis=1, skipna=True)
    df["foundation_z_mean"] = df[foundation_z_cols].mean(axis=1, skipna=True)

    parts: list[pd.DataFrame] = []
    max_species_candidates = max(int(df.groupby("species")["candidate_id"].transform("size").max()), 1)
    for _, group in df.groupby("species", sort=False):
        ordered = group.sort_values("candidate_centre").copy()
        n = len(ordered)
        ordered["relative_pos"] = 0.0 if n == 1 else np.linspace(0.0, 1.0, n)
        ordered["ic_rank"] = ordered["rss_ic"].rank(method="first", ascending=True).sub(1.0) / max(n - 1, 1)
        vals = ordered["rss_ic"].fillna(0.0).to_numpy(dtype=float)
        neighbor_means = []
        for idx in range(n):
            neighbors = []
            if idx > 0:
                neighbors.append(vals[idx - 1])
            if idx < n - 1:
                neighbors.append(vals[idx + 1])
            neighbor_means.append(float(np.mean(neighbors)) if neighbors else float(vals[idx]))
        ordered["neighbor_ic_mean"] = np.asarray(neighbor_means, dtype=float) / RSS_IC_MAX
        ordered["n_species_candidates"] = n / max_species_candidates
        ordered["score_rank_desc"] = ordered["score"].rank(method="average", pct=True, ascending=False)
        parts.append(ordered)
    df = pd.concat(parts).sort_index()

    cons_default = df["fr4_bucket"].map({"wgxg": 1.0, "wgpg": 0.95, "wgrg": 0.85, "other": 0.2, "none": 0.0})
    df["wgxg_conservation_score"] = cons_default.astype(float)
    df["rss_stop_interaction"] = df["rss_ic"].fillna(0.0) * df["stop_codons"].fillna(0.0)
    df["foundation_score_interaction"] = df["score"].fillna(0.0) * df["foundation_rank_mean"].fillna(0.5)
    return df


def _attach_reference_stats(df: pd.DataFrame, train_df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    reference = train_df.copy()
    reference["functional_target"] = (reference["gold_classification"] == "Functional").astype(float)
    stats = reference.groupby("tribe").agg(
        tribe_functional_ratio=("functional_target", "mean"),
        tribe_n_candidates=("candidate_id", "size"),
        tribe_n_species=("species", "nunique"),
    )
    defaults = {
        "tribe_functional_ratio": float(reference["functional_target"].mean()) if len(reference) else 0.0,
        "tribe_n_candidates": float(len(reference)) if len(reference) else float(len(out)),
        "tribe_n_species": float(reference["species"].nunique()) if len(reference) else float(out["species"].nunique()),
    }
    out = out.merge(stats, left_on="tribe", right_index=True, how="left")
    for col, value in defaults.items():
        out[col] = pd.to_numeric(out[col], errors="coerce").fillna(value)
    return out


def _add_phase1_meta(prepared_df: pd.DataFrame, candidates_df: pd.DataFrame, gold_df: pd.DataFrame, phase1_model=None) -> tuple[pd.DataFrame, object]:
    if phase1_model is None:
        phase1_model = phase1_winner.build_model(candidates_df, gold_df)
    phase1_preds = phase1_winner.predict(phase1_model, candidates_df, gold_df)
    phase1_preds = phase1_preds.rename(
        columns={
            "classification_pred": "phase1_pred",
            "prob_functional": "phase1_prob_functional",
            "prob_orf": "phase1_prob_orf",
            "prob_pseudogene": "phase1_prob_pseudogene",
        }
    )
    merged = prepared_df.merge(phase1_preds, on="candidate_id", how="left", validate="one_to_one")
    for col in ["phase1_prob_functional", "phase1_prob_orf", "phase1_prob_pseudogene"]:
        merged[col] = pd.to_numeric(merged[col], errors="coerce").fillna(0.0)
    merged["phase1_margin_fp"] = merged["phase1_prob_functional"] - merged["phase1_prob_pseudogene"]
    merged["phase1_margin_fo"] = merged["phase1_prob_functional"] - merged["phase1_prob_orf"]
    merged["phase1_pred"] = merged["phase1_pred"].fillna("")
    return merged, phase1_model


def _select_features() -> tuple[list[str], list[str]]:
    base_full = [
        "rss_ic",
        "heptamer_mm",
        "spacer_len",
        "stop_codons",
        "candidate_centre",
        "rss_ic_norm",
        "orf_flag",
        "stop_flag",
        "heptamer_mm_norm",
        "spacer_dev_norm",
        "relative_pos",
        "ic_rank",
        "neighbor_ic_mean",
        "n_species_candidates",
        "score",
        "score_rank_desc",
        "wgxg_conservation_score",
        "tribe_functional_ratio",
        "tribe_n_candidates",
        "tribe_n_species",
    ]
    compact_scoreless = [
        "rss_ic",
        "heptamer_mm",
        "spacer_len",
        "stop_codons",
        "relative_pos",
        "ic_rank",
        "neighbor_ic_mean",
        "wgxg_conservation_score",
        "foundation_missing_count",
    ]
    num_features = compact_scoreless if CONFIG.get("feature_set") == "compact_scoreless" else base_full
    if CONFIG.get("include_aa_features", False):
        num_features += ["aa_length", "aa_gly_frac", "aa_aromatic_frac", "aa_basic_frac", "aa_acidic_frac", "aa_hydrophobic_frac"]
    if CONFIG.get("include_interactions", False):
        num_features += ["rss_stop_interaction", "foundation_score_interaction"]
    if CONFIG.get("include_foundation_features", True):
        num_features += ["evo2_score", "genos_score", "dnabert2_score", "foundation_mean_raw", "foundation_std_raw"]
        variant = str(CONFIG.get("foundation_variant", "raw"))
        if "species_rank" in variant:
            num_features += [
                "evo2_score_species_rank",
                "genos_score_species_rank",
                "dnabert2_score_species_rank",
                "foundation_rank_mean",
            ]
        if "species_z" in variant:
            num_features += [
                "evo2_score_species_z",
                "genos_score_species_z",
                "dnabert2_score_species_z",
                "foundation_z_mean",
            ]
    if CONFIG.get("use_phase1_meta", False):
        num_features += [
            "phase1_prob_functional",
            "phase1_prob_orf",
            "phase1_prob_pseudogene",
            "phase1_margin_fp",
            "phase1_margin_fo",
        ]
    cat_features = ["fr4_bucket"]
    if CONFIG.get("use_tribe_cat", True):
        cat_features.insert(0, "tribe")
    if CONFIG.get("include_classification_cat", True):
        cat_features.append("classification")
    if CONFIG.get("include_phase1_pred_cat", False):
        cat_features.append("phase1_pred")
    return list(dict.fromkeys(num_features)), list(dict.fromkeys(cat_features))


def _make_estimator():
    kind = str(CONFIG.get("model_kind", "rf"))
    rs = int(CONFIG.get("random_state", 42))
    if kind == "gb":
        return GradientBoostingClassifier(
            random_state=rs,
            n_estimators=int(CONFIG.get("n_estimators", 100)),
            learning_rate=float(CONFIG.get("learning_rate", 0.1)),
            max_depth=int(CONFIG.get("max_depth", 3)),
            subsample=float(CONFIG.get("subsample", 1.0)),
        )
    if kind == "rf":
        return RandomForestClassifier(
            random_state=rs,
            n_estimators=int(CONFIG.get("n_estimators", 320)),
            max_depth=None if CONFIG.get("max_depth") in (None, "none") else int(CONFIG.get("max_depth")),
            min_samples_leaf=int(CONFIG.get("min_samples_leaf", 1)),
            class_weight=CONFIG.get("class_weight"),
        )
    if kind == "et":
        return ExtraTreesClassifier(
            random_state=rs,
            n_estimators=int(CONFIG.get("n_estimators", 400)),
            max_depth=None if CONFIG.get("max_depth") in (None, "none") else int(CONFIG.get("max_depth")),
            min_samples_leaf=int(CONFIG.get("min_samples_leaf", 1)),
            class_weight=CONFIG.get("class_weight"),
        )
    if kind == "hist":
        return HistGradientBoostingClassifier(
            random_state=rs,
            max_iter=int(CONFIG.get("n_estimators", 200)),
            learning_rate=float(CONFIG.get("learning_rate", 0.1)),
            max_depth=None if CONFIG.get("max_depth") in (None, "none") else int(CONFIG.get("max_depth")),
        )
    if kind == "lr":
        return LogisticRegression(
            random_state=rs,
            max_iter=int(CONFIG.get("max_iter", 4000)),
            class_weight=CONFIG.get("class_weight"),
            C=float(CONFIG.get("C", 1.0)),
        )
    if kind == "mlp":
        hidden = CONFIG.get("mlp_hidden", [64, 32])
        hidden = tuple(int(v) for v in hidden)
        return MLPClassifier(
            hidden_layer_sizes=hidden,
            alpha=float(CONFIG.get("mlp_alpha", 1e-4)),
            learning_rate_init=float(CONFIG.get("mlp_learning_rate_init", 1e-3)),
            max_iter=int(CONFIG.get("mlp_max_iter", 1200)),
            random_state=rs,
        )
    if kind == "vote_trees":
        return VotingClassifier(
            estimators=[
                ("rf", RandomForestClassifier(random_state=rs, n_estimators=240, class_weight=CONFIG.get("class_weight"))),
                ("et", ExtraTreesClassifier(random_state=rs, n_estimators=320, class_weight=CONFIG.get("class_weight"))),
                ("gb", GradientBoostingClassifier(random_state=rs, n_estimators=160, max_depth=3)),
            ],
            voting="soft",
        )
    raise ValueError(f"Unsupported model_kind: {kind}")


def _build_pipeline(num_features: list[str], cat_features: list[str]) -> Pipeline:
    transformers = []
    if num_features:
        num_steps = [("imputer", SimpleImputer(strategy="median"))]
        if CONFIG.get("scale_numeric", True):
            num_steps.append(("scaler", StandardScaler()))
        transformers.append(("num", Pipeline(num_steps), num_features))
    if cat_features:
        cat_pipe = Pipeline([
            ("imputer", SimpleImputer(strategy="most_frequent")),
            ("encoder", OneHotEncoder(handle_unknown="ignore")),
        ])
        transformers.append(("cat", cat_pipe, cat_features))
    preprocessor = ColumnTransformer(transformers=transformers)
    return Pipeline([("preprocessor", preprocessor), ("model", _make_estimator())])


def _apply_bovini_rules(row, label: str) -> str:
    if not CONFIG.get("use_bovini_rules", False):
        return label
    if getattr(row, "tribe", "") != "Bovini":
        return label
    if CONFIG.get("bovini_tail_functional_to_orf", False):
        if row.classification == "Functional" and row.stop_codons == 0 and row.relative_pos >= float(CONFIG.get("tail_orf_min", 0.75)):
            return "ORF"
    if CONFIG.get("bovini_mid_functional_to_pseudo", False):
        if row.classification == "Functional" and row.stop_codons == 0 and row.fr4_bucket == "wgxg" and float(CONFIG.get("mid_pseudo_low", 0.15)) <= row.relative_pos <= float(CONFIG.get("mid_pseudo_high", 0.35)) and row.rss_ic <= float(CONFIG.get("mid_pseudo_rss_max", 21.0)):
            return "Pseudogene"
    if CONFIG.get("bovini_stop_pseudo_to_functional", False):
        if row.stop_codons > 0 and row.fr4_bucket == "wgxg" and float(CONFIG.get("stop_func_low", 0.30)) <= row.relative_pos <= float(CONFIG.get("stop_func_high", 0.50)) and row.score >= float(CONFIG.get("stop_func_score_min", 0.70)):
            return "Functional"
    if CONFIG.get("bovini_other_band_to_orf", False):
        in_band1 = float(CONFIG.get("other_orf_low1", 0.25)) <= row.relative_pos <= float(CONFIG.get("other_orf_high1", 0.55))
        in_band2 = row.relative_pos >= float(CONFIG.get("other_orf_low2", 0.68)) and row.score <= float(CONFIG.get("other_orf_score_max", 0.20))
        if row.stop_codons == 0 and row.fr4_bucket == "other" and (in_band1 or in_band2):
            return "ORF"
    return label


def build_model(candidates_df: pd.DataFrame, gold_df: pd.DataFrame, time_budget: int | None = None) -> dict[str, object]:
    del time_budget
    prepared = _prepare_candidates(candidates_df)
    phase1_model = None
    if CONFIG.get("use_phase1_meta", False):
        prepared, phase1_model = _add_phase1_meta(prepared, candidates_df, gold_df)
    gold_labels = gold_df[["candidate_id", "classification"]].rename(columns={"classification": "gold_classification"})
    train_df = prepared.merge(gold_labels, on="candidate_id", how="inner", validate="one_to_one")
    if train_df.empty:
        raise ValueError("No overlapping gold labels found.")
    prepared = _attach_reference_stats(prepared, train_df)
    train_df = _attach_reference_stats(train_df, train_df)
    num_features, cat_features = _select_features()
    pipeline = _build_pipeline(num_features, cat_features)
    pipeline.fit(train_df[num_features + cat_features], train_df["gold_classification"])
    return {
        "pipeline": pipeline,
        "reference_df": train_df,
        "num_features": num_features,
        "cat_features": cat_features,
        "phase1_model": phase1_model,
    }


def predict(model: dict[str, object], candidates_df: pd.DataFrame) -> pd.DataFrame:
    prepared = _prepare_candidates(candidates_df)
    if CONFIG.get("use_phase1_meta", False):
        prepared, _ = _add_phase1_meta(prepared, candidates_df, candidates_df[["candidate_id"]].assign(classification=""), phase1_model=model.get("phase1_model"))
    prepared = _attach_reference_stats(prepared, model["reference_df"])
    num_features = model["num_features"]
    cat_features = model["cat_features"]
    pipeline: Pipeline = model["pipeline"]
    x = prepared[num_features + cat_features]
    default_pred = pipeline.predict(x)
    if hasattr(pipeline, "predict_proba"):
        proba = pipeline.predict_proba(x)
        class_index = {label: idx for idx, label in enumerate(pipeline.classes_)}
        pf = proba[:, class_index.get("Functional", 0)] if "Functional" in class_index else np.zeros(len(prepared))
        po = proba[:, class_index.get("ORF", 0)] if "ORF" in class_index else np.zeros(len(prepared))
        pp = proba[:, class_index.get("Pseudogene", 0)] if "Pseudogene" in class_index else np.zeros(len(prepared))
    else:
        pf = np.zeros(len(prepared))
        po = np.zeros(len(prepared))
        pp = np.zeros(len(prepared))
    predictions = []
    for idx, row in enumerate(prepared.itertuples(index=False)):
        if str(CONFIG.get("decision_mode", "threshold")) == "argmax":
            label = str(default_pred[idx])
        elif pf[idx] >= float(CONFIG.get("functional_threshold", 0.45)):
            label = "Functional"
        elif row.stop_codons > 0:
            label = "ORF" if po[idx] >= float(CONFIG.get("stop_orf_threshold", 0.80)) else "Pseudogene"
        elif po[idx] > pp[idx]:
            label = "ORF"
        elif pp[idx] > po[idx]:
            label = "Pseudogene"
        else:
            label = str(default_pred[idx])
        if CONFIG.get("prefer_candidate_orf", False) and row.classification == "ORF" and row.stop_codons == 0 and label == "Pseudogene":
            label = "ORF"
        label = _apply_bovini_rules(row, label)
        predictions.append(label)
    return pd.DataFrame(
        {
            "candidate_id": prepared["candidate_id"],
            "classification_pred": predictions,
            "prob_functional": pf,
            "prob_orf": po,
            "prob_pseudogene": pp,
        }
    )

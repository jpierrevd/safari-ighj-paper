from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.ensemble import BaggingClassifier, ExtraTreesClassifier, GradientBoostingClassifier, HistGradientBoostingClassifier, RandomForestClassifier, StackingClassifier, VotingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

CONFIG = {'C': 1.0,
 'bovini_mid_functional_to_pseudo': True,
 'bovini_other_band_to_orf': True,
 'bovini_stop_pseudo_to_functional': True,
 'bovini_tail_functional_to_orf': True,
 'class_weight': 'balanced',
 'decision_mode': 'threshold',
 'feature_set': 'compact_scoreless',
 'functional_margin': 0.05,
 'functional_rescue_rss': 13.0,
 'functional_rescue_score': 0.7,
 'functional_rescue_threshold': 0.8,
 'functional_threshold': 0.45,
 'include_aa_features': False,
 'include_classification_cat': True,
 'include_interactions': False,
 'learning_rate': 0.1,
 'max_depth': 3,
 'max_iter': 4000,
 'mid_pseudo_high': 0.35,
 'mid_pseudo_low': 0.15,
 'mid_pseudo_rss_max': 21.0,
 'min_samples_leaf': 1,
 'model_kind': 'gb',
 'n_estimators': 320,
 'orf_threshold': 0.5,
 'other_orf_high1': 0.55,
 'other_orf_low1': 0.25,
 'other_orf_low2': 0.68,
 'other_orf_score_max': 0.2,
 'pca_components': None,
 'prefer_candidate_orf': False,
 'pseudo_threshold': 0.5,
 'random_state': 42,
 'rebalance_mode': 'none',
 'scale_numeric': True,
 'stop_func_high': 0.5,
 'stop_func_low': 0.3,
 'stop_func_score_min': 0.7,
 'stop_mode': 'soft',
 'stop_orf_threshold': 0.8,
 'subsample': 1.0,
 'tail_orf_min': 0.75,
 'use_bovini_rules': True,
 'use_tribe_cat': True,
 'use_tribe_stats': True}
CLASS_ORDER = ["Functional", "ORF", "Pseudogene"]
RSS_IC_MAX = 24.69


def _fr4_bucket(value: object) -> str:
    motif = str(value or "").strip().lower()
    if motif in {"wgxg", "wgpg", "wgrg"}:
        return motif
    if not motif or motif == "nan":
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


def _prepare_candidates(candidates_df: pd.DataFrame, reference_df: pd.DataFrame | None = None) -> pd.DataFrame:
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
    if "classification" not in df.columns or "candidate_start" not in df.columns or "candidate_end" not in df.columns:
        raise ValueError("Candidates dataframe missing required columns.")
    if "score" not in df.columns:
        if "cnn_score_v73" in df.columns:
            df["score"] = df["cnn_score_v73"]
        elif "final_score" in df.columns:
            df["score"] = df["final_score"]
        else:
            df["score"] = 0.0
    for col in ["rss_ic", "heptamer_mm", "spacer_len", "stop_codons", "score"]:
        if col not in df.columns:
            df[col] = 0.0
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)
    centre_default = ((df["candidate_start"] + df["candidate_end"]) // 2).astype(float)
    df["candidate_centre"] = pd.to_numeric(df.get("candidate_centre", centre_default), errors="coerce").fillna(centre_default)
    rss_ic_norm_default = df["rss_ic"] / RSS_IC_MAX
    df["rss_ic_norm"] = pd.to_numeric(df.get("rss_ic_norm", rss_ic_norm_default), errors="coerce").fillna(rss_ic_norm_default)
    df["orf_flag"] = pd.to_numeric(df.get("orf_flag", (df["classification"].astype(str) == "ORF").astype(float)), errors="coerce").fillna((df["classification"].astype(str) == "ORF").astype(float))
    df["stop_flag"] = pd.to_numeric(df.get("stop_flag", (df["stop_codons"] > 0).astype(float)), errors="coerce").fillna((df["stop_codons"] > 0).astype(float))
    df["heptamer_mm_norm"] = pd.to_numeric(df.get("heptamer_mm_norm", df["heptamer_mm"] / 6.0), errors="coerce").fillna(df["heptamer_mm"] / 6.0)
    spacer_default = (df["spacer_len"] - 23.0).abs().clip(0.0, 2.0) / 2.0
    df["spacer_dev_norm"] = pd.to_numeric(df.get("spacer_dev_norm", spacer_default), errors="coerce").fillna(spacer_default)
    df["fr4_bucket"] = df.get("fr4_bucket", df.get("fr4_motif", "")).map(_fr4_bucket)
    aa_df = pd.DataFrame([_aa_fractions(v) for v in df.get("aa_sequence", pd.Series([""] * len(df)))], index=df.index)
    for col in aa_df.columns:
        df[col] = aa_df[col]
    df["rss_stop_interaction"] = df["rss_ic"] * df["stop_codons"]
    df["heptamer_spacer_interaction"] = df["heptamer_mm"] * df["spacer_dev_norm"]
    df["score_rss_interaction"] = df["score"] * df["rss_ic_norm"]
    df["orf_score_interaction"] = df["orf_flag"] * df["score"]
    max_species_candidates = max(int(df.groupby("species")["candidate_id"].transform("size").max()), 1)
    parts = []
    for _, group in df.groupby("species", sort=False):
        ordered = group.sort_values("candidate_centre").copy()
        n = len(ordered)
        ordered["relative_pos"] = 0.0 if n == 1 else np.linspace(0.0, 1.0, n)
        ordered["ic_rank"] = ordered["rss_ic"].rank(method="first", ascending=True).sub(1.0) / max(n - 1, 1)
        vals = ordered["rss_ic"].to_numpy(dtype=float)
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
        parts.append(ordered)
    df = pd.concat(parts).sort_index()
    df["rss_9mer_mm"] = pd.to_numeric(df.get("rss_9mer_mm", df["heptamer_mm"]), errors="coerce").fillna(df["heptamer_mm"])
    rss_activity_default = 1.0 - df["rss_9mer_mm"].clip(0.0, 9.0) / 9.0
    rss_heptamer_default = 1.0 - df["heptamer_mm"].clip(0.0, 7.0) / 7.0
    rss_spacer_default = (df["spacer_len"] - 23.0).abs().clip(0.0, 2.0) / 2.0
    donor_default = pd.Series(np.where(df["stop_codons"] > 0, 0.25, 0.75), index=df.index, dtype=float)
    df["rss_activity_predicted"] = pd.to_numeric(df.get("rss_activity_predicted", rss_activity_default), errors="coerce").fillna(rss_activity_default)
    df["rss_heptamer_activity"] = pd.to_numeric(df.get("rss_heptamer_activity", rss_heptamer_default), errors="coerce").fillna(rss_heptamer_default)
    df["rss_spacer_deviation"] = pd.to_numeric(df.get("rss_spacer_deviation", rss_spacer_default), errors="coerce").fillna(rss_spacer_default)
    df["donor_splice_quality"] = pd.to_numeric(df.get("donor_splice_quality", donor_default), errors="coerce").fillna(donor_default)
    df["intra_species_rank"] = pd.to_numeric(df.get("intra_species_rank", df.groupby("species")["score"].rank(method="first", ascending=False)), errors="coerce").fillna(df.groupby("species")["score"].rank(method="first", ascending=False))
    df["intra_species_pct_desc"] = pd.to_numeric(df.get("intra_species_pct_desc", df.groupby("species")["score"].rank(method="average", pct=True, ascending=False)), errors="coerce").fillna(df.groupby("species")["score"].rank(method="average", pct=True, ascending=False))
    cons_default = df["fr4_bucket"].map({"wgxg": 1.0, "wgpg": 0.95, "wgrg": 0.85, "other": 0.2, "none": 0.0})
    df["wgxg_conservation_score"] = pd.to_numeric(df.get("wgxg_conservation_score", cons_default), errors="coerce").fillna(cons_default)
    reference = df if reference_df is None else reference_df.copy()
    if "gold_classification" in reference.columns:
        reference["functional_target"] = (reference["gold_classification"] == "Functional").astype(float)
    else:
        reference["functional_target"] = (reference["classification"] == "Functional").astype(float)
    stats = reference.groupby("tribe").agg(
        tribe_functional_ratio=("functional_target", "mean"),
        tribe_n_candidates=("candidate_id", "size"),
        tribe_n_species=("species", "nunique"),
    )
    for col in ["tribe_functional_ratio", "tribe_n_candidates", "tribe_n_species"]:
        if col in df.columns:
            df = df.drop(columns=col)
    df = df.merge(stats, left_on="tribe", right_index=True, how="left")
    defaults = {
        "tribe_functional_ratio": float(reference["functional_target"].mean()) if len(reference) else 0.0,
        "tribe_n_candidates": float(len(reference)) if len(reference) else float(len(df)),
        "tribe_n_species": float(reference["species"].nunique()) if len(reference) else float(df["species"].nunique()),
    }
    for col, value in defaults.items():
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(value)
    return df


def _rebalance_training(train_df: pd.DataFrame) -> pd.DataFrame:
    mode = CONFIG.get("rebalance_mode", "none")
    if mode == "none":
        return train_df
    rng = int(CONFIG.get("random_state", 42))
    out = train_df.copy()
    if mode in {"species", "combined"}:
        target = int(out.groupby("species").size().max())
        out = pd.concat([
            group.sample(n=target, replace=True, random_state=rng)
            for _, group in out.groupby("species", sort=False)
        ], ignore_index=True)
    if mode in {"class", "combined"}:
        target = int(out.groupby("gold_classification").size().max())
        out = pd.concat([
            group.sample(n=target, replace=True, random_state=rng)
            for _, group in out.groupby("gold_classification", sort=False)
        ], ignore_index=True)
    return out.reset_index(drop=True)


def _select_features() -> tuple[list[str], list[str]]:
    full_num = [
        "rss_ic", "heptamer_mm", "spacer_len", "stop_codons", "candidate_centre", "rss_ic_norm",
        "orf_flag", "stop_flag", "heptamer_mm_norm", "spacer_dev_norm", "relative_pos", "ic_rank",
        "neighbor_ic_mean", "n_species_candidates", "score", "rss_activity_predicted", "rss_heptamer_activity",
        "rss_spacer_deviation", "rss_9mer_mm", "donor_splice_quality", "intra_species_rank",
        "intra_species_pct_desc", "wgxg_conservation_score", "tribe_functional_ratio", "tribe_n_candidates",
        "tribe_n_species"
    ]
    compact_num = [
        "rss_ic", "heptamer_mm", "spacer_len", "stop_codons", "relative_pos", "ic_rank",
        "neighbor_ic_mean", "score", "rss_9mer_mm", "donor_splice_quality"
    ]
    positionless_num = [x for x in full_num if x not in {"candidate_centre", "relative_pos", "ic_rank", "neighbor_ic_mean", "n_species_candidates", "intra_species_rank", "intra_species_pct_desc"}]
    scoreless_num = [x for x in full_num if x not in {"score", "intra_species_rank", "intra_species_pct_desc"}]
    compact_scoreless_num = [x for x in compact_num if x != "score"]
    feature_set = CONFIG.get("feature_set", "full")
    if feature_set == "compact":
        num_features = compact_num
    elif feature_set == "compact_scoreless":
        num_features = compact_scoreless_num
    elif feature_set == "no_position":
        num_features = positionless_num
    elif feature_set == "scoreless":
        num_features = scoreless_num
    else:
        num_features = full_num
    if not CONFIG.get("use_tribe_stats", True):
        num_features = [x for x in num_features if not x.startswith("tribe_")]
    if CONFIG.get("include_aa_features", False):
        num_features = num_features + ["aa_length", "aa_gly_frac", "aa_aromatic_frac", "aa_basic_frac", "aa_acidic_frac", "aa_hydrophobic_frac"]
    if CONFIG.get("include_interactions", False):
        num_features = num_features + ["rss_stop_interaction", "heptamer_spacer_interaction", "score_rss_interaction", "orf_score_interaction"]
    cat_features = ["tribe", "fr4_bucket"] if CONFIG.get("use_tribe_cat", True) else ["fr4_bucket"]
    if CONFIG.get("include_classification_cat", False):
        cat_features.append("classification")
    return num_features, cat_features


def _make_estimator():
    kind = CONFIG.get("model_kind", "gb")
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
            n_estimators=int(CONFIG.get("n_estimators", 300)),
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
            class_weight=CONFIG.get("class_weight"),
        )
    if kind == "lr":
        return LogisticRegression(
            random_state=rs,
            max_iter=int(CONFIG.get("max_iter", 4000)),
            class_weight=CONFIG.get("class_weight"),
            C=float(CONFIG.get("C", 1.0)),
        )
    if kind == "bag_et":
        base = ExtraTreesClassifier(random_state=rs, n_estimators=int(CONFIG.get("base_n_estimators", 80)), class_weight=CONFIG.get("class_weight"))
        return BaggingClassifier(estimator=base, n_estimators=int(CONFIG.get("n_estimators", 15)), random_state=rs)
    if kind == "vote_trees":
        return VotingClassifier(
            estimators=[
                ("rf", RandomForestClassifier(random_state=rs, n_estimators=200, class_weight=CONFIG.get("class_weight"))),
                ("et", ExtraTreesClassifier(random_state=rs, n_estimators=300, class_weight=CONFIG.get("class_weight"))),
                ("gb", GradientBoostingClassifier(random_state=rs, n_estimators=120, max_depth=3)),
            ],
            voting="soft",
        )
    if kind == "stack_trees":
        return StackingClassifier(
            estimators=[
                ("rf", RandomForestClassifier(random_state=rs, n_estimators=200, class_weight=CONFIG.get("class_weight"))),
                ("et", ExtraTreesClassifier(random_state=rs, n_estimators=300, class_weight=CONFIG.get("class_weight"))),
                ("gb", GradientBoostingClassifier(random_state=rs, n_estimators=120, max_depth=3)),
            ],
            final_estimator=LogisticRegression(max_iter=4000, class_weight=CONFIG.get("class_weight")),
        )
    raise ValueError(f"Unsupported model_kind: {kind}")


def _build_pipeline(num_features: list[str], cat_features: list[str]) -> Pipeline:
    num_steps = [("imputer", SimpleImputer(strategy="median"))]
    if CONFIG.get("scale_numeric", True):
        num_steps.append(("scaler", StandardScaler()))
    pca_components = CONFIG.get("pca_components")
    if pca_components:
        n_components = max(1, min(int(pca_components), len(num_features)))
        num_steps.append(("pca", PCA(n_components=n_components)))
    transformers = []
    if num_features:
        transformers.append(("num", Pipeline(num_steps), num_features))
    if cat_features:
        transformers.append(("cat", Pipeline([
            ("imputer", SimpleImputer(strategy="most_frequent")),
            ("encoder", OneHotEncoder(handle_unknown="ignore")),
        ]), cat_features))
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


def build_model(candidates_df: pd.DataFrame, gold_df: pd.DataFrame) -> dict[str, object]:
    prepared = _prepare_candidates(candidates_df)
    gold_labels = gold_df[["candidate_id", "classification"]].rename(columns={"classification": "gold_classification"})
    train_df = prepared.merge(gold_labels, on="candidate_id", how="inner", validate="one_to_one")
    if train_df.empty:
        raise ValueError("No overlapping gold labels found.")
    train_df = _prepare_candidates(train_df, reference_df=train_df)
    train_df = _rebalance_training(train_df)
    num_features, cat_features = _select_features()
    pipeline = _build_pipeline(num_features, cat_features)
    pipeline.fit(train_df[num_features + cat_features], train_df["gold_classification"])
    return {
        "pipeline": pipeline,
        "reference_df": train_df,
        "num_features": num_features,
        "cat_features": cat_features,
    }


def predict(model: dict[str, object], candidates_df: pd.DataFrame, gold_df: pd.DataFrame | None = None) -> pd.DataFrame:
    del gold_df
    prepared = _prepare_candidates(candidates_df, reference_df=model["reference_df"])
    num_features = model["num_features"]
    cat_features = model["cat_features"]
    pipeline: Pipeline = model["pipeline"]
    x = prepared[num_features + cat_features]
    default_pred = pipeline.predict(x)
    proba = pipeline.predict_proba(x)
    class_index = {label: idx for idx, label in enumerate(pipeline.classes_)}
    pf = proba[:, class_index.get("Functional", 0)] if "Functional" in class_index else np.zeros(len(prepared))
    po = proba[:, class_index.get("ORF", 0)] if "ORF" in class_index else np.zeros(len(prepared))
    pp = proba[:, class_index.get("Pseudogene", 0)] if "Pseudogene" in class_index else np.zeros(len(prepared))
    predictions = []
    for idx, row in enumerate(prepared.itertuples(index=False)):
        if pf[idx] >= float(CONFIG.get("functional_threshold", 0.45)):
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
    return pd.DataFrame({
        "candidate_id": prepared["candidate_id"],
        "classification_pred": predictions,
        "prob_functional": pf,
        "prob_orf": po,
        "prob_pseudogene": pp,
    })

# SAFARI-Score: Phase 2 ML Classifier

Machine learning classifier for SAFARI-IGHJ candidate reclassification.

## Files

| File | Description |
|---|---|
| `safari_dl.py` | Core ML model (sklearn Random Forest + Bovini rules) |
| `phase1_winner_safari_score.py` | Phase 1 scoring module (RSS IC, FR4, ORF features) |
| `prepare_dl.py` | Evaluation harness (composite F1, LOOCV) |
| `predict_bos_taurus.py` | Bos taurus positive control predictions |
| `gold_standard_v2.tsv` | Ground truth labels (8 species from IMGT) |
| `candidates_all_with_embeddings.tsv` | 147 candidates with features (21 species) |
| `sklearn_predictions_all147.tsv` | ML predictions for all candidates |
| `loocv_report.tsv` | Leave-one-species-out cross-validation results |

## Usage

```python
import safari_dl

# Build model
gold_df = pd.read_csv("gold_standard_v2.tsv", sep="\t")
candidates_df = pd.read_csv("candidates_all_with_embeddings.tsv", sep="\t")
model = safari_dl.build_model(candidates_df, gold_df)

# Predict new candidates
predictions = safari_dl.predict(model, new_candidates_df)
```

## Performance

- Composite score: 0.881
- Gold standard F1: 1.0
- LOOCV mean accuracy: 0.898
- Bos taurus positive control: 100% concordance with Digger on shared predictions

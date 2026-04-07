"""
Phase A-2: ペプチド結合性分類モデル 学習スクリプト

data/peptide_dataset.csv を読み込み、RandomForest と LightGBM で二値分類モデルを学習する。
評価: AUC-ROC, Precision, Recall, F1
最良モデルを models/peptide_classifier.joblib に保存する。

使用方法:
    uv run python scripts/train_classifier.py
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from lightgbm import LGBMClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    ConfusionMatrixDisplay,
    auc,
    classification_report,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

DATASET_PATH = ROOT / "data" / "peptide_dataset.csv"
MODEL_DIR = ROOT / "models"
MODEL_PATH = MODEL_DIR / "peptide_classifier.joblib"

# 学習に使う特徴量カラム（sequenceとlabelを除く全カラム）
FEATURE_COLS = [
    "length",
    "net_charge",
    "avg_hydrophobicity",
    "gravy",
    "instability_index",
    "isoelectric_point",
    "aromaticity",
    "helix_fraction",
    "turn_fraction",
    "sheet_fraction",
    "molecular_weight",
]


# ─────────────────────────────────────────────
# データ読み込み
# ─────────────────────────────────────────────

def load_dataset() -> tuple[pd.DataFrame, np.ndarray]:
    """データセットを読み込み、特徴量DataFrame X とラベル y を返す。"""
    if not DATASET_PATH.exists():
        raise FileNotFoundError(
            f"{DATASET_PATH} が見つかりません。先に prepare_dataset.py を実行してください。"
        )
    df = pd.read_csv(DATASET_PATH)

    missing = [c for c in FEATURE_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"データセットに必要なカラムがありません: {missing}")

    X = df[FEATURE_COLS].astype(np.float32)
    y = df["label"].values.astype(np.int32)
    return X, y


# ─────────────────────────────────────────────
# モデル定義
# ─────────────────────────────────────────────

def build_models() -> dict[str, object]:
    """評価するモデルの辞書を返す。"""
    return {
        "RandomForest": Pipeline([
            ("scaler", StandardScaler()),
            ("clf", RandomForestClassifier(
                n_estimators=300,
                max_depth=None,
                min_samples_leaf=2,
                random_state=42,
                n_jobs=-1,
            )),
        ]),
        "LightGBM": Pipeline([
            ("scaler", StandardScaler()),
            ("clf", LGBMClassifier(
                n_estimators=300,
                learning_rate=0.05,
                num_leaves=31,
                random_state=42,
                n_jobs=-1,
                verbose=-1,
            )),
        ]),
    }


# ─────────────────────────────────────────────
# 評価
# ─────────────────────────────────────────────

def evaluate_models(
    models: dict[str, object],
    X: np.ndarray,
    y: np.ndarray,
) -> dict[str, float]:
    """
    5-fold 交差検証で各モデルの AUC-ROC を評価し、結果を表示する。

    Returns:
        モデル名 → 平均 AUC-ROC のdict
    """
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    results: dict[str, float] = {}

    print("\n" + "=" * 50)
    print("5-Fold Cross Validation (AUC-ROC)")
    print("=" * 50)

    for name, model in models.items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            scores = cross_val_score(model, X, y, cv=cv, scoring="roc_auc", n_jobs=-1)
        mean_auc = scores.mean()
        results[name] = mean_auc
        print(f"  {name:15s}: {mean_auc:.4f} ± {scores.std():.4f}")

    return results


def evaluate_best_model(model, X: np.ndarray, y: np.ndarray, name: str) -> None:
    """最良モデルをホールドアウトセットで詳細評価する。"""
    from sklearn.model_selection import train_test_split

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, stratify=y, random_state=42
    )
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    y_proba = model.predict_proba(X_test)[:, 1]

    print(f"\n{'=' * 50}")
    print(f"Holdout Evaluation: {name}")
    print("=" * 50)
    print(f"  AUC-ROC : {roc_auc_score(y_test, y_proba):.4f}")
    print(f"\n{classification_report(y_test, y_pred, target_names=['non-binder', 'binder'])}")

    # 特徴量重要度（RandomForest / LightGBM 共通）
    clf = model.named_steps["clf"]
    if hasattr(clf, "feature_importances_"):
        importances = clf.feature_importances_
        print("特徴量重要度 (上位5):")
        ranked = sorted(zip(FEATURE_COLS, importances), key=lambda x: -x[1])
        for feat, imp in ranked[:5]:
            print(f"  {feat:30s}: {imp:.4f}")


# ─────────────────────────────────────────────
# メイン
# ─────────────────────────────────────────────

def main():
    print("=== Phase A-2: モデル学習 ===")

    # データ読み込み
    print(f"\n[1/4] データセットを読み込み中: {DATASET_PATH}")
    X, y = load_dataset()
    print(f"  → サンプル数: {len(y)}  陽性: {y.sum()}  陰性: {(y == 0).sum()}")
    print(f"  → 特徴量: {FEATURE_COLS}")

    # モデル評価
    print("\n[2/4] モデルを交差検証で評価中...")
    models = build_models()
    results = evaluate_models(models, X, y)

    # 最良モデルを選択
    best_name = max(results, key=results.get)
    best_model = models[best_name]
    print(f"\n→ 最良モデル: {best_name} (AUC={results[best_name]:.4f})")

    # 詳細評価（ホールドアウト）
    print("\n[3/4] 最良モデルを詳細評価中...")
    evaluate_best_model(best_model, X, y, best_name)

    # 全データで再学習 → 保存
    print("\n[4/4] 全データで再学習してモデルを保存中...")
    best_model.fit(X, y)

    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {"model": best_model, "feature_cols": FEATURE_COLS, "model_name": best_name},
        MODEL_PATH,
    )
    print(f"  → 保存完了: {MODEL_PATH}")
    print(f"\n完了!")


if __name__ == "__main__":
    main()

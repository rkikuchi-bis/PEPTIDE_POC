"""
Phase A-2: ML スコアリングモジュール

学習済みモデル（models/peptide_classifier.joblib）を使って
ペプチド配列の結合確率スコアを返す。

Streamlit に依存しない純粋な Python モジュール。
"""
from __future__ import annotations

import warnings
from pathlib import Path

import joblib
import numpy as np

# モデルファイルのデフォルトパス
_DEFAULT_MODEL_PATH = Path(__file__).resolve().parent.parent / "models" / "peptide_classifier.joblib"

_STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# キャッシュ（アプリ起動中はモデルを一度だけロード）
_cache: dict = {}


def _load_model(model_path: Path = _DEFAULT_MODEL_PATH) -> dict:
    """モデルをキャッシュしてロードする。"""
    key = str(model_path)
    if key not in _cache:
        if not model_path.exists():
            raise FileNotFoundError(
                f"モデルファイルが見つかりません: {model_path}\n"
                "先に `uv run python scripts/train_classifier.py` を実行してください。"
            )
        _cache[key] = joblib.load(model_path)
    return _cache[key]


def _sequence_to_features(sequence: str) -> np.ndarray | None:
    """
    配列を特徴ベクトルに変換する。

    Returns:
        shape (1, n_features) の numpy 配列、または変換不可の場合 None
    """
    # 循環インポート回避のため遅延インポート
    from core.rescorer import _protparam_features

    clean = "".join(aa for aa in sequence.upper() if aa in _STANDARD_AA)
    if len(clean) < 3:
        return None

    feats = _protparam_features(clean)
    length = len(clean)
    net_charge = sum(1 for aa in clean if aa in "KRH") - sum(1 for aa in clean if aa in "DE")
    avg_hydrophobicity = sum(1 for aa in clean if aa in "AILMFWVY") / length

    # train_classifier.py の FEATURE_COLS と同じ順序
    feature_values = [
        length,
        net_charge,
        avg_hydrophobicity,
        feats["gravy"],
        feats["instability_index"],
        feats["isoelectric_point"],
        feats["aromaticity"],
        feats["helix_fraction"],
        feats["turn_fraction"],
        feats["sheet_fraction"],
        feats["molecular_weight"],
    ]
    return np.array([feature_values], dtype=np.float32)


def score_with_ml(sequence: str, model_path: Path = _DEFAULT_MODEL_PATH) -> float:
    """
    ペプチド配列の ML 結合確率スコアを返す。

    Args:
        sequence:   ペプチド配列（ワンレター）
        model_path: モデルファイルパス（省略時はデフォルト）

    Returns:
        0.0〜1.0 の結合確率スコア。モデル未ロード時は 0.5（中立）を返す。
    """
    try:
        artifact = _load_model(model_path)
        model = artifact["model"]
        feature_cols = artifact["feature_cols"]

        import pandas as pd
        X = _sequence_to_features(sequence)
        if X is None:
            return 0.5

        # DataFrame として渡すことで feature_names 警告を回避
        X_df = pd.DataFrame(X, columns=feature_cols)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            proba = model.predict_proba(X_df)[0][1]
        return float(round(proba, 4))

    except FileNotFoundError:
        # モデル未学習の場合は中立スコアを返してアプリを止めない
        return 0.5
    except Exception:
        return 0.5


def is_model_available(model_path: Path = _DEFAULT_MODEL_PATH) -> bool:
    """モデルファイルが存在するか確認する。"""
    return model_path.exists()

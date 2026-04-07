"""
Phase A-1: BioPython ProtParam による科学的特徴量ベースの再スコアリング

旧実装との対応：
  charge_match_score    → 等電点 pI を用いた電荷マッチ（Kyte-Doolittle準拠）
  hydrophobic_match_score → GRAVY（Grand Average of hYdropathicity）による疎水性マッチ
  complexity_score      → 芳香族性・二次構造傾向・配列多様性による複雑性
  stability_score       → 不安定性指数（Guruprasad et al., 1990）を追加

新規追加カラム：
  gravy, instability_index, isoelectric_point, aromaticity,
  helix_fraction, turn_fraction, sheet_fraction, stability_score
"""
from __future__ import annotations

import math
from typing import Dict

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from core.ml_scorer import is_model_available, score_with_ml

# ProteinAnalysis が受け付ける20種標準アミノ酸
_STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

_PHYSIOLOGICAL_PH = 7.4  # ヒト生理的pH基準


# ─────────────────────────────────────────────
# BioPython ProtParam 特徴量計算
# ─────────────────────────────────────────────

def _protparam_features(sequence: str) -> Dict[str, float]:
    """
    BioPython ProteinAnalysis から科学的特徴量を一括計算する。

    返り値キー：
        gravy              - Grand Average of hYdropathicity（Kyte-Doolittle）
        instability_index  - 不安定性指数（Guruprasad et al., 1990）; <40 = 安定
        isoelectric_point  - 等電点 pI
        aromaticity        - F/W/Y の割合
        helix_fraction     - αヘリックス傾向
        turn_fraction      - ターン傾向
        sheet_fraction     - βシート傾向
        molecular_weight   - 分子量 [Da]
    """
    # 非標準アミノ酸を除いてから渡す（安全のため）
    clean_seq = "".join(aa for aa in sequence.upper() if aa in _STANDARD_AA)
    if not clean_seq:
        return _fallback_features()

    try:
        pa = ProteinAnalysis(clean_seq)
        helix, turn, sheet = pa.secondary_structure_fraction()
        return {
            "gravy": pa.gravy(),
            "instability_index": pa.instability_index(),
            "isoelectric_point": pa.isoelectric_point(),
            "aromaticity": pa.aromaticity(),
            "helix_fraction": helix,
            "turn_fraction": turn,
            "sheet_fraction": sheet,
            "molecular_weight": pa.molecular_weight(),
        }
    except Exception:
        return _fallback_features()


def _fallback_features() -> Dict[str, float]:
    return {
        "gravy": 0.0,
        "instability_index": 50.0,
        "isoelectric_point": 7.0,
        "aromaticity": 0.0,
        "helix_fraction": 0.33,
        "turn_fraction": 0.33,
        "sheet_fraction": 0.33,
        "molecular_weight": 0.0,
    }


# ─────────────────────────────────────────────
# 個別スコア関数
# ─────────────────────────────────────────────

def _score_charge_match(isoelectric_point: float, pocket_charge: str) -> float:
    """
    等電点（pI）を用いたポケット電荷マッチスコア。
    生理的pH 7.4 基準で：
        pI > 7.4 → ペプチドは正電荷優勢（カチオン性）
        pI < 7.4 → ペプチドは負電荷優勢（アニオン性）
    """
    deviation = isoelectric_point - _PHYSIOLOGICAL_PH  # 正 = カチオン性

    if pocket_charge == "negative":
        # 負ポケット → カチオン性ペプチドが有利
        score = 0.5 + min(deviation * 0.055, 0.45)
    elif pocket_charge == "positive":
        # 正ポケット → アニオン性ペプチドが有利
        score = 0.5 - min(deviation * 0.055, 0.45)
    else:
        # 中性ポケット → pI が 7.4 付近（偏りが少ない）を好む
        score = 1.0 - min(abs(deviation) * 0.05, 0.45)

    return max(0.0, min(1.0, score))


def _score_hydrophobic_match(gravy: float, pocket_hydrophobicity: str) -> float:
    """
    GRAVY（Grand Average of hYdropathicity）を用いた疎水性マッチスコア。
    GRAVY > 0 = 疎水性, GRAVY < 0 = 親水性
    シグモイド変換で 0〜1 に正規化してからポケット傾向と照合。
    """
    # GRAVY を 0〜1 にシグモイド変換（GRAVY=0 → 0.5, +2 → 0.88, -2 → 0.12）
    normalized = 1.0 / (1.0 + math.exp(-gravy))

    if pocket_hydrophobicity == "high":
        score = 0.25 + normalized * 0.75
    elif pocket_hydrophobicity == "low":
        score = 0.25 + (1.0 - normalized) * 0.75
    else:
        # 中程度 → 極端でない（0.5 付近）を好む
        score = 1.0 - abs(normalized - 0.5) * 0.9

    return max(0.0, min(1.0, score))


def _score_stability(instability_index: float) -> float:
    """
    不安定性指数（Guruprasad et al., 1990）に基づく安定性スコア。
    < 40: 安定タンパク質（医薬用途に有利）
    > 40: 不安定な傾向
    """
    if instability_index <= 25:
        return 1.0
    elif instability_index <= 40:
        # 25〜40: なだらかに減点
        return 1.0 - (instability_index - 25) / 15 * 0.25
    elif instability_index <= 60:
        # 40〜60: さらに減点
        return 0.75 - (instability_index - 40) / 20 * 0.40
    else:
        return max(0.15, 0.35 - (instability_index - 60) / 100 * 0.20)


def _score_complexity(
    aromaticity: float,
    helix_fraction: float,
    sequence: str,
) -> float:
    """
    芳香族性・αヘリックス傾向・配列多様性による複雑性スコア。

    - 芳香族残基（F/W/Y）は π-π 相互作用・スタッキングで結合に寄与
    - αヘリックス傾向は多くのタンパク質ポケットへの結合で有利
    - 配列多様性（ユニーク残基割合）は繰り返し配列を減点
    """
    length = len(sequence)
    unique_frac = len(set(sequence)) / length if length > 0 else 0.0

    # 芳香族：0〜25% を [0, 1] に。多すぎると下がらない（上限で飽和）
    aromatic_score = min(aromaticity * 4.0, 1.0)

    # ヘリックス傾向：35% 付近を理想とする山形
    helix_score = max(0.0, 1.0 - abs(helix_fraction - 0.35) * 2.5)

    # 配列多様性：50% 以上のユニーク残基で満点
    diversity_score = min(unique_frac * 2.0, 1.0)

    score = 0.35 * aromatic_score + 0.35 * helix_score + 0.30 * diversity_score
    return max(0.0, min(1.0, score))


def _score_length_preference(
    length: int,
    preferred_min: int = 8,
    preferred_max: int = 12,
) -> float:
    """希望長さ範囲内なら満点。範囲外は距離に応じて減点。"""
    if preferred_min <= length <= preferred_max:
        return 1.0
    if length < preferred_min:
        return max(0.4, 1.0 - (preferred_min - length) * 0.12)
    return max(0.4, 1.0 - (length - preferred_max) * 0.10)


# ─────────────────────────────────────────────
# 注釈生成
# ─────────────────────────────────────────────

def _build_notes(
    pocket_charge: str,
    pocket_hydrophobicity: str,
    isoelectric_point: float,
    gravy: float,
    instability_index: float,
    aromaticity: float,
) -> str:
    notes = []

    # 電荷マッチ
    if pocket_charge == "negative" and isoelectric_point > 8.5:
        notes.append(f"cationic match (pI={isoelectric_point:.1f})")
    elif pocket_charge == "positive" and isoelectric_point < 6.0:
        notes.append(f"anionic match (pI={isoelectric_point:.1f})")
    else:
        notes.append(f"pI={isoelectric_point:.1f}")

    # 疎水性マッチ
    if pocket_hydrophobicity == "high" and gravy > 0.3:
        notes.append(f"hydrophobic match (GRAVY={gravy:.2f})")
    elif pocket_hydrophobicity == "low" and gravy < -0.3:
        notes.append(f"hydrophilic match (GRAVY={gravy:.2f})")
    else:
        notes.append(f"GRAVY={gravy:.2f}")

    # 安定性
    if instability_index < 40:
        notes.append(f"stable (II={instability_index:.0f})")
    else:
        notes.append(f"unstable (II={instability_index:.0f})")

    # 芳香族
    if aromaticity >= 0.15:
        notes.append(f"aromatic-rich ({aromaticity:.0%})")

    return "; ".join(notes)


# ─────────────────────────────────────────────
# メイン関数
# ─────────────────────────────────────────────

def rescore_candidates(
    df: pd.DataFrame,
    pocket_charge: str,
    pocket_hydrophobicity: str,
    preferred_len_min: int = 8,
    preferred_len_max: int = 12,
) -> pd.DataFrame:
    """
    Phase A-1: BioPython ProtParam 特徴量ベースの再スコアリング。

    スコア構成：
        charge_match_score      (0.30) - pI × ポケット電荷マッチ
        hydrophobic_match_score (0.25) - GRAVY × ポケット疎水性マッチ
        length_pref_score       (0.20) - 希望長さ範囲
        stability_score         (0.15) - 不安定性指数
        complexity_score        (0.10) - 芳香族性・ヘリックス・多様性
    """
    out = df.copy()

    rows_feats = []
    charge_scores, hydro_scores, length_scores = [], [], []
    stability_scores, complexity_scores, rescoring_scores = [], [], []
    notes_list = []

    for _, row in out.iterrows():
        seq = row["sequence"]
        length = int(row["length"])

        feats = _protparam_features(seq)
        rows_feats.append(feats)

        charge_score = _score_charge_match(feats["isoelectric_point"], pocket_charge)
        hydro_score = _score_hydrophobic_match(feats["gravy"], pocket_hydrophobicity)
        length_score = _score_length_preference(length, preferred_len_min, preferred_len_max)
        stability_score = _score_stability(feats["instability_index"])
        complexity_score = _score_complexity(
            feats["aromaticity"], feats["helix_fraction"], seq
        )

        rescoring_score = (
            0.30 * charge_score
            + 0.25 * hydro_score
            + 0.20 * length_score
            + 0.15 * stability_score
            + 0.10 * complexity_score
        )

        charge_scores.append(round(charge_score, 4))
        hydro_scores.append(round(hydro_score, 4))
        length_scores.append(round(length_score, 4))
        stability_scores.append(round(stability_score, 4))
        complexity_scores.append(round(complexity_score, 4))
        rescoring_scores.append(round(rescoring_score, 4))
        notes_list.append(
            _build_notes(
                pocket_charge,
                pocket_hydrophobicity,
                feats["isoelectric_point"],
                feats["gravy"],
                feats["instability_index"],
                feats["aromaticity"],
            )
        )

    # ProtParam 特徴量カラム
    out["gravy"] = [round(f["gravy"], 4) for f in rows_feats]
    out["instability_index"] = [round(f["instability_index"], 2) for f in rows_feats]
    out["isoelectric_point"] = [round(f["isoelectric_point"], 2) for f in rows_feats]
    out["aromaticity"] = [round(f["aromaticity"], 4) for f in rows_feats]
    out["helix_fraction"] = [round(f["helix_fraction"], 4) for f in rows_feats]
    out["turn_fraction"] = [round(f["turn_fraction"], 4) for f in rows_feats]
    out["sheet_fraction"] = [round(f["sheet_fraction"], 4) for f in rows_feats]
    out["molecular_weight"] = [round(f["molecular_weight"], 1) for f in rows_feats]

    # スコアカラム
    out["charge_match_score"] = charge_scores
    out["hydrophobic_match_score"] = hydro_scores
    out["length_pref_score"] = length_scores
    out["stability_score"] = stability_scores
    out["complexity_score"] = complexity_scores
    out["rescoring_score"] = rescoring_scores
    out["rescoring_notes"] = notes_list

    # ML スコア（モデル未学習の場合は 0.5 の中立値）
    ml_available = is_model_available()
    out["ml_score"] = [
        round(score_with_ml(seq), 4) for seq in out["sequence"]
    ]
    out["ml_score_available"] = ml_available

    # 最終統合スコア
    # ML モデルが利用可能な場合は ml_score を加味した重みに変更
    if ml_available:
        out["final_score"] = (
            0.20 * out["gen_score"]
            + 0.20 * out["property_score"]
            + 0.30 * out["rescoring_score"]
            + 0.30 * out["ml_score"]
        ).round(4)
    else:
        out["final_score"] = (
            0.30 * out["gen_score"]
            + 0.30 * out["property_score"]
            + 0.40 * out["rescoring_score"]
        ).round(4)

    return out

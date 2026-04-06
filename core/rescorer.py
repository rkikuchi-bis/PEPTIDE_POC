from __future__ import annotations

from typing import Dict
import pandas as pd

# 再スコア用の簡易アミノ酸カテゴリ
POSITIVE = set("KRH")
NEGATIVE = set("DE")
HYDROPHOBIC = set("AILMFWVY")
POLAR = set("STNQGP")

AROMATIC = set("FWY")
SMALL = set("AGST")
SPECIAL = set("CP")

def _calc_composition_features(sequence: str) -> Dict[str, float]:
    length = len(sequence)
    if length == 0:
        return {
            "frac_positive": 0.0,
            "frac_negative": 0.0,
            "frac_hydrophobic": 0.0,
            "frac_polar": 0.0,
            "frac_aromatic": 0.0,
            "frac_small": 0.0,
            "frac_special": 0.0,
        }

    return {
        "frac_positive": sum(1 for aa in sequence if aa in POSITIVE) / length,
        "frac_negative": sum(1 for aa in sequence if aa in NEGATIVE) / length,
        "frac_hydrophobic": sum(1 for aa in sequence if aa in HYDROPHOBIC) / length,
        "frac_polar": sum(1 for aa in sequence if aa in POLAR) / length,
        "frac_aromatic": sum(1 for aa in sequence if aa in AROMATIC) / length,
        "frac_small": sum(1 for aa in sequence if aa in SMALL) / length,
        "frac_special": sum(1 for aa in sequence if aa in SPECIAL) / length,
    }

def _score_charge_match(sequence: str, pocket_charge: str) -> float:
    feats = _calc_composition_features(sequence)
    pos = feats["frac_positive"]
    neg = feats["frac_negative"]

    if pocket_charge == "negative":
        # 負ポケットには正電荷がやや有利
        score = 0.4 + min(pos * 1.2, 0.5) - min(neg * 0.6, 0.2)
    elif pocket_charge == "positive":
        # 正ポケットには負電荷がやや有利
        score = 0.4 + min(neg * 1.2, 0.5) - min(pos * 0.6, 0.2)
    else:
        # 中性なら偏りすぎない方を少し好む
        imbalance = abs(pos - neg)
        score = 0.7 - min(imbalance * 0.8, 0.3)

    return max(min(score, 1.0), 0.0)

def _score_hydrophobic_match(sequence: str, pocket_hydrophobicity: str) -> float:
    feats = _calc_composition_features(sequence)
    hyd = feats["frac_hydrophobic"]
    polar = feats["frac_polar"]

    if pocket_hydrophobicity == "high":
        score = 0.4 + min(hyd * 1.1, 0.5) - min(polar * 0.3, 0.1)
    elif pocket_hydrophobicity == "low":
        score = 0.4 + min(polar * 0.8, 0.4) - min(hyd * 0.6, 0.2)
    else:
        # 中間は極端すぎない方がよい
        score = 0.8 - min(abs(hyd - 0.45) * 1.2, 0.35)

    return max(min(score, 1.0), 0.0)

def _score_length_preference(length: int, preferred_min: int = 8, preferred_max: int = 12) -> float:
    if preferred_min <= length <= preferred_max:
        return 1.0
    if length < preferred_min:
        return max(0.4, 1.0 - (preferred_min - length) * 0.12)
    return max(0.4, 1.0 - (length - preferred_max) * 0.10)

def _score_sequence_complexity(sequence: str) -> float:
    """
    極端な単調性を軽く減点。
    """
    length = len(sequence)
    if length == 0:
        return 0.0

    unique_frac = len(set(sequence)) / length
    aromatic_frac = sum(1 for aa in sequence if aa in AROMATIC) / length
    special_frac = sum(1 for aa in sequence if aa in SPECIAL) / length

    score = 0.55
    score += min(unique_frac * 0.35, 0.25)
    score += min(aromatic_frac * 0.25, 0.10)
    score -= min(special_frac * 0.30, 0.10)

    return max(min(score, 1.0), 0.0)

def _build_notes(sequence: str, pocket_charge: str, pocket_hydrophobicity: str) -> str:
    feats = _calc_composition_features(sequence)
    notes = []

    if pocket_charge == "negative" and feats["frac_positive"] >= 0.20:
        notes.append("good positive-charge complementarity")
    elif pocket_charge == "positive" and feats["frac_negative"] >= 0.15:
        notes.append("good negative-charge complementarity")
    elif pocket_charge == "neutral":
        notes.append("balanced charge profile")

    if pocket_hydrophobicity == "high" and feats["frac_hydrophobic"] >= 0.35:
        notes.append("good hydrophobic match")
    elif pocket_hydrophobicity == "low" and feats["frac_polar"] >= 0.30:
        notes.append("good polar match")
    elif pocket_hydrophobicity == "medium":
        notes.append("moderate hydrophobic balance")

    if len(set(sequence)) / len(sequence) >= 0.6:
        notes.append("reasonable sequence diversity")

    return "; ".join(notes) if notes else "no strong signal"

def rescore_candidates(
    df: pd.DataFrame,
    pocket_charge: str,
    pocket_hydrophobicity: str,
    preferred_len_min: int = 8,
    preferred_len_max: int = 12,
) -> pd.DataFrame:
    """
    軽量な再スコア。
    後で docking / Boltz-2 / classifier に置き換える想定。
    """
    out = df.copy()

    charge_scores = []
    hydro_scores = []
    length_scores = []
    complexity_scores = []
    rescoring_scores = []
    notes_list = []

    for _, row in out.iterrows():
        seq = row["sequence"]
        length = int(row["length"])

        charge_score = _score_charge_match(seq, pocket_charge)
        hydro_score = _score_hydrophobic_match(seq, pocket_hydrophobicity)
        length_score = _score_length_preference(length, preferred_len_min, preferred_len_max)
        complexity_score = _score_sequence_complexity(seq)

        rescoring_score = (
            0.35 * charge_score +
            0.30 * hydro_score +
            0.20 * length_score +
            0.15 * complexity_score
        )

        charge_scores.append(round(charge_score, 4))
        hydro_scores.append(round(hydro_score, 4))
        length_scores.append(round(length_score, 4))
        complexity_scores.append(round(complexity_score, 4))
        rescoring_scores.append(round(rescoring_score, 4))
        notes_list.append(_build_notes(seq, pocket_charge, pocket_hydrophobicity))

    out["charge_match_score"] = charge_scores
    out["hydrophobic_match_score"] = hydro_scores
    out["length_pref_score"] = length_scores
    out["complexity_score"] = complexity_scores
    out["rescoring_score"] = rescoring_scores
    out["rescoring_notes"] = notes_list

    # 最終統合スコア
    out["final_score"] = (
        0.30 * out["gen_score"] +
        0.30 * out["property_score"] +
        0.40 * out["rescoring_score"]
    ).round(4)

    return out
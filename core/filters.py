from __future__ import annotations

from typing import Dict

import pandas as pd

HYDRO_SCALE: Dict[str, float] = {
    "A": 0.70, "C": 0.78, "D": 0.11, "E": 0.11, "F": 0.81,
    "G": 0.46, "H": 0.14, "I": 1.00, "K": 0.07, "L": 0.92,
    "M": 0.74, "N": 0.11, "P": 0.32, "Q": 0.11, "R": 0.00,
    "S": 0.41, "T": 0.42, "V": 0.96, "W": 0.40, "Y": 0.36,
}

CHARGE_MAP: Dict[str, int] = {
    "K": 1, "R": 1, "H": 1,
    "D": -1, "E": -1,
}

def calc_length(sequence: str) -> int:
    return len(sequence)

def calc_net_charge(sequence: str) -> int:
    return sum(CHARGE_MAP.get(res, 0) for res in sequence)

def calc_avg_hydrophobicity(sequence: str) -> float:
    if not sequence:
        return 0.0
    return sum(HYDRO_SCALE.get(res, 0.0) for res in sequence) / len(sequence)

def has_excessive_repeat(sequence: str, max_repeat_residue: int = 2) -> bool:
    if not sequence:
        return False

    current = 1
    for i in range(1, len(sequence)):
        if sequence[i] == sequence[i - 1]:
            current += 1
            if current > max_repeat_residue:
                return True
        else:
            current = 1
    return False

def _canonical_signature(sequence: str) -> str:
    return sequence

def add_basic_properties(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["length"] = out["sequence"].apply(calc_length)
    out["net_charge"] = out["sequence"].apply(calc_net_charge)
    out["avg_hydrophobicity"] = out["sequence"].apply(calc_avg_hydrophobicity)
    out["duplicate_key"] = out["sequence"].apply(_canonical_signature)
    return out

def _compute_property_score(row: pd.Series) -> float:
    score = 1.0

    if row["length"] < 6 or row["length"] > 16:
        score -= 0.2

    if abs(row["net_charge"]) > 4:
        score -= 0.2

    if row["avg_hydrophobicity"] > 0.70:
        score -= 0.2

    if not row["repeat_ok"]:
        score -= 0.3

    score += min(float(row.get("gen_score", 0.0)) * 0.2, 0.2)

    return max(min(score, 1.0), 0.0)

def apply_filters(
    df: pd.DataFrame,
    min_len: int,
    max_len: int,
    max_abs_charge: int,
    max_hydrophobicity: float,
    max_repeat_residue: int = 2,
    remove_near_duplicates: bool = True,
) -> pd.DataFrame:
    out = df.copy()

    out["repeat_ok"] = ~out["sequence"].apply(
        lambda x: has_excessive_repeat(x, max_repeat_residue=max_repeat_residue)
    )

    out = out[(out["length"] >= min_len) & (out["length"] <= max_len)]
    out = out[out["net_charge"].abs() <= max_abs_charge]
    out = out[out["avg_hydrophobicity"] <= max_hydrophobicity]
    out = out[out["repeat_ok"]]

    if remove_near_duplicates:
        out["duplicate_flag"] = out.duplicated(subset=["duplicate_key"], keep="first")
        out = out[~out["duplicate_flag"]]
    else:
        out["duplicate_flag"] = False

    out["property_score"] = out.apply(_compute_property_score, axis=1)

    return out.reset_index(drop=True)
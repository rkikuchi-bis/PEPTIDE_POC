from __future__ import annotations

from typing import List
import pandas as pd


def hamming_distance(seq1: str, seq2: str) -> int:
    """
    同じ長さの配列同士の Hamming 距離。
    """
    if len(seq1) != len(seq2):
        raise ValueError("Hamming distance requires sequences of equal length.")
    return sum(a != b for a, b in zip(seq1, seq2))


def normalized_hamming_distance(seq1: str, seq2: str) -> float:
    """
    Hamming 距離を 0-1 に正規化。
    長さが異なる場合は別物として扱うため 1.0 を返す。
    """
    if len(seq1) != len(seq2):
        return 1.0
    if len(seq1) == 0:
        return 0.0
    return hamming_distance(seq1, seq2) / len(seq1)


def _compute_min_distance_to_kept(
    sequence: str,
    kept_sequences: List[str],
) -> float:
    """
    既採用配列群に対する最小 normalized Hamming distance を返す。
    比較対象がない場合は 1.0。
    """
    comparable_distances = []

    for kept_seq in kept_sequences:
        if len(sequence) != len(kept_seq):
            continue
        dist = normalized_hamming_distance(sequence, kept_seq)
        comparable_distances.append(dist)

    if not comparable_distances:
        return 1.0

    return min(comparable_distances)


def diversify_candidates(
    df: pd.DataFrame,
    sequence_col: str = "sequence",
    score_col: str = "final_score",
    min_normalized_distance: float = 0.30,
    max_candidates: int | None = None,
) -> pd.DataFrame:
    """
    高スコア順に候補を見て、既採用候補と近すぎる配列を間引く。

    ルール:
    - score_col の高い順に走査
    - 同じ長さの配列について normalized Hamming distance を計算
    - 最小距離が min_normalized_distance 未満なら除外
    - 長さが違う配列は別物として扱う
    """
    if df.empty:
        out = df.copy()
        out["diversity_kept"] = pd.Series(dtype=bool)
        out["diversity_min_distance"] = pd.Series(dtype=float)
        return out

    work = df.sort_values(by=score_col, ascending=False).reset_index(drop=True).copy()

    kept_rows = []
    kept_sequences: List[str] = []
    kept_min_distances: List[float] = []

    for _, row in work.iterrows():
        seq = row[sequence_col]
        min_dist = _compute_min_distance_to_kept(seq, kept_sequences)

        keep = min_dist >= min_normalized_distance

        if keep:
            row_copy = row.copy()
            row_copy["diversity_kept"] = True
            row_copy["diversity_min_distance"] = round(min_dist, 4)

            kept_rows.append(row_copy)
            kept_sequences.append(seq)
            kept_min_distances.append(min_dist)

            if max_candidates is not None and len(kept_rows) >= max_candidates:
                break

    if not kept_rows:
        out = work.iloc[0:0].copy()
        out["diversity_kept"] = pd.Series(dtype=bool)
        out["diversity_min_distance"] = pd.Series(dtype=float)
        return out

    out = pd.DataFrame(kept_rows).reset_index(drop=True)
    return out
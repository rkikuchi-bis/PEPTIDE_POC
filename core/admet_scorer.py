"""
admet_scorer.py
ペプチド生物活性ヒューリスティックスコアリングモジュール。

Peptide Ranker (UCD Dublin) は 2025年時点でサービス停止中のため、
既存の ProtParam 特徴量（gravy, instability_index, aromaticity, length,
net_charge）を組み合わせたローカルヒューリスティックで近似する。

スコアは 0〜1（高いほど生物活性ペプチドらしい）。
参考: Mooney et al. (2012) Peptide Ranker; Wang et al. (2022) bioactive peptide properties.

計算根拠:
  - 疎水性 (GRAVY): 中程度（-1.5〜0.5）が最も生物活性ペプチドに多い
  - 安定性 (instability_index): < 40 = 安定; 不安定ペプチドは in vivo で機能しにくい
  - 長さ: 6〜15残基が生物活性ペプチドの大多数を占める
  - 芳香族性 (aromaticity): F/W/Y は受容体との π-stacking に寄与
  - 電荷 (net_charge): 絶対値が小さいほど膜透過性・全般的な結合能が高い
"""

import math
from typing import Optional


def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def _bioactivity_heuristic(row) -> float:
    """
    1行分の ProtParam 特徴量から生物活性スコア (0〜1) を計算する。
    """

    # 1. 疎水性 (GRAVY): ピーク -0.3、幅 σ=0.9 のガウス分布
    gravy = float(row.get("gravy") or 0.0)
    gravy_score = math.exp(-((gravy - (-0.3)) ** 2) / (2 * 0.9 ** 2))

    # 2. 安定性 (instability_index): 40 を境にシグモイド減衰
    instability = float(row.get("instability_index") or 50.0)
    stability_score = _sigmoid(-0.15 * (instability - 35.0))

    # 3. 長さ: ピーク 10残基、幅 σ=4 のガウス分布
    length = int(row.get("length") or 10)
    length_score = math.exp(-((length - 10) ** 2) / (2 * 4.0 ** 2))

    # 4. 芳香族性 (aromaticity): 0〜0.15 の範囲を線形に評価
    aromaticity = float(row.get("aromaticity") or 0.0)
    aromaticity_score = min(1.0, aromaticity / 0.15)

    # 5. 電荷 (net_charge): 絶対値が大きいほど減点
    net_charge = abs(int(row.get("net_charge") or 0))
    charge_score = _sigmoid(-0.5 * (net_charge - 2.0))

    # 重み付け平均
    score = (
        0.30 * gravy_score
        + 0.30 * stability_score
        + 0.20 * length_score
        + 0.10 * aromaticity_score
        + 0.10 * charge_score
    )
    return round(min(1.0, max(0.0, score)), 4)


def score_top_candidates(df, top_n: int = 20):
    """
    DataFrame の上位 top_n 候補に生物活性ヒューリスティックスコアを付与する。

    Returns:
        bioactivity_score カラムを追加した DataFrame（コピー）。
        ProtParam 特徴量がない行は NaN。
    """
    import pandas as pd

    df = df.copy()
    df["bioactivity_score"] = float("nan")

    required = {"gravy", "instability_index", "length", "aromaticity", "net_charge"}
    if not required.issubset(df.columns):
        return df  # ProtParam 特徴量がなければスキップ

    target_idx = df.index[:top_n].tolist()
    for idx in target_idx:
        row = df.loc[idx]
        df.at[idx, "bioactivity_score"] = _bioactivity_heuristic(row)

    return df

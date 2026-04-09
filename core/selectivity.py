"""
Direction B: Selectivity（選択性スコアリング）— Phase C-1

LigandForge との差別化:
  LigandForge は「1ターゲットに結合するペプチドを生成」するだけ。
  このモジュールは「ターゲットAには結合するがオフターゲットBには結合しにくい」
  ペプチドを定量的に評価する。

設計方針:
  - ターゲットの rescoring_score（charge match + hydrophobicity match + length + stability）
    はポケット特性に依存するため、ターゲット固有の指標として使える
  - ML / MPNN スコアは配列内因性なので両ターゲット共通 → 差し引きで消える
  - selectivity_score = rescoring_score_target - rescoring_score_offtarget
  - 正値ほどターゲット選択的、負値ほどオフターゲット側に有利
"""
from __future__ import annotations

import pandas as pd

from core.rescorer import rescore_candidates


def compute_selectivity(
    result_df: pd.DataFrame,
    pocket_charge_offtarget: str,
    pocket_hydrophobicity_offtarget: str,
    offtarget_label: str = "Off-target",
    preferred_len_min: int = 8,
    preferred_len_max: int = 12,
    structure_text_offtarget: str | None = None,
    file_format_offtarget: str | None = None,
    pocket_centroid_offtarget: tuple[float, float, float] | None = None,
) -> pd.DataFrame:
    """
    diversity 後の候補に対してオフターゲット向けの rescoring を実行し、
    選択性スコアを計算して result_df に追加する。

    追加カラム:
        offtarget_rescoring_score  - オフターゲット向け rescoring score
        selectivity_score          - rescoring_target - rescoring_offtarget
                                     正値 = ターゲット選択的
        offtarget_label            - オフターゲットのラベル文字列

    注:
        ML / ProteinMPNN スコアは配列内因性のため共通。
        ポケット特性に依存する rescoring_score の差で選択性を評価する。
    """
    if result_df.empty:
        return result_df

    # rescore_candidates に必要な最小カラムを渡す
    # （gen_score / property_score は rescoring_score の計算には不要だが
    #    final_score 計算で参照されるため渡す）
    required_cols = ["sequence", "length", "net_charge", "avg_hydrophobicity",
                     "gen_score", "property_score"]
    available_cols = [c for c in required_cols if c in result_df.columns]
    input_df = result_df[available_cols].copy()

    # オフターゲット向けに再スコアリング（ProteinMPNN / ESMFold は省略: 時間節約）
    # PEPFOLD_MAX_SEQS=0 で ESMFold を無効化する
    import os
    orig = os.environ.get("PEPFOLD_MAX_SEQS")
    os.environ["PEPFOLD_MAX_SEQS"] = "0"
    try:
        offtarget_df = rescore_candidates(
            input_df,
            pocket_charge=pocket_charge_offtarget,
            pocket_hydrophobicity=pocket_hydrophobicity_offtarget,
            preferred_len_min=preferred_len_min,
            preferred_len_max=preferred_len_max,
            structure_text=structure_text_offtarget,
            file_format=file_format_offtarget,
            pocket_centroid=pocket_centroid_offtarget,
        )
    finally:
        if orig is None:
            os.environ.pop("PEPFOLD_MAX_SEQS", None)
        else:
            os.environ["PEPFOLD_MAX_SEQS"] = orig

    out = result_df.copy()
    out["offtarget_rescoring_score"] = offtarget_df["rescoring_score"].values
    out["selectivity_score"] = (
        out["rescoring_score"] - out["offtarget_rescoring_score"]
    ).round(4)
    out["offtarget_label"] = offtarget_label

    return out

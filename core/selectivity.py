"""
Direction B: Selectivity（選択性スコアリング）— Phase C-1 / C-2

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


def compute_docking_selectivity(
    result_df: pd.DataFrame,
    receptor_pdbqt_offtarget: str,
    box_center_offtarget: tuple[float, float, float],
    box_size_offtarget: tuple[float, float, float] = (20.0, 20.0, 20.0),
    top_n: int = 15,
    exhaustiveness: int = 8,
) -> pd.DataFrame:
    """
    Phase C-2: ドッキングスコア差分による選択性。

    ターゲットへの docking_score（Phase B-1 で計算済み）と
    オフターゲットへの docking_score の差から選択性を定量化する。

    Vina スコアは負値が強い結合を示す（単位: kcal/mol）。
    docking_selectivity_score = docking_score_offtarget - docking_score_target
        正値 = ターゲット選択的（ターゲットへの結合がより強い）
        負値 = オフターゲット選択的（避けるべき候補）

    前提: result_df に docking_score カラムが存在すること（Phase B-1 実行後）。

    追加カラム:
        docking_score_offtarget      - オフターゲット Vina スコア [kcal/mol]
        docking_selectivity_score    - docking_score_offtarget - docking_score_target
    """
    from core.docking import dock_peptide, is_vina_available, FLEXIBLE_DOCKING_MAX_LENGTH

    if result_df.empty:
        return result_df

    if "docking_score" not in result_df.columns:
        return result_df

    if not is_vina_available():
        return result_df

    out = result_df.copy()
    out["docking_score_offtarget"] = float("nan")

    for idx in out.index:
        target_score = out.at[idx, "docking_score"]
        if pd.isna(target_score):
            continue  # ターゲットドッキング未実施の行はスキップ

        seq = out.at[idx, "sequence"]
        use_flexible = len(seq) <= FLEXIBLE_DOCKING_MAX_LENGTH

        score_ot = dock_peptide(
            sequence=seq,
            receptor_pdbqt_path=receptor_pdbqt_offtarget,
            box_center=box_center_offtarget,
            box_size=box_size_offtarget,
            exhaustiveness=exhaustiveness,
            flexible=use_flexible,
        )
        if score_ot is not None:
            out.at[idx, "docking_score_offtarget"] = score_ot

    # 差分計算: offtarget - target（正値 = ターゲット選択的）
    # 剛体ドッキング失敗（正値スコア）を含む行は信頼性がないため NaN にする
    has_ot = out["docking_score_offtarget"].notna()
    valid = (
        has_ot
        & (out["docking_score"] <= 0)
        & (out["docking_score_offtarget"] <= 0)
    )
    out["docking_selectivity_score"] = float("nan")
    out.loc[valid, "docking_selectivity_score"] = (
        out.loc[valid, "docking_score_offtarget"]
        - out.loc[valid, "docking_score"]
    ).round(3)

    return out

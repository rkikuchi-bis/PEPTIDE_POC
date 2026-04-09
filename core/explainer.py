"""
Direction A: Explainability（説明可能性）

各ペプチド候補について、スコア根拠を日本語の自然言語で説明する。
LigandForge 等のブラックボックスモデルとの差別化ポイント。

使い方:
    from core.explainer import explain_candidate
    explanation = explain_candidate(row, pocket_charge, pocket_hydrophobicity)
"""
from __future__ import annotations

import pandas as pd


def explain_candidate(
    row: pd.Series,
    pocket_charge: str,
    pocket_hydrophobicity: str,
) -> str:
    """
    ペプチド候補の推薦理由を日本語で自然言語生成する。

    Args:
        row: result_df の1行（Series）
        pocket_charge: "negative" / "positive" / "neutral"
        pocket_hydrophobicity: "high" / "medium" / "low"

    Returns:
        日本語説明文（複数文を「。」で連結）
    """
    parts: list[str] = []

    # ── 1. 総合評価 ───────────────────────────────────────────
    score = float(row.get("final_score", 0))
    seq = str(row.get("sequence", ""))
    length = int(row.get("length", len(seq)))

    if score >= 0.75:
        parts.append(f"配列 {seq}（{length}残基）は複数の評価軸で高い適合性を示す有力候補です")
    elif score >= 0.60:
        parts.append(f"配列 {seq}（{length}残基）はバランスよく高評価を示す候補です")
    else:
        parts.append(f"配列 {seq}（{length}残基）は一部の特性に優れた候補です")

    # ── 2. 電荷適合 ───────────────────────────────────────────
    pI = float(row.get("isoelectric_point", 7.0))
    charge_score = float(row.get("charge_match_score", 0.5))
    if charge_score >= 0.75:
        if pocket_charge == "negative":
            parts.append(
                f"等電点 pI={pI:.1f} で生理的 pH 7.4 では正電荷優勢（カチオン性）となり、"
                f"負電荷ポケットへの静電引力が期待できます"
            )
        elif pocket_charge == "positive":
            parts.append(
                f"等電点 pI={pI:.1f} で生理的 pH 7.4 では負電荷優勢（アニオン性）となり、"
                f"正電荷ポケットへの静電引力が期待できます"
            )
        else:
            parts.append(
                f"等電点 pI={pI:.1f} で中性付近の電荷バランスを持ち、中性ポケットに適合します"
            )
    elif charge_score < 0.40:
        parts.append(
            f"等電点 pI={pI:.1f} はポケット電荷（{pocket_charge}）との適合が低く、"
            f"静電的な結合力は期待しにくい状況です"
        )
    else:
        parts.append(f"等電点 pI={pI:.1f}（電荷マッチスコア: {charge_score:.2f}）")

    # ── 3. 疎水性適合 ──────────────────────────────────────────
    gravy = float(row.get("gravy", 0.0))
    hydro_score = float(row.get("hydrophobic_match_score", 0.5))
    gravy_label = "疎水性" if gravy > 0 else "親水性"
    if hydro_score >= 0.70:
        if pocket_hydrophobicity == "high" and gravy > 0.2:
            parts.append(
                f"GRAVY={gravy:.2f} の疎水性配列で、疎水性ポケット内部への埋没と"
                f"疎水効果による安定的結合が期待できます"
            )
        elif pocket_hydrophobicity == "low" and gravy < -0.2:
            parts.append(
                f"GRAVY={gravy:.2f} の親水性配列で、親水性ポケットへの水素結合・"
                f"静電相互作用が主な結合力になります"
            )
        else:
            parts.append(
                f"GRAVY={gravy:.2f} でポケット疎水性（{pocket_hydrophobicity}）と"
                f"良好にマッチしています"
            )
    elif hydro_score < 0.40:
        parts.append(
            f"GRAVY={gravy:.2f}（{gravy_label}）はポケット疎水性（{pocket_hydrophobicity}）と"
            f"ミスマッチです"
        )
    else:
        parts.append(f"GRAVY={gravy:.2f}（疎水性マッチスコア: {hydro_score:.2f}）")

    # ── 4. 安定性 ─────────────────────────────────────────────
    ii = float(row.get("instability_index", 50.0))
    if ii < 40:
        parts.append(
            f"不安定性指数 II={ii:.0f}（基準値 <40）で安定と予測され、"
            f"薬剤候補として代謝安定性が期待できます"
        )
    elif ii > 60:
        parts.append(
            f"不安定性指数 II={ii:.0f}（>60）で不安定な傾向があり、"
            f"生体内での分解リスクに注意が必要です"
        )

    # ── 5. 芳香族・二次構造 ────────────────────────────────────
    arom = float(row.get("aromaticity", 0.0))
    helix = float(row.get("helix_fraction", 0.0))
    if arom >= 0.15:
        aromatic_aas = [aa for aa in seq if aa in "FWY"]
        aa_str = "".join(aromatic_aas) if aromatic_aas else ""
        parts.append(
            f"芳香族残基（{aa_str}）を {arom:.0%} 含み、π-π スタッキングや"
            f"疎水性相互作用による結合への寄与が期待されます"
        )
    if helix >= 0.45:
        parts.append(
            f"αヘリックス傾向 {helix:.0%} で安定したヘリックス構造を形成しやすく、"
            f"多くのタンパク質ポケットへの結合に有利な立体配置が期待できます"
        )

    # ── 6. ML 結合確率 ─────────────────────────────────────────
    ml = row.get("ml_score")
    ml_available = bool(row.get("ml_score_available", False))
    if ml_available and ml is not None and not pd.isna(ml):
        ml = float(ml)
        if ml >= 0.70:
            parts.append(
                f"機械学習モデル（LightGBM、PDB 陽性例学習済み）の結合確率 {ml:.1%} で、"
                f"既知結合ペプチドの特徴パターンに近い配列です"
            )
        elif ml < 0.35:
            parts.append(
                f"機械学習モデルの結合確率 {ml:.1%} で、既知結合パターンからの乖離が見られます"
            )

    # ── 7. ProteinMPNN スコア ──────────────────────────────────
    mpnn = row.get("proteinmpnn_score")
    mpnn_available = bool(row.get("proteinmpnn_score_available", False))
    if mpnn_available and mpnn is not None and not pd.isna(mpnn):
        mpnn = float(mpnn)
        receptor_cond = bool(row.get("proteinmpnn_receptor_conditioned", False))
        esmfold_used = bool(row.get("esmfold_backbone_used", False))
        mode_label = (
            "ESMFold 骨格・受容体条件付き" if receptor_cond and esmfold_used
            else "受容体条件付き" if receptor_cond
            else "構造フリー"
        )
        if mpnn >= 0.70:
            parts.append(
                f"ProteinMPNN 設計可能性スコア {mpnn:.3f}（{mode_label}）で"
                f"高い設計適合性が評価されています"
            )
        elif mpnn < 0.35:
            parts.append(
                f"ProteinMPNN スコア {mpnn:.3f}（{mode_label}）で設計可能性は低めです"
            )

    # ── 8. ドッキングスコア ────────────────────────────────────
    dock = row.get("docking_score")
    if dock is not None and not pd.isna(dock):
        dock = float(dock)
        docking_mode = str(row.get("docking_mode", "rigid"))
        if dock < -7.0:
            parts.append(
                f"AutoDock Vina スコア {dock:.1f} kcal/mol（{docking_mode}ドッキング）で"
                f"強い結合親和性が示唆されます"
            )
        elif dock < -5.0:
            parts.append(
                f"AutoDock Vina スコア {dock:.1f} kcal/mol（{docking_mode}ドッキング）で"
                f"中程度の結合親和性が示唆されます"
            )
        elif dock > 0:
            parts.append(
                f"Vina スコア {dock:.1f} kcal/mol（剛体ドッキング）は正値で、"
                f"絶対値より相対ランキングでの参照を推奨します"
            )

    # ── 9. 選択性スコア ────────────────────────────────────────
    sel = row.get("selectivity_score")
    if sel is not None and not pd.isna(sel):
        sel = float(sel)
        offtarget_label = str(row.get("offtarget_label", "オフターゲット"))
        if sel >= 0.15:
            parts.append(
                f"選択性スコア Δ={sel:.3f} で {offtarget_label} よりもターゲットへの"
                f"選択的結合が期待できます（差別化の重要な指標）"
            )
        elif sel <= -0.10:
            parts.append(
                f"選択性スコア Δ={sel:.3f} で {offtarget_label} への非選択的結合が懸念されます"
            )
        else:
            parts.append(
                f"選択性スコア Δ={sel:.3f}（ターゲットとオフターゲットの差が小さい）"
            )

    return "。".join(parts) + "。"

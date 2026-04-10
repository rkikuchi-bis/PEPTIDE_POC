"""
Direction A: Explainability（説明可能性）

各ペプチド候補について、スコア根拠を日英両言語の自然言語で説明する。
LigandForge 等のブラックボックスモデルとの差別化ポイント。

使い方:
    from core.explainer import explain_candidate
    ja, en = explain_candidate(row, pocket_charge, pocket_hydrophobicity)
"""
from __future__ import annotations

import pandas as pd


def explain_candidate(
    row: pd.Series,
    pocket_charge: str,
    pocket_hydrophobicity: str,
) -> tuple[str, str]:
    """
    ペプチド候補の推薦理由を日英両言語で自然言語生成する。

    Returns:
        (ja_text, en_text) のタプル
    """
    ja_parts: list[str] = []
    en_parts: list[str] = []

    # ── 1. 総合評価 ───────────────────────────────────────────
    score = float(row.get("final_score", 0))
    seq = str(row.get("sequence", ""))
    length = int(row.get("length", len(seq)))

    if score >= 0.75:
        ja_parts.append(f"配列 {seq}（{length}残基）は複数の評価軸で高い適合性を示す有力候補です")
        en_parts.append(f"Sequence {seq} ({length} residues) is a strong candidate with high fitness across multiple scoring axes")
    elif score >= 0.60:
        ja_parts.append(f"配列 {seq}（{length}残基）はバランスよく高評価を示す候補です")
        en_parts.append(f"Sequence {seq} ({length} residues) shows well-balanced scores across evaluations")
    else:
        ja_parts.append(f"配列 {seq}（{length}残基）は一部の特性に優れた候補です")
        en_parts.append(f"Sequence {seq} ({length} residues) excels in selected properties")

    # ── 2. 電荷適合 ───────────────────────────────────────────
    pI = float(row.get("isoelectric_point", 7.0))
    charge_score = float(row.get("charge_match_score", 0.5))
    if charge_score >= 0.75:
        if pocket_charge == "negative":
            ja_parts.append(
                f"等電点 pI={pI:.1f} で生理的 pH 7.4 では正電荷優勢（カチオン性）となり、"
                f"負電荷ポケットへの静電引力が期待できます"
            )
            en_parts.append(
                f"pI={pI:.1f} yields a cationic peptide at physiological pH 7.4, "
                f"enabling electrostatic attraction to the negatively charged pocket"
            )
        elif pocket_charge == "positive":
            ja_parts.append(
                f"等電点 pI={pI:.1f} で生理的 pH 7.4 では負電荷優勢（アニオン性）となり、"
                f"正電荷ポケットへの静電引力が期待できます"
            )
            en_parts.append(
                f"pI={pI:.1f} renders the peptide anionic at pH 7.4, "
                f"promoting electrostatic binding to the positively charged pocket"
            )
        else:
            ja_parts.append(
                f"等電点 pI={pI:.1f} で中性付近の電荷バランスを持ち、中性ポケットに適合します"
            )
            en_parts.append(
                f"pI={pI:.1f} provides near-neutral charge balance, well-suited to a neutral pocket"
            )
    elif charge_score < 0.40:
        ja_parts.append(
            f"等電点 pI={pI:.1f} はポケット電荷（{pocket_charge}）との適合が低く、"
            f"静電的な結合力は期待しにくい状況です"
        )
        en_parts.append(
            f"pI={pI:.1f} shows poor electrostatic compatibility with the {pocket_charge} pocket; "
            f"charge-driven binding is unlikely"
        )
    else:
        ja_parts.append(f"等電点 pI={pI:.1f}（電荷マッチスコア: {charge_score:.2f}）")
        en_parts.append(f"pI={pI:.1f} (charge match score: {charge_score:.2f})")

    # ── 3. 疎水性適合 ──────────────────────────────────────────
    gravy = float(row.get("gravy", 0.0))
    hydro_score = float(row.get("hydrophobic_match_score", 0.5))
    gravy_label = "疎水性" if gravy > 0 else "親水性"
    gravy_label_en = "hydrophobic" if gravy > 0 else "hydrophilic"
    if hydro_score >= 0.70:
        if pocket_hydrophobicity == "high" and gravy > 0.2:
            ja_parts.append(
                f"GRAVY={gravy:.2f} の疎水性配列で、疎水性ポケット内部への埋没と"
                f"疎水効果による安定的結合が期待できます"
            )
            en_parts.append(
                f"GRAVY={gravy:.2f} indicates a hydrophobic sequence expected to bury into the "
                f"hydrophobic pocket and form stable hydrophobic interactions"
            )
        elif pocket_hydrophobicity == "low" and gravy < -0.2:
            ja_parts.append(
                f"GRAVY={gravy:.2f} の親水性配列で、親水性ポケットへの水素結合・"
                f"静電相互作用が主な結合力になります"
            )
            en_parts.append(
                f"GRAVY={gravy:.2f} indicates a hydrophilic sequence; hydrogen bonds and "
                f"electrostatic interactions will dominate binding to the polar pocket"
            )
        else:
            ja_parts.append(
                f"GRAVY={gravy:.2f} でポケット疎水性（{pocket_hydrophobicity}）と"
                f"良好にマッチしています"
            )
            en_parts.append(
                f"GRAVY={gravy:.2f} is well-matched to the pocket hydrophobicity ({pocket_hydrophobicity})"
            )
    elif hydro_score < 0.40:
        ja_parts.append(
            f"GRAVY={gravy:.2f}（{gravy_label}）はポケット疎水性（{pocket_hydrophobicity}）と"
            f"ミスマッチです"
        )
        en_parts.append(
            f"GRAVY={gravy:.2f} ({gravy_label_en}) is mismatched with the {pocket_hydrophobicity} pocket hydrophobicity"
        )
    else:
        ja_parts.append(f"GRAVY={gravy:.2f}（疎水性マッチスコア: {hydro_score:.2f}）")
        en_parts.append(f"GRAVY={gravy:.2f} (hydrophobic match score: {hydro_score:.2f})")

    # ── 4. 安定性 ─────────────────────────────────────────────
    ii = float(row.get("instability_index", 50.0))
    if ii < 40:
        ja_parts.append(
            f"不安定性指数 II={ii:.0f}（基準値 <40）で安定と予測され、"
            f"薬剤候補として代謝安定性が期待できます"
        )
        en_parts.append(
            f"Instability index II={ii:.0f} (threshold <40) predicts a stable peptide "
            f"with good metabolic stability as a drug candidate"
        )
    elif ii > 60:
        ja_parts.append(
            f"不安定性指数 II={ii:.0f}（>60）で不安定な傾向があり、"
            f"生体内での分解リスクに注意が必要です"
        )
        en_parts.append(
            f"Instability index II={ii:.0f} (>60) suggests a tendency toward instability; "
            f"in vivo degradation risk should be considered"
        )

    # ── 5. 芳香族・二次構造 ────────────────────────────────────
    arom = float(row.get("aromaticity", 0.0))
    helix = float(row.get("helix_fraction", 0.0))
    if arom >= 0.15:
        aromatic_aas = [aa for aa in seq if aa in "FWY"]
        aa_str = "".join(aromatic_aas) if aromatic_aas else ""
        ja_parts.append(
            f"芳香族残基（{aa_str}）を {arom:.0%} 含み、π-π スタッキングや"
            f"疎水性相互作用による結合への寄与が期待されます"
        )
        en_parts.append(
            f"Contains {arom:.0%} aromatic residues ({aa_str}), suggesting contributions from "
            f"π-π stacking and hydrophobic interactions"
        )
    if helix >= 0.45:
        ja_parts.append(
            f"αヘリックス傾向 {helix:.0%} で安定したヘリックス構造を形成しやすく、"
            f"多くのタンパク質ポケットへの結合に有利な立体配置が期待できます"
        )
        en_parts.append(
            f"Helix propensity of {helix:.0%} favors stable α-helical conformation, "
            f"an advantageous geometry for many protein pocket interactions"
        )

    # ── 6. ML 結合確率 ─────────────────────────────────────────
    ml = row.get("ml_score")
    ml_available = bool(row.get("ml_score_available", False))
    if ml_available and ml is not None and not pd.isna(ml):
        ml = float(ml)
        if ml >= 0.70:
            ja_parts.append(
                f"機械学習モデル（LightGBM、PDB 陽性例学習済み）の結合確率 {ml:.1%} で、"
                f"既知結合ペプチドの特徴パターンに近い配列です"
            )
            en_parts.append(
                f"ML model (LightGBM trained on PDB binding data) predicts {ml:.1%} binding probability, "
                f"closely matching patterns of known binding peptides"
            )
        elif ml < 0.35:
            ja_parts.append(
                f"機械学習モデルの結合確率 {ml:.1%} で、既知結合パターンからの乖離が見られます"
            )
            en_parts.append(
                f"ML binding probability of {ml:.1%} deviates from known binding patterns"
            )

    # ── 7. ProteinMPNN スコア ──────────────────────────────────
    mpnn = row.get("proteinmpnn_score")
    mpnn_available = bool(row.get("proteinmpnn_score_available", False))
    if mpnn_available and mpnn is not None and not pd.isna(mpnn):
        mpnn = float(mpnn)
        receptor_cond = bool(row.get("proteinmpnn_receptor_conditioned", False))
        esmfold_used = bool(row.get("esmfold_backbone_used", False))
        mode_label_ja = (
            "ESMFold 骨格・受容体条件付き" if receptor_cond and esmfold_used
            else "受容体条件付き" if receptor_cond
            else "構造フリー"
        )
        mode_label_en = (
            "ESMFold backbone + receptor-conditioned" if receptor_cond and esmfold_used
            else "receptor-conditioned" if receptor_cond
            else "structure-free"
        )
        if mpnn >= 0.70:
            ja_parts.append(
                f"ProteinMPNN 設計可能性スコア {mpnn:.3f}（{mode_label_ja}）で"
                f"高い設計適合性が評価されています"
            )
            en_parts.append(
                f"ProteinMPNN designability score {mpnn:.3f} ({mode_label_en}) "
                f"indicates high structural design compatibility"
            )
        elif mpnn < 0.35:
            ja_parts.append(
                f"ProteinMPNN スコア {mpnn:.3f}（{mode_label_ja}）で設計可能性は低めです"
            )
            en_parts.append(
                f"ProteinMPNN score {mpnn:.3f} ({mode_label_en}) indicates limited designability"
            )

    # ── 8. ドッキングスコア ────────────────────────────────────
    dock = row.get("docking_score")
    if dock is not None and not pd.isna(dock):
        dock = float(dock)
        docking_mode = str(row.get("docking_mode", "rigid"))
        docking_mode_en = docking_mode
        if dock < -7.0:
            ja_parts.append(
                f"AutoDock Vina スコア {dock:.1f} kcal/mol（{docking_mode}ドッキング）で"
                f"強い結合親和性が示唆されます"
            )
            en_parts.append(
                f"AutoDock Vina score {dock:.1f} kcal/mol ({docking_mode_en} docking) "
                f"suggests strong binding affinity"
            )
        elif dock < -5.0:
            ja_parts.append(
                f"AutoDock Vina スコア {dock:.1f} kcal/mol（{docking_mode}ドッキング）で"
                f"中程度の結合親和性が示唆されます"
            )
            en_parts.append(
                f"AutoDock Vina score {dock:.1f} kcal/mol ({docking_mode_en} docking) "
                f"indicates moderate binding affinity"
            )
        elif dock > 0:
            ja_parts.append(
                f"Vina スコア {dock:.1f} kcal/mol（剛体ドッキング）は正値で、"
                f"絶対値より相対ランキングでの参照を推奨します"
            )
            en_parts.append(
                f"Vina score {dock:.1f} kcal/mol (rigid docking) is positive; "
                f"use for relative ranking rather than absolute affinity"
            )

    # ── 9. 選択性スコア ────────────────────────────────────────
    sel = row.get("selectivity_score")
    if sel is not None and not pd.isna(sel):
        sel = float(sel)
        offtarget_label = str(row.get("offtarget_label", "オフターゲット"))
        if sel >= 0.15:
            ja_parts.append(
                f"選択性スコア Δ={sel:.3f} で {offtarget_label} よりもターゲットへの"
                f"選択的結合が期待できます（差別化の重要な指標）"
            )
            en_parts.append(
                f"Selectivity score Δ={sel:.3f} indicates preferential binding to the target "
                f"over {offtarget_label} — a key differentiator"
            )
        elif sel <= -0.10:
            ja_parts.append(
                f"選択性スコア Δ={sel:.3f} で {offtarget_label} への非選択的結合が懸念されます"
            )
            en_parts.append(
                f"Selectivity score Δ={sel:.3f} raises concern for non-selective binding to {offtarget_label}"
            )
        else:
            ja_parts.append(
                f"選択性スコア Δ={sel:.3f}（ターゲットとオフターゲットの差が小さい）"
            )
            en_parts.append(
                f"Selectivity score Δ={sel:.3f} (small difference between target and off-target)"
            )

    ja_text = "。".join(ja_parts) + "。"
    en_text = ". ".join(en_parts) + "."
    return ja_text, en_text

import pandas as pd
import streamlit as st


def render_results(result_df, pdb_summary=None):
    if result_df is not None:
        if pdb_summary is not None:
            st.subheader("Selected structure summary")
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Chain", pdb_summary["selected_chain"] if pdb_summary.get("selected_chain") else "-")
            c2.metric("Mode", pdb_summary.get("source_mode", "-"))
            c3.metric("Auto charge", pdb_summary["pocket_charge_guess"])
            c4.metric("Auto hydrophobicity", pdb_summary["pocket_hydrophobicity_guess"])

        st.subheader("Results")

        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Remaining candidates", len(result_df))
        col2.metric("Top final score", f"{result_df['final_score'].max():.3f}" if len(result_df) > 0 else "0.000")
        col3.metric("Median length", int(result_df["length"].median()) if len(result_df) > 0 else 0)
        if "motif_compare_score" in result_df.columns and len(result_df) > 0 and result_df["motif_compare_score"].notna().any():
            col4.metric("Top motif compare", f"{result_df['motif_compare_score'].fillna(0).max():.3f}")
        else:
            col4.metric("Top motif compare", "-")

        display_cols = [
            # ── 順位・配列 ──
            "rank",
            "sequence",
            # ── 統合スコア ──
            "final_score",
            "ml_score",
            "rescoring_score",
            "gen_score",
            "property_score",
            # ── 基本物性 ──
            "length",
            "net_charge",
            "avg_hydrophobicity",
            # ── Phase A-1: ProtParam 科学的特徴量 ──
            "isoelectric_point",
            "gravy",
            "instability_index",
            "aromaticity",
            "molecular_weight",
            "helix_fraction",
            "turn_fraction",
            "sheet_fraction",
            # ── スコア内訳 ──
            "charge_match_score",
            "hydrophobic_match_score",
            "stability_score",
            "complexity_score",
            # ── 多様性 ──
            "diversity_kept",
            "diversity_min_distance",
            # ── モチーフ比較 ──
            "best_known_motif",
            "known_sequence_identity_max",
            "known_kmer_jaccard_max",
            "motif_compare_score",
        ]
        display_cols = [c for c in display_cols if c in result_df.columns]

        st.dataframe(
            result_df[display_cols],
            use_container_width=True,
            height=500,
        )

        st.subheader("Top candidate details")
        if len(result_df) > 0:
            selected_rank = st.selectbox("Select rank", result_df["rank"].tolist(), index=0)
            row = result_df[result_df["rank"] == selected_rank].iloc[0]

            st.markdown(f"**Sequence:** `{row['sequence']}`")
            st.write(f"- Final score: {row['final_score']:.3f}")
            if "ml_score" in row.index and pd.notna(row["ml_score"]):
                ml_label = "ML binding probability (Phase A-2)"
                st.write(f"- {ml_label}: {row['ml_score']:.3f}")
            st.write(f"- Rescoring score: {row['rescoring_score']:.3f}")
            st.write(f"- Generated score: {row['gen_score']:.3f}")
            st.write(f"- Property score: {row['property_score']:.3f}")
            st.write(f"- Length: {int(row['length'])}")
            st.write(f"- Net charge: {int(row['net_charge'])}")
            st.write(f"- Avg hydrophobicity: {row['avg_hydrophobicity']:.3f}")
            st.write(f"- Charge match score: {row['charge_match_score']:.3f}")
            st.write(f"- Hydrophobic match score: {row['hydrophobic_match_score']:.3f}")

            if "length_pref_score" in row.index and pd.notna(row["length_pref_score"]):
                st.write(f"- Length preference score: {row['length_pref_score']:.3f}")
            if "stability_score" in row.index and pd.notna(row["stability_score"]):
                st.write(f"- Stability score: {row['stability_score']:.3f}")
            if "complexity_score" in row.index and pd.notna(row["complexity_score"]):
                st.write(f"- Complexity score: {row['complexity_score']:.3f}")

            # ── Phase A-1: ProtParam 特徴量 ──
            protparam_fields = [
                ("isoelectric_point", "Isoelectric point (pI)", ".2f"),
                ("gravy", "GRAVY (hydropathicity)", ".3f"),
                ("instability_index", "Instability index", ".1f"),
                ("aromaticity", "Aromaticity", ".3f"),
                ("molecular_weight", "Molecular weight [Da]", ".1f"),
                ("helix_fraction", "Helix fraction", ".3f"),
                ("turn_fraction", "Turn fraction", ".3f"),
                ("sheet_fraction", "Sheet fraction", ".3f"),
            ]
            has_protparam = any(k in row.index and pd.notna(row[k]) for k, *_ in protparam_fields)
            if has_protparam:
                st.markdown("**ProtParam features (Phase A-1)**")
                for key, label, fmt in protparam_fields:
                    if key in row.index and pd.notna(row[key]):
                        st.write(f"- {label}: {row[key]:{fmt}}")

            st.write(f"- Diversity kept: {row['diversity_kept']}")
            if "diversity_min_distance" in row.index and pd.notna(row["diversity_min_distance"]):
                st.write(f"- Diversity min distance: {row['diversity_min_distance']:.3f}")

            if "aromatic_ratio" in row.index:
                st.write(
                    f"- Aromatic ratio: {float(row['aromatic_ratio']):.3f}"
                    if pd.notna(row["aromatic_ratio"])
                    else "- Aromatic ratio: -"
                )
            if "charge_pattern" in row.index:
                st.write(f"- Charge pattern: {row['charge_pattern']}")
            if "best_known_motif" in row.index:
                st.write(f"- Best known motif: {row['best_known_motif']}")
            if "known_aromatic_gap_min" in row.index and pd.notna(row["known_aromatic_gap_min"]):
                st.write(f"- Known aromatic gap min: {row['known_aromatic_gap_min']:.3f}")
            if "known_charge_pattern_similarity_max" in row.index and pd.notna(row["known_charge_pattern_similarity_max"]):
                st.write(f"- Known charge pattern similarity max: {row['known_charge_pattern_similarity_max']:.3f}")
            if "known_sequence_identity_max" in row.index and pd.notna(row["known_sequence_identity_max"]):
                st.write(f"- Known sequence identity max: {row['known_sequence_identity_max']:.3f}")
            if "known_kmer_jaccard_max" in row.index and pd.notna(row["known_kmer_jaccard_max"]):
                st.write(f"- Known k-mer jaccard max: {row['known_kmer_jaccard_max']:.3f}")
            if "motif_compare_score" in row.index and pd.notna(row["motif_compare_score"]):
                st.write(f"- Motif compare score: {row['motif_compare_score']:.3f}")
            if "rescoring_notes" in row.index:
                st.write(f"- Notes: {row['rescoring_notes']}")
    else:
        st.info("左の設定を調整して、Generate and Filter を押してください。")

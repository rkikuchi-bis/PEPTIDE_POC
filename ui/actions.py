"""
Post-pipeline action panels.

Each function handles one optional processing phase.
Conditions are checked internally; app.py calls these unconditionally.
"""
from __future__ import annotations

import os
import tempfile

import streamlit as st


def render_pipeline_run(params: dict) -> None:
    """Main generation pipeline + automatic Phase C-1 selectivity. Triggers on Run button."""
    if not params.get("run_button"):
        return

    with st.spinner("Generating candidates..."):
        from core.pipeline import run_pipeline
        from core.pdb_utils import (
            detect_structure_format,
            get_pocket_ca_centroid,
            parse_structure_text,
        )

        structure_bytes = params.get("structure_bytes")
        structure_filename = params.get("structure_filename", "")
        pdb_summary = params.get("pdb_summary")
        structure_text = file_format = pocket_centroid = None

        if structure_bytes is not None and pdb_summary is not None:
            try:
                file_format = detect_structure_format(structure_filename)
                structure_text = structure_bytes.decode("utf-8", errors="ignore")
                structure_obj = parse_structure_text(structure_text, file_format=file_format)
                pocket_centroid = get_pocket_ca_centroid(structure_obj, pdb_summary)
            except Exception:
                structure_text = file_format = pocket_centroid = None

        result_df, saved_path = run_pipeline(
            num_candidates=params["num_candidates"],
            min_len=params["min_len"],
            max_len=params["max_len"],
            pocket_charge=params["pocket_charge"],
            pocket_hydrophobicity=params["pocket_hydrophobicity"],
            avoid_cysteine=params["avoid_cysteine"],
            max_abs_charge=params["max_abs_charge"],
            max_hydrophobicity=params["max_hydrophobicity"],
            max_repeat_residue=params["max_repeat_residue"],
            preferred_len_min=params["preferred_len_min"],
            preferred_len_max=params["preferred_len_max"],
            use_diversity_filter=params["use_diversity_filter"],
            min_diversity_distance=params["min_diversity_distance"],
            max_diverse_candidates=params["max_diverse_candidates"],
            known_sequences=params["known_sequences"],
            target_name=params["target_name"],
            structure_text=structure_text,
            file_format=file_format,
            pocket_centroid=pocket_centroid,
            seed_sequence=params.get("seed_sequence"),
            variant_strategies=params.get("variant_strategies"),
        )

        st.session_state.update({
            "result_df": result_df,
            "docking_done": False,
            "docking_selectivity_done": False,
            "admet_done": False,
            "viewer_selected_rank": 1,
            "pocket_centroid": pocket_centroid,
            "last_run_structure_filename": params.get("structure_filename"),
            "pocket_charge": params["pocket_charge"],
            "pocket_hydrophobicity": params["pocket_hydrophobicity"],
        })

        n = len(result_df) if result_df is not None else 0
        mpnn_mode = (
            "receptor-conditioned ProteinMPNN"
            if pocket_centroid is not None
            else "ProteinMPNN (structure-free)"
        )
        seed = params.get("seed_sequence")
        gen_mode = f"variant generation (seed: {seed})" if seed else "random generation"
        st.success(
            f"Done ✅  {n} candidates generated. "
            f"Generation: {gen_mode}.  Scoring: {mpnn_mode}.  CSV: `{saved_path}`"
        )

    # Phase C-1: Selectivity (runs automatically when configured)
    sel_params = params.get("selectivity_params")
    if sel_params is not None and st.session_state["result_df"] is not None:
        _run_selectivity(params, sel_params)

    st.rerun()


def _run_selectivity(params: dict, sel_params: dict) -> None:
    """Compute Phase C-1 rescoring-based selectivity scores."""
    from core.pdb_utils import detect_structure_format, get_pocket_ca_centroid, parse_structure_text
    from core.selectivity import compute_selectivity

    with st.spinner(f"Computing selectivity vs. {sel_params['offtarget_label']}..."):
        structure_text_ot = file_format_ot = pocket_centroid_ot = None

        ot_bytes = sel_params.get("offtarget_structure_bytes")
        ot_filename = sel_params.get("offtarget_structure_filename")
        pdb_summary_ot = sel_params.get("pdb_summary_offtarget")

        if ot_bytes is not None and pdb_summary_ot is not None:
            try:
                file_format_ot = detect_structure_format(ot_filename)
                structure_text_ot = ot_bytes.decode("utf-8", errors="ignore")
                ot_obj = parse_structure_text(structure_text_ot, file_format=file_format_ot)
                pocket_centroid_ot = get_pocket_ca_centroid(ot_obj, pdb_summary_ot)
            except Exception:
                structure_text_ot = file_format_ot = pocket_centroid_ot = None

        result_df = compute_selectivity(
            st.session_state["result_df"],
            pocket_charge_offtarget=sel_params["pocket_charge_offtarget"],
            pocket_hydrophobicity_offtarget=sel_params["pocket_hydrophobicity_offtarget"],
            offtarget_label=sel_params["offtarget_label"],
            preferred_len_min=params["preferred_len_min"],
            preferred_len_max=params["preferred_len_max"],
            structure_text_offtarget=structure_text_ot,
            file_format_offtarget=file_format_ot,
            pocket_centroid_offtarget=pocket_centroid_ot,
        )

        lam = sel_params.get("selectivity_lambda", 0.0)
        result_df["selective_final_score"] = (
            result_df["final_score"] + lam * result_df["selectivity_score"]
        ).round(4)

        if lam > 0.0:
            result_df = result_df.sort_values(
                "selective_final_score", ascending=False
            ).reset_index(drop=True)
            result_df["rank"] = range(1, len(result_df) + 1)

        st.session_state["result_df"] = result_df

        top_sel = result_df["selectivity_score"].max()
        lam_msg = (
            f" (λ={lam:.2f} applied to ranking)"
            if lam > 0.0
            else " (λ=0: not reflected in ranking)"
        )
        st.info(
            f"Selectivity computed vs. {sel_params['offtarget_label']}. "
            f"Top: {top_sel:.3f}{lam_msg}"
        )


def render_docking_button(params: dict) -> None:
    """Phase B-1: Show AutoDock Vina button and run docking on click."""
    if not (
        params.get("enable_docking")
        and params.get("docking_params") is not None
        and st.session_state["result_df"] is not None
        and not st.session_state["docking_done"]
    ):
        return

    structure_bytes = params.get("structure_bytes")
    structure_filename = params.get("structure_filename", "receptor.pdb")

    if structure_bytes is None:
        st.info("Docking requires a structure file. Please upload one.")
        return

    if not st.button("Run Docking on top candidates", type="secondary"):
        return

    from core.docking import dock_top_candidates, prepare_receptor_pdbqt

    dp = params["docking_params"]
    suffix = "." + structure_filename.rsplit(".", 1)[-1]

    with st.spinner(
        f"Preparing receptor and docking top {dp['top_n']} candidates "
        f"(exhaustiveness={dp['exhaustiveness']})..."
    ):
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
            tmp.write(structure_bytes)
            tmp_path = tmp.name

        try:
            receptor_pdbqt = prepare_receptor_pdbqt(tmp_path)
            if receptor_pdbqt is None:
                st.error("Receptor PDBQT preparation failed. Is obabel installed?")
                return
            result_df = dock_top_candidates(
                st.session_state["result_df"],
                receptor_pdbqt_path=receptor_pdbqt,
                box_center=dp["box_center"],
                box_size=dp["box_size"],
                top_n=dp["top_n"],
                exhaustiveness=dp["exhaustiveness"],
            )
            st.session_state["result_df"] = result_df
            st.session_state["docking_done"] = True
            docked = result_df["docking_score"].notna().sum()
            st.success(f"Docking complete. {docked} candidates scored.")
        finally:
            os.unlink(tmp_path)


def render_docking_selectivity_button(params: dict) -> None:
    """Phase C-2: Show docking-based selectivity button and run on click."""
    sel_params = params.get("selectivity_params")
    has_ot_structure = (
        sel_params is not None
        and sel_params.get("offtarget_structure_bytes") is not None
        and sel_params.get("pdb_summary_offtarget") is not None
    )
    if not (
        st.session_state["docking_done"]
        and has_ot_structure
        and st.session_state["result_df"] is not None
        and not st.session_state["docking_selectivity_done"]
    ):
        return

    from core.docking import is_vina_available

    if not is_vina_available():
        return

    if not st.button("🎯 Run Docking Selectivity (Phase C-2)", type="secondary"):
        return

    from core.docking import compute_pocket_center, prepare_receptor_pdbqt
    from core.selectivity import compute_docking_selectivity

    dp = params.get("docking_params") or {}
    ot_bytes = sel_params["offtarget_structure_bytes"]
    ot_filename = sel_params["offtarget_structure_filename"]
    pdb_summary_ot = sel_params["pdb_summary_offtarget"]
    offtarget_label = sel_params.get("offtarget_label", "Off-target")
    suffix = "." + ot_filename.rsplit(".", 1)[-1]

    with st.spinner(f"Preparing off-target receptor and docking vs. {offtarget_label}..."):
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
            tmp.write(ot_bytes)
            tmp_path = tmp.name

        try:
            receptor_pdbqt_ot = prepare_receptor_pdbqt(tmp_path)
            if receptor_pdbqt_ot is None:
                st.error("Off-target receptor PDBQT preparation failed.")
                return
            box_center_ot = compute_pocket_center(ot_bytes, ot_filename, pdb_summary_ot)
            if box_center_ot is None:
                st.error("Could not compute off-target pocket center coordinates.")
                return
            result_df = compute_docking_selectivity(
                st.session_state["result_df"],
                receptor_pdbqt_offtarget=receptor_pdbqt_ot,
                box_center_offtarget=box_center_ot,
                box_size_offtarget=dp.get("box_size", (20.0, 20.0, 20.0)),
                exhaustiveness=dp.get("exhaustiveness", 8),
            )
            st.session_state["result_df"] = result_df
            st.session_state["docking_selectivity_done"] = True
            n_scored = result_df["docking_selectivity_score"].notna().sum()
            top_sel = result_df["docking_selectivity_score"].max()
            st.success(
                f"Phase C-2 complete ✅  {n_scored} candidates scored. "
                f"Top docking selectivity: {top_sel:.2f} kcal/mol (vs {offtarget_label})"
            )
        finally:
            os.unlink(tmp_path)


def render_bioactivity_button() -> None:
    """Heuristic bioactivity scoring for top candidates."""
    if st.session_state["result_df"] is None or st.session_state["admet_done"]:
        return

    top_n = min(20, len(st.session_state["result_df"]))
    if not st.button(f"🧪 Score Bioactivity (top {top_n})", type="secondary"):
        return

    from core.admet_scorer import score_top_candidates

    with st.spinner("Scoring bioactivity..."):
        result_df = score_top_candidates(st.session_state["result_df"], top_n=top_n)
        st.session_state["result_df"] = result_df
        st.session_state["admet_done"] = True
        n_scored = result_df["bioactivity_score"].notna().sum()
        st.success(f"Bioactivity scoring complete ✅  {n_scored} candidates scored.")

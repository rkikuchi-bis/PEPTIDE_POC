import os
import tempfile

import streamlit as st

from core.pipeline import run_pipeline
from core.pdb_utils import get_pocket_ca_centroid
from ui.sidebar import render_sidebar
from ui.results import render_results
from ui.structure_viewer import render_viewer_section


# =========================
# App config
# =========================
st.set_page_config(page_title="Peptide Discovery PoC", layout="wide")

st.title("Peptide Discovery PoC")
st.caption(
    "Input → Candidate Generation → Filter → Rescoring → Diversity Control → Motif Comparison → Output  "
    "｜  入力 → 候補生成 → フィルタ → 再スコア → 多様性制御 → 既知モチーフ比較 → 表示"
)


# =========================
# Session state init
# =========================
defaults = {
    "result_df": None,
    "docking_done": False,
    "rcsb_results": [],
    "rcsb_selected_index": 0,
    "downloaded_structure_bytes": None,
    "downloaded_structure_name": None,
    "downloaded_structure_pdb_id": None,
    "rcsb_last_query": "",
    "structure_source": "Upload local file",
    "pocket_charge": "negative",
    "pocket_hydrophobicity": "medium",
}
for key, val in defaults.items():
    if key not in st.session_state:
        st.session_state[key] = val


# =========================
# Sidebar
# =========================
params = render_sidebar()


# =========================
# 3D Structure Viewer
# =========================
render_viewer_section(
    structure_bytes=params["structure_bytes"],
    structure_filename=params["structure_filename"],
    pdb_summary=params["pdb_summary"],
)


# =========================
# Run generation pipeline
# =========================
if params["run_button"]:
    with st.spinner("Generating candidates..."):
        # 受容体構造を条件付けとした ProteinMPNN スコアリング用に情報を準備
        structure_bytes = params.get("structure_bytes")
        structure_filename = params.get("structure_filename", "")
        pdb_summary = params.get("pdb_summary")
        structure_text = None
        file_format = None
        pocket_centroid = None

        if structure_bytes is not None and pdb_summary is not None:
            from core.pdb_utils import detect_structure_format, parse_structure_text
            try:
                file_format = detect_structure_format(structure_filename)
                structure_text = structure_bytes.decode("utf-8", errors="ignore")
                structure_obj = parse_structure_text(structure_text, file_format=file_format)
                pocket_centroid = get_pocket_ca_centroid(structure_obj, pdb_summary)
            except Exception:
                structure_text = None
                file_format = None
                pocket_centroid = None

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
        )
        st.session_state["result_df"] = result_df
        st.session_state["docking_done"] = False
        # ポケット特性を session_state に保存（説明生成に使用）
        st.session_state["pocket_charge"] = params["pocket_charge"]
        st.session_state["pocket_hydrophobicity"] = params["pocket_hydrophobicity"]

        # 受容体条件付きスコアリングの使用状況を表示
        if pocket_centroid is not None:
            cx, cy, cz = pocket_centroid
            st.success(
                f"Done. Receptor-conditioned ProteinMPNN enabled "
                f"(centroid: {cx:.1f}, {cy:.1f}, {cz:.1f}). "
                f"Saved CSV to: {saved_path}"
            )
        else:
            st.success(f"Done. Saved CSV to: {saved_path}")

    # =========================
    # Selectivity (Phase C-1)
    # =========================
    selectivity_params = params.get("selectivity_params")
    if selectivity_params is not None and st.session_state["result_df"] is not None:
        with st.spinner(
            f"Computing selectivity vs. {selectivity_params['offtarget_label']}..."
        ):
            from core.selectivity import compute_selectivity
            from core.pdb_utils import detect_structure_format, parse_structure_text

            # オフターゲット構造が指定されていれば pocket centroid を計算
            structure_text_offtarget = None
            file_format_offtarget = None
            pocket_centroid_offtarget = None
            ot_bytes = selectivity_params.get("offtarget_structure_bytes")
            ot_filename = selectivity_params.get("offtarget_structure_filename")
            pdb_summary_ot = selectivity_params.get("pdb_summary_offtarget")
            if ot_bytes is not None and pdb_summary_ot is not None:
                try:
                    file_format_offtarget = detect_structure_format(ot_filename)
                    structure_text_offtarget = ot_bytes.decode("utf-8", errors="ignore")
                    ot_structure_obj = parse_structure_text(structure_text_offtarget, file_format=file_format_offtarget)
                    pocket_centroid_offtarget = get_pocket_ca_centroid(ot_structure_obj, pdb_summary_ot)
                except Exception:
                    structure_text_offtarget = None
                    file_format_offtarget = None
                    pocket_centroid_offtarget = None

            result_df = compute_selectivity(
                st.session_state["result_df"],
                pocket_charge_offtarget=selectivity_params["pocket_charge_offtarget"],
                pocket_hydrophobicity_offtarget=selectivity_params["pocket_hydrophobicity_offtarget"],
                offtarget_label=selectivity_params["offtarget_label"],
                preferred_len_min=params["preferred_len_min"],
                preferred_len_max=params["preferred_len_max"],
                structure_text_offtarget=structure_text_offtarget,
                file_format_offtarget=file_format_offtarget,
                pocket_centroid_offtarget=pocket_centroid_offtarget,
            )

            # λ > 0 のとき、選択性をランキングに反映する
            lam = selectivity_params.get("selectivity_lambda", 0.0)
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
            lam_msg = f"、λ={lam:.2f} でランキングに反映済み" if lam > 0.0 else "（λ=0: ランキングへの反映なし）"
            st.info(
                f"Selectivity computed vs. {selectivity_params['offtarget_label']}. "
                f"Top selectivity score: {top_sel:.3f}{lam_msg}"
            )


# =========================
# Docking (Phase B-1)
# =========================
if (
    params.get("enable_docking")
    and params.get("docking_params") is not None
    and st.session_state["result_df"] is not None
    and not st.session_state["docking_done"]
):
    dp = params["docking_params"]
    structure_bytes = params.get("structure_bytes")
    structure_filename = params.get("structure_filename", "receptor.pdb")

    if structure_bytes is not None:
        if st.button("Run Docking on top candidates", type="secondary"):
            from core.docking import dock_top_candidates, prepare_receptor_pdbqt

            with st.spinner(
                f"Preparing receptor and docking top {dp['top_n']} candidates "
                f"(exhaustiveness={dp['exhaustiveness']})..."
            ):
                # 構造バイトを一時ファイルに書き出し
                suffix = "." + structure_filename.rsplit(".", 1)[-1]
                with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
                    tmp.write(structure_bytes)
                    tmp_pdb_path = tmp.name

                try:
                    receptor_pdbqt = prepare_receptor_pdbqt(tmp_pdb_path)
                    if receptor_pdbqt is None:
                        st.error("受容体 PDBQT の準備に失敗しました。")
                    else:
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
                    os.unlink(tmp_pdb_path)
    else:
        st.info("ドッキングには構造ファイルのアップロードが必要です。")


# =========================
# Result display
# =========================
render_results(
    st.session_state["result_df"],
    params["pdb_summary"],
    pocket_charge=st.session_state["pocket_charge"],
    pocket_hydrophobicity=st.session_state["pocket_hydrophobicity"],
)

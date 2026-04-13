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
st.set_page_config(page_title="Peptide Discovery PoC", layout="wide", page_icon="🧬")

st.title("🧬 Peptide Discovery PoC")
st.caption(
    "タンパク質ポケット構造を入力として、結合候補ペプチドを生成・スコアリング・可視化する創薬支援ツール。"
    "  |  Generates, scores, and visualizes peptide candidates from a protein pocket structure."
)


# =========================
# Session state init
# =========================
defaults = {
    "result_df": None,
    "docking_done": False,
    "admet_done": False,
    "viewer_selected_rank": 1,
    "rcsb_results": [],
    "rcsb_selected_index": 0,
    "downloaded_structure_bytes": None,
    "downloaded_structure_name": None,
    "downloaded_structure_pdb_id": None,
    "rcsb_last_query": "",
    "structure_source": "Upload local file",
    "pocket_charge": "negative",
    "pocket_hydrophobicity": "medium",
    "pocket_centroid": None,
    "last_run_structure_filename": None,  # centroid が有効な構造のファイル名
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
# 構造が変わっていたら古いcentroid（別構造のもの）を使わない
_current_filename = params.get("structure_filename")
_last_run_filename = st.session_state.get("last_run_structure_filename")
_valid_centroid = (
    st.session_state["pocket_centroid"]
    if _current_filename is not None and _current_filename == _last_run_filename
    else None
)

render_viewer_section(
    structure_bytes=params["structure_bytes"],
    structure_filename=params["structure_filename"],
    pdb_summary=params["pdb_summary"],
    result_df=st.session_state["result_df"] if _valid_centroid is not None else None,
    pocket_centroid=_valid_centroid,
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
            seed_sequence=params.get("seed_sequence"),
            variant_strategies=params.get("variant_strategies"),
        )
        st.session_state["result_df"] = result_df
        st.session_state["docking_done"] = False
        st.session_state["admet_done"] = False
        st.session_state["viewer_selected_rank"] = 1
        st.session_state["pocket_centroid"] = pocket_centroid
        st.session_state["last_run_structure_filename"] = params.get("structure_filename")
        # ポケット特性を session_state に保存（説明生成に使用）
        st.session_state["pocket_charge"] = params["pocket_charge"]
        st.session_state["pocket_hydrophobicity"] = params["pocket_hydrophobicity"]

        # 完了メッセージ
        n_candidates = len(result_df) if result_df is not None else 0
        mpnn_mode = "受容体条件付き ProteinMPNN" if pocket_centroid is not None else "ProteinMPNN (構造フリー)"
        _seed = params.get("seed_sequence")
        _gen_mode = f"バリアント生成（シード: {_seed}）" if _seed else "ランダム生成"
        st.success(
            f"完了 ✅  {n_candidates} 件の候補を生成しました。"
            f"  生成: {_gen_mode}。  スコアリング: {mpnn_mode}。  CSV: `{saved_path}`"
        )

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

    # ビューアの再描画（Selectivity含む全処理後に実行）
    # ビューアはページ上部で描画済みのため、pocket_centroid を反映するには再レンダリングが必要
    st.rerun()


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
# ADMET (Peptide Ranker)
# =========================
if (
    st.session_state["result_df"] is not None
    and not st.session_state["admet_done"]
):
    top_n_admet = min(20, len(st.session_state["result_df"]))
    if st.button(
        f"🧪 Score Bioactivity (top {top_n_admet})",
        type="secondary",
    ):
        from core.admet_scorer import score_top_candidates

        with st.spinner("Scoring bioactivity..."):
            result_df = score_top_candidates(
                st.session_state["result_df"],
                top_n=top_n_admet,
            )
            st.session_state["result_df"] = result_df
            st.session_state["admet_done"] = True
            n_scored = result_df["bioactivity_score"].notna().sum()
            st.success(
                f"生物活性スコア計算完了 ✅  上位 {n_scored} 件をスコアリングしました。"
            )


# =========================
# Result display
# =========================
render_results(
    st.session_state["result_df"],
    params["pdb_summary"],
    pocket_charge=st.session_state["pocket_charge"],
    pocket_hydrophobicity=st.session_state["pocket_hydrophobicity"],
)

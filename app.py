import streamlit as st

from core.pipeline import run_pipeline
from ui.sidebar import render_sidebar
from ui.results import render_results
from ui.structure_viewer import render_viewer_section


# =========================
# App config
# =========================
st.set_page_config(page_title="Peptide Discovery PoC", layout="wide")

st.title("Peptide Discovery PoC")
st.caption(
    "入力 → 候補生成 → フィルタ → 再スコア → 多様性制御 → 既知モチーフ比較 → 表示 の最小版"
)


# =========================
# Session state init
# =========================
defaults = {
    "result_df": None,
    "rcsb_results": [],
    "rcsb_selected_index": 0,
    "downloaded_structure_bytes": None,
    "downloaded_structure_name": None,
    "downloaded_structure_pdb_id": None,
    "rcsb_last_query": "",
    "structure_source": "Upload local file",
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
        )
        st.session_state["result_df"] = result_df
        st.success(f"Done. Saved CSV to: {saved_path}")


# =========================
# Result display
# =========================
render_results(st.session_state["result_df"], params["pdb_summary"])

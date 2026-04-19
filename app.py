import streamlit as st

from ui.actions import (
    render_bioactivity_button,
    render_docking_button,
    render_docking_selectivity_button,
    render_pipeline_run,
)
from ui.results import render_results
from ui.sidebar import render_sidebar
from ui.structure_viewer import render_viewer_section

# ─── Page config ─────────────────────────────────────────────────────────────
st.set_page_config(page_title="Peptide Discovery PoC", layout="wide", page_icon="🧬")

st.title("🧬 Peptide Discovery PoC")
st.caption(
    "Generates, scores, and visualizes peptide candidates from a protein pocket structure.  "
    "|  タンパク質ポケット構造を入力として、結合候補ペプチドを生成・スコアリング・可視化します。"
)

# ─── Session state defaults ───────────────────────────────────────────────────
_DEFAULTS: dict = {
    "result_df": None,
    "docking_done": False,
    "docking_selectivity_done": False,
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
    "last_run_structure_filename": None,
}
for _key, _val in _DEFAULTS.items():
    if _key not in st.session_state:
        st.session_state[_key] = _val

# ─── Sidebar ──────────────────────────────────────────────────────────────────
params = render_sidebar()

# ─── 3D Structure Viewer ──────────────────────────────────────────────────────
# Only pass the pocket centroid if it was computed for the currently loaded structure.
_current_fn = params.get("structure_filename")
_last_fn = st.session_state.get("last_run_structure_filename")
_valid_centroid = (
    st.session_state["pocket_centroid"]
    if _current_fn is not None and _current_fn == _last_fn
    else None
)
render_viewer_section(
    structure_bytes=params["structure_bytes"],
    structure_filename=params["structure_filename"],
    pdb_summary=params["pdb_summary"],
    result_df=st.session_state["result_df"] if _valid_centroid is not None else None,
    pocket_centroid=_valid_centroid,
)

# ─── Pipeline + optional post-processing ─────────────────────────────────────
render_pipeline_run(params)
render_docking_button(params)
render_docking_selectivity_button(params)
render_bioactivity_button()

# ─── Results table ────────────────────────────────────────────────────────────
render_results(
    st.session_state["result_df"],
    params["pdb_summary"],
    pocket_charge=st.session_state["pocket_charge"],
    pocket_hydrophobicity=st.session_state["pocket_hydrophobicity"],
)

from pathlib import Path
from datetime import datetime
import json
from urllib import request, error

import pandas as pd
import streamlit as st

from core.generator import generate_candidates
from core.filters import apply_filters, add_basic_properties
from core.rescorer import rescore_candidates
from core.utils import save_run_csv
from core.diversity import diversify_candidates
from core.pdb_utils import (
    load_structure_and_summary,
    summarize_structure_region,
    summarize_structure_ligand_pocket,
    get_chain_residue_numbers,
    get_ligand_options,
)
from core.motif_compare import (
    compare_candidates_to_known,
    get_default_motif_presets,
    parse_known_sequences,
)


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
if "result_df" not in st.session_state:
    st.session_state["result_df"] = None

if "rcsb_results" not in st.session_state:
    st.session_state["rcsb_results"] = []

if "rcsb_selected_index" not in st.session_state:
    st.session_state["rcsb_selected_index"] = 0

if "downloaded_structure_bytes" not in st.session_state:
    st.session_state["downloaded_structure_bytes"] = None

if "downloaded_structure_name" not in st.session_state:
    st.session_state["downloaded_structure_name"] = None

if "downloaded_structure_pdb_id" not in st.session_state:
    st.session_state["downloaded_structure_pdb_id"] = None

if "rcsb_last_query" not in st.session_state:
    st.session_state["rcsb_last_query"] = ""

if "structure_source" not in st.session_state:
    st.session_state["structure_source"] = "Upload local file"


# =========================
# Helpers
# =========================
class InMemoryUploadedStructure:
    def __init__(self, data: bytes, name: str):
        self._data = data
        self.name = name

    def getvalue(self):
        return self._data


def http_post_json(url: str, payload: dict, timeout: int = 20) -> dict:
    body = json.dumps(payload).encode("utf-8")
    req = request.Request(
        url,
        data=body,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    with request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def http_get_json(url: str, timeout: int = 20) -> dict:
    req = request.Request(url, headers={"Accept": "application/json"})
    with request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def http_get_bytes(url: str, timeout: int = 30) -> bytes:
    req = request.Request(url, headers={"Accept": "*/*"})
    with request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def safe_nested_get(d: dict, path: list, default=None):
    cur = d
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def normalize_text(x) -> str:
    return str(x).strip().lower()


def calc_target_match_score(target_label: str, title: str) -> float:
    target = normalize_text(target_label)
    title_l = normalize_text(title)

    if not target:
        return 0.0

    if target in title_l:
        return 1.0

    target_tokens = [x for x in target.replace("-", " ").split() if x]
    if target_tokens and any(tok in title_l for tok in target_tokens):
        return 0.6

    return 0.0


def calc_ligand_score(title: str) -> float:
    t = normalize_text(title)
    strong_keywords = ["inhibitor", "ligand", "bound", "complex"]
    if any(k in t for k in strong_keywords):
        return 1.0
    return 0.0


def calc_resolution_score(resolution) -> float:
    try:
        r = float(resolution)
    except Exception:
        return 0.3

    if r <= 2.0:
        return 1.0
    elif r <= 2.5:
        return 0.8
    elif r <= 3.0:
        return 0.6
    elif r <= 3.5:
        return 0.4
    else:
        return 0.2


def calc_method_score(method: str) -> float:
    m = normalize_text(method)

    if "x-ray" in m or "diffraction" in m:
        return 1.0
    elif "cryo" in m:
        return 0.8
    elif "nmr" in m:
        return 0.5
    elif m.strip():
        return 0.4
    else:
        return 0.3


def calc_domain_score(title: str) -> float:
    t = normalize_text(title)
    score = 0.5

    if "kinase" in t:
        score += 0.3
    if "inhibitor" in t:
        score += 0.2
    if "complex" in t:
        score += 0.1
    if "extracellular" in t:
        score -= 0.3

    return max(0.0, min(1.0, score))


def calc_structure_priority_score(record: dict, target_label: str) -> dict:
    title = record.get("title", "")
    method = record.get("method", "")
    resolution = record.get("resolution", None)

    target_match_score = calc_target_match_score(target_label, title)
    ligand_score = calc_ligand_score(title)
    resolution_score = calc_resolution_score(resolution)
    method_score = calc_method_score(method)
    domain_score = calc_domain_score(title)

    structure_priority_score = (
        0.35 * target_match_score
        + 0.25 * ligand_score
        + 0.20 * resolution_score
        + 0.10 * method_score
        + 0.10 * domain_score
    )

    out = dict(record)
    out["target_match_score"] = round(target_match_score, 3)
    out["ligand_score"] = round(ligand_score, 3)
    out["resolution_score"] = round(resolution_score, 3)
    out["method_score"] = round(method_score, 3)
    out["domain_score"] = round(domain_score, 3)
    out["structure_priority_score"] = round(structure_priority_score, 3)
    return out


def fetch_rcsb_entry_metadata(pdb_id: str) -> dict:
    pdb_id = str(pdb_id).lower()
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    try:
        data = http_get_json(url)
    except Exception:
        return {
            "pdb_id": pdb_id.upper(),
            "title": "",
            "method": "",
            "resolution": None,
            "deposit_date": "",
        }

    title = safe_nested_get(data, ["struct", "title"], "")
    exptl = data.get("exptl", [])
    method = ""
    if isinstance(exptl, list) and len(exptl) > 0:
        method = exptl[0].get("method", "")

    resolution = None
    rcsb_entry_info = data.get("rcsb_entry_info", {})
    if isinstance(rcsb_entry_info, dict):
        res_list = rcsb_entry_info.get("resolution_combined")
        if isinstance(res_list, list) and len(res_list) > 0:
            resolution = res_list[0]

    deposit_date = ""
    accession_info = data.get("rcsb_accession_info", {})
    if isinstance(accession_info, dict):
        deposit_date = accession_info.get("deposit_date", "")

    return {
        "pdb_id": pdb_id.upper(),
        "title": title,
        "method": method,
        "resolution": resolution,
        "deposit_date": deposit_date,
    }


def search_rcsb_structures(query_text: str, target_label: str, rows: int = 10) -> list[dict]:
    """
    Search RCSB by free text, enrich with entry metadata,
    then add structure priority score and sort descending.
    Handles both compact(string) and object-style result_set items safely.
    """
    q = str(query_text).strip()
    if not q:
        return []

    payload = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": q,
            },
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": rows,
            },
            "results_verbosity": "minimal",
            "sort": [
                {
                    "sort_by": "score",
                    "direction": "desc",
                }
            ],
        },
    }

    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    data = http_post_json(search_url, payload)

    result_set = data.get("result_set", [])
    pdb_ids = []

    for item in result_set:
        identifier = None

        if isinstance(item, str):
            identifier = item
        elif isinstance(item, dict):
            identifier = item.get("identifier")

        if identifier:
            pdb_ids.append(str(identifier).upper())

    out = []
    for pdb_id in pdb_ids:
        meta = fetch_rcsb_entry_metadata(pdb_id)
        meta = calc_structure_priority_score(meta, target_label=target_label)
        out.append(meta)

    out = sorted(
        out,
        key=lambda x: x.get("structure_priority_score", 0),
        reverse=True,
    )
    return out


def download_rcsb_mmcif(pdb_id: str) -> tuple[bytes, str]:
    pid = str(pdb_id).lower()
    url = f"https://files.rcsb.org/download/{pid}.cif"
    file_bytes = http_get_bytes(url)
    filename = f"{pid}.cif"
    return file_bytes, filename


def format_resolution(x):
    if x is None or x == "":
        return "-"
    try:
        return f"{float(x):.2f} Å"
    except Exception:
        return str(x)


def build_rcsb_label(record: dict) -> str:
    pdb_id = record.get("pdb_id", "-")
    title = record.get("title", "") or "-"
    resolution = format_resolution(record.get("resolution"))
    score = record.get("structure_priority_score", "-")
    return f"{pdb_id} | score={score} | {title} | {resolution}"


def reset_downloaded_structure():
    st.session_state["downloaded_structure_bytes"] = None
    st.session_state["downloaded_structure_name"] = None
    st.session_state["downloaded_structure_pdb_id"] = None


# =========================
# Main state
# =========================
pdb_summary = None
pdb_parse_error = None
structure = None
available_chains = []
selected_chain = None
uploaded_structure = None


# =========================
# Sidebar
# =========================
with st.sidebar:
    st.header("Input")

    target_name = st.text_input("Target label", value="EGFR")
    num_candidates = st.slider("Number of candidates", 100, 5000, 1000, 100)
    min_len = st.slider("Min peptide length", 4, 20, 7)
    max_len = st.slider("Max peptide length", 4, 30, 12)

    st.markdown("### Structure input")

    structure_source = st.radio(
        "Structure source",
        ["Upload local file", "Search RCSB by target"],
        index=0 if st.session_state["structure_source"] == "Upload local file" else 1,
        key="structure_source",
    )

    if structure_source == "Upload local file":
        uploaded_structure = st.file_uploader(
            "Upload PDB or mmCIF file",
            type=["pdb", "cif", "mmcif"],
            accept_multiple_files=False,
            key="uploaded_structure_local",
        )

        if uploaded_structure is not None:
            reset_downloaded_structure()

    else:
        rcsb_query = st.text_input(
            "Search keyword",
            value=target_name if not st.session_state["rcsb_last_query"] else st.session_state["rcsb_last_query"],
            key="rcsb_query_text",
        )

        search_col1, search_col2 = st.columns([1, 1])
        with search_col1:
            search_button = st.button("Search RCSB", key="search_rcsb_button", use_container_width=True)
        with search_col2:
            clear_button = st.button("Clear search", key="clear_rcsb_button", use_container_width=True)

        if clear_button:
            st.session_state["rcsb_results"] = []
            st.session_state["rcsb_selected_index"] = 0
            st.session_state["rcsb_last_query"] = ""
            reset_downloaded_structure()

        if search_button:
            try:
                st.session_state["rcsb_last_query"] = rcsb_query
                results = search_rcsb_structures(
                    query_text=rcsb_query,
                    target_label=target_name,
                    rows=10,
                )
                st.session_state["rcsb_results"] = results
                st.session_state["rcsb_selected_index"] = 0
                reset_downloaded_structure()
                if not results:
                    st.info("RCSBで候補が見つかりませんでした。検索語を変えて試してください。")
            except error.HTTPError as e:
                st.error(f"RCSB検索でHTTPエラー: {e}")
            except Exception as e:
                st.error(f"RCSB検索でエラー: {type(e).__name__}: {e}")

        rcsb_results = st.session_state["rcsb_results"]

        if rcsb_results:
            result_labels = [build_rcsb_label(x) for x in rcsb_results]
            selected_label = st.selectbox(
                "Search results",
                options=result_labels,
                index=min(st.session_state["rcsb_selected_index"], len(result_labels) - 1),
                key="rcsb_result_selectbox",
            )
            st.session_state["rcsb_selected_index"] = result_labels.index(selected_label)
            selected_record = rcsb_results[st.session_state["rcsb_selected_index"]]

            st.caption(
                f"priority={selected_record.get('structure_priority_score', '-')}, "
                f"target={selected_record.get('target_match_score', '-')}, "
                f"ligand={selected_record.get('ligand_score', '-')}, "
                f"resolution={selected_record.get('resolution_score', '-')}, "
                f"method={selected_record.get('method_score', '-')}, "
                f"domain={selected_record.get('domain_score', '-')}"
            )

            if st.button("Use selected structure", key="use_selected_rcsb_structure", use_container_width=True):
                try:
                    file_bytes, filename = download_rcsb_mmcif(selected_record["pdb_id"])
                    st.session_state["downloaded_structure_bytes"] = file_bytes
                    st.session_state["downloaded_structure_name"] = filename
                    st.session_state["downloaded_structure_pdb_id"] = selected_record["pdb_id"]
                    st.success(f"Loaded structure: {selected_record['pdb_id']}")
                except error.HTTPError as e:
                    st.error(f"構造ダウンロードでHTTPエラー: {e}")
                except Exception as e:
                    st.error(f"構造ダウンロードでエラー: {type(e).__name__}: {e}")

        if st.session_state["downloaded_structure_bytes"] is not None:
            uploaded_structure = InMemoryUploadedStructure(
                st.session_state["downloaded_structure_bytes"],
                st.session_state["downloaded_structure_name"],
            )
            st.caption(
                f"Current selected structure: "
                f"{st.session_state['downloaded_structure_pdb_id']} "
                f"({st.session_state['downloaded_structure_name']})"
            )

    use_structure_features = st.checkbox(
        "Use selected structure to auto-fill pocket bias",
        value=True,
    )

    pocket_mode = st.radio(
        "Pocket source",
        ["Manual region", "Ligand neighborhood"],
        index=0,
    )

    if uploaded_structure is not None:
        try:
            structure, initial_summary = load_structure_and_summary(
                uploaded_structure.getvalue(),
                filename=uploaded_structure.name,
                structure_id=uploaded_structure.name,
            )

            available_chains = initial_summary["chains"]

            if available_chains:
                selected_chain = st.selectbox("Select chain", available_chains, index=0)

                if pocket_mode == "Manual region":
                    available_residue_numbers = get_chain_residue_numbers(structure, selected_chain)

                    if available_residue_numbers:
                        min_res = min(available_residue_numbers)
                        max_res = max(available_residue_numbers)

                        residue_range = st.slider(
                            "Residue range",
                            min_value=min_res,
                            max_value=max_res,
                            value=(min_res, max_res),
                        )

                        residue_start, residue_end = residue_range

                        _, pdb_summary = summarize_structure_region(
                            uploaded_structure.getvalue(),
                            filename=uploaded_structure.name,
                            structure_id=uploaded_structure.name,
                            selected_chain_id=selected_chain,
                            residue_start=residue_start,
                            residue_end=residue_end,
                        )
                    else:
                        _, pdb_summary = summarize_structure_region(
                            uploaded_structure.getvalue(),
                            filename=uploaded_structure.name,
                            structure_id=uploaded_structure.name,
                            selected_chain_id=selected_chain,
                        )

                else:
                    ligand_options = get_ligand_options(structure, selected_chain_id=selected_chain)

                    if ligand_options:
                        ligand_labels = [x["label"] for x in ligand_options]
                        selected_ligand_label = st.selectbox("Select ligand", ligand_labels, index=0)
                        selected_ligand = next(x for x in ligand_options if x["label"] == selected_ligand_label)

                        pocket_radius = st.slider(
                            "Pocket radius (Å)",
                            min_value=3.0,
                            max_value=12.0,
                            value=6.0,
                            step=0.5,
                        )

                        _, pdb_summary = summarize_structure_ligand_pocket(
                            uploaded_structure.getvalue(),
                            filename=uploaded_structure.name,
                            structure_id=uploaded_structure.name,
                            ligand_chain_id=selected_ligand["chain_id"],
                            ligand_residue_id=selected_ligand["residue_id"],
                            radius=pocket_radius,
                        )
                    else:
                        st.warning("No ligand-like HETATM residues found in the selected chain.")
                        pdb_summary = initial_summary
            else:
                pdb_summary = initial_summary

        except Exception as e:
            pdb_parse_error = str(e)

    charge_options = ["neutral", "negative", "positive"]
    hydrophobicity_options = ["low", "medium", "high"]

    default_pocket_charge = "negative"
    default_pocket_hydrophobicity = "medium"

    if pdb_summary is not None and use_structure_features:
        default_pocket_charge = pdb_summary["pocket_charge_guess"]
        default_pocket_hydrophobicity = pdb_summary["pocket_hydrophobicity_guess"]

    st.markdown("### Pocket bias")
    pocket_charge = st.selectbox(
        "Pocket charge tendency",
        charge_options,
        index=charge_options.index(default_pocket_charge),
    )
    pocket_hydrophobicity = st.selectbox(
        "Pocket hydrophobicity",
        hydrophobicity_options,
        index=hydrophobicity_options.index(default_pocket_hydrophobicity),
    )

    st.markdown("### Constraints")
    avoid_cysteine = st.checkbox("Avoid Cys (C)", value=True)
    max_abs_charge = st.slider("Max absolute net charge", 0, 8, 4)
    max_hydrophobicity = st.slider("Max average hydrophobicity", 0.0, 1.0, 0.65, 0.01)
    max_repeat_residue = st.slider("Max repeated same residue in a row", 1, 5, 2)

    st.markdown("### Rescoring preferences")
    preferred_len_min = st.slider("Preferred length (min)", 4, 20, 8)
    preferred_len_max = st.slider("Preferred length (max)", 4, 30, 12)

    st.markdown("### Diversity control")
    use_diversity_filter = st.checkbox("Enable diversity filter", value=True)
    min_diversity_distance = st.slider(
        "Min normalized sequence distance",
        min_value=0.0,
        max_value=1.0,
        value=0.30,
        step=0.05,
    )
    max_diverse_candidates = st.slider(
        "Max diverse candidates to keep",
        min_value=10,
        max_value=500,
        value=100,
        step=10,
    )

    st.markdown("### Known ligand / motif comparison")
    preset_dict = get_default_motif_presets()
    preset_labels = list(preset_dict.keys())

    selected_presets = st.multiselect(
        "Preset motifs",
        options=preset_labels,
        default=[],
        help="PoC向け参照モチーフ。厳密な既知 binder 保証ではなく説明性向上用です。",
    )

    custom_known_sequences_text = st.text_area(
        "Custom known sequences (one per line / comma / semicolon / space separated)",
        value="",
        height=120,
        help="例: KLAKLAKKLAKLAK",
    )

    run_button = st.button("Generate and Filter", type="primary")

    if pdb_parse_error is not None:
        st.error(f"Failed to parse structure: {pdb_parse_error}")
    elif pdb_summary is not None:
        st.markdown("### Parsed structure summary")
        st.write(f"Format: {pdb_summary.get('file_format', '-')}")
        if pdb_summary.get("selected_chain") is not None:
            st.write(f"Selected chain: {pdb_summary['selected_chain']}")
        if pdb_summary.get("source_mode") is not None:
            st.write(f"Source mode: {pdb_summary['source_mode']}")
        if pdb_summary.get("residue_start") is not None and pdb_summary.get("residue_end") is not None:
            st.write(f"Residue range: {pdb_summary['residue_start']} - {pdb_summary['residue_end']}")
        if pdb_summary.get("ligand_names"):
            st.write(f"Ligand(s): {', '.join(pdb_summary['ligand_names'])}")
        if pdb_summary.get("search_radius") is not None:
            st.write(f"Pocket radius: {pdb_summary['search_radius']} Å")
        st.write(f"Chains: {', '.join(pdb_summary['chains']) if pdb_summary['chains'] else '(none)'}")
        st.write(f"Residues counted: {pdb_summary['residue_count']}")
        st.write(f"Atoms counted: {pdb_summary['atom_count']}")
        st.write(f"Auto charge guess: {pdb_summary['pocket_charge_guess']}")
        st.write(f"Auto hydrophobicity guess: {pdb_summary['pocket_hydrophobicity_guess']}")
    else:
        st.caption("構造を選択すると parsed structure summary が表示されます。")


# =========================
# Run generation pipeline
# =========================
if run_button:
    with st.spinner("Generating candidates..."):
        preset_sequences = [preset_dict[label] for label in selected_presets]
        custom_sequences = parse_known_sequences(custom_known_sequences_text)

        known_sequences = []
        seen = set()
        for seq in preset_sequences + custom_sequences:
            if seq not in seen:
                seen.add(seq)
                known_sequences.append(seq)

        candidates = generate_candidates(
            n=num_candidates,
            min_len=min_len,
            max_len=max_len,
            pocket_charge=pocket_charge,
            pocket_hydrophobicity=pocket_hydrophobicity,
            avoid_residues=["C"] if avoid_cysteine else [],
        )

        df = pd.DataFrame(candidates)
        df = add_basic_properties(df)

        filtered_df = apply_filters(
            df,
            min_len=min_len,
            max_len=max_len,
            max_abs_charge=max_abs_charge,
            max_hydrophobicity=max_hydrophobicity,
            max_repeat_residue=max_repeat_residue,
            remove_near_duplicates=True,
        )

        rescored_df = rescore_candidates(
            filtered_df,
            pocket_charge=pocket_charge,
            pocket_hydrophobicity=pocket_hydrophobicity,
            preferred_len_min=preferred_len_min,
            preferred_len_max=preferred_len_max,
        )

        rescored_df = rescored_df.sort_values(
            by=["final_score", "rescoring_score", "property_score", "gen_score"],
            ascending=False,
        ).reset_index(drop=True)

        if use_diversity_filter:
            rescored_df = diversify_candidates(
                rescored_df,
                sequence_col="sequence",
                score_col="final_score",
                min_normalized_distance=min_diversity_distance,
                max_candidates=max_diverse_candidates,
            )
        else:
            rescored_df = rescored_df.copy()
            rescored_df["diversity_kept"] = True
            if "diversity_min_distance" not in rescored_df.columns:
                rescored_df["diversity_min_distance"] = 1.0

        rescored_df = compare_candidates_to_known(
            rescored_df,
            known_sequences=known_sequences,
        )

        rescored_df = rescored_df.reset_index(drop=True)
        rescored_df.insert(0, "rank", range(1, len(rescored_df) + 1))
        st.session_state["result_df"] = rescored_df

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path("outputs")
        output_file = output_dir / f"{target_name}_{timestamp}.csv"
        saved_path = save_run_csv(rescored_df, str(output_file))

        st.success(f"Done. Saved CSV to: {saved_path}")


# =========================
# Result display
# =========================
result_df = st.session_state["result_df"]

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
        "rank",
        "sequence",
        "final_score",
        "rescoring_score",
        "gen_score",
        "property_score",
        "length",
        "net_charge",
        "avg_hydrophobicity",
        "charge_match_score",
        "hydrophobic_match_score",
        "diversity_kept",
        "diversity_min_distance",
        "aromatic_ratio",
        "charge_pattern",
        "best_known_motif",
        "known_aromatic_gap_min",
        "known_charge_pattern_similarity_max",
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
        if "complexity_score" in row.index and pd.notna(row["complexity_score"]):
            st.write(f"- Complexity score: {row['complexity_score']:.3f}")

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
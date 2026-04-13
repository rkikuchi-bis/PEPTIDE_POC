from urllib import error

import streamlit as st

from core.rcsb_client import (
    search_rcsb_structures,
    download_rcsb_mmcif,
    build_rcsb_label,
    InMemoryUploadedStructure,
)
from core.pdb_utils import (
    load_structure_and_summary,
    summarize_structure_region,
    summarize_structure_ligand_pocket,
    get_chain_residue_numbers,
    get_ligand_options,
    get_recommended_chain,
)
from core.motif_compare import (
    get_default_motif_presets,
    parse_known_sequences,
)


def reset_downloaded_structure():
    st.session_state["downloaded_structure_bytes"] = None
    st.session_state["downloaded_structure_name"] = None
    st.session_state["downloaded_structure_pdb_id"] = None


def render_sidebar() -> dict:
    pdb_summary = None
    pdb_parse_error = None
    structure = None
    available_chains = []
    selected_chain = None
    uploaded_structure = None

    # Simple mode のデフォルト値（後で上書きされる場合あり）
    num_candidates = 500
    min_len = 6
    max_len = 12
    avoid_cysteine = True
    max_abs_charge = 4
    max_hydrophobicity = 0.65
    max_repeat_residue = 2
    preferred_len_min = 6
    preferred_len_max = 10
    use_diversity_filter = True
    min_diversity_distance = 0.30
    max_diverse_candidates = 100
    known_sequences = []
    enable_docking = False
    docking_params = None
    pocket_charge = "negative"
    pocket_hydrophobicity = "medium"
    charge_options = ["neutral", "negative", "positive"]
    hydrophobicity_options = ["low", "medium", "high"]

    with st.sidebar:

        # ── Mode toggle ──────────────────────────────
        mode = st.radio(
            "Mode",
            ["Simple", "Expert"],
            horizontal=True,
            key="sidebar_mode",
            help=(
                "Simple: 構造ファイルをアップロードして Run するだけ。\n"
                "Expert: 全パラメータを手動で調整できます。"
            ),
        )
        st.divider()

        # ── Project name ─────────────────────────────
        target_name = st.text_input(
            "Project name",
            value="",
            placeholder="例: MDM2_inhibitor（空白可）",
        )
        st.caption("結果CSVのファイル名に使われます（空白可）")

        # ── Structure input（Simple / Expert 共通）───
        st.markdown("### Structure")

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
                value=(
                    target_name
                    if target_name and not st.session_state["rcsb_last_query"]
                    else st.session_state["rcsb_last_query"]
                ),
                key="rcsb_query_text",
            )
            search_col1, search_col2 = st.columns([1, 1])
            with search_col1:
                search_button = st.button("Search RCSB", key="search_rcsb_button", width="stretch")
            with search_col2:
                clear_button = st.button("Clear search", key="clear_rcsb_button", width="stretch")

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
                        st.info("No results found on RCSB.")
                except error.HTTPError as e:
                    st.error(f"RCSB search HTTP error: {e}")
                except Exception as e:
                    st.error(f"RCSB search error: {type(e).__name__}: {e}")

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
                if st.button("Use selected structure", key="use_selected_rcsb_structure", width="stretch"):
                    try:
                        file_bytes, filename = download_rcsb_mmcif(selected_record["pdb_id"])
                        st.session_state["downloaded_structure_bytes"] = file_bytes
                        st.session_state["downloaded_structure_name"] = filename
                        st.session_state["downloaded_structure_pdb_id"] = selected_record["pdb_id"]
                        st.success(f"Loaded: {selected_record['pdb_id']}")
                    except error.HTTPError as e:
                        st.error(f"Download error: {e}")
                    except Exception as e:
                        st.error(f"Download error: {type(e).__name__}: {e}")

            if st.session_state["downloaded_structure_bytes"] is not None:
                uploaded_structure = InMemoryUploadedStructure(
                    st.session_state["downloaded_structure_bytes"],
                    st.session_state["downloaded_structure_name"],
                )
                st.caption(
                    f"Current: {st.session_state['downloaded_structure_pdb_id']} "
                    f"({st.session_state['downloaded_structure_name']})"
                )

        # ── Structure parsing ────────────────────────
        if uploaded_structure is not None:
            try:
                structure, initial_summary = load_structure_and_summary(
                    uploaded_structure.getvalue(),
                    filename=uploaded_structure.name,
                    structure_id=uploaded_structure.name,
                )
                available_chains = initial_summary["chains"]

                if available_chains:

                    if mode == "Simple":
                        # チェーン自動推奨（最多AA残基）
                        recommended_chain = get_recommended_chain(structure)
                        default_chain_idx = (
                            available_chains.index(recommended_chain)
                            if recommended_chain in available_chains
                            else 0
                        )
                        selected_chain = st.selectbox(
                            "Chain",
                            available_chains,
                            index=default_chain_idx,
                        )
                        if recommended_chain and recommended_chain == selected_chain:
                            st.caption(f"Chain {recommended_chain}（自動推奨: 最多残基）")

                        # ポケット自動選択
                        ligand_options = get_ligand_options(structure, selected_chain_id=selected_chain)
                        if ligand_options:
                            lig_labels = [x["label"] for x in ligand_options]
                            selected_lig_label = st.selectbox("Ligand", lig_labels, index=0)
                            selected_lig = next(
                                x for x in ligand_options if x["label"] == selected_lig_label
                            )
                            _, pdb_summary = summarize_structure_ligand_pocket(
                                uploaded_structure.getvalue(),
                                filename=uploaded_structure.name,
                                structure_id=uploaded_structure.name,
                                ligand_chain_id=selected_lig["chain_id"],
                                ligand_residue_id=selected_lig["residue_id"],
                                radius=6.0,
                            )
                            st.caption(f"Pocket: Ligand neighborhood ({selected_lig_label}, 6Å)")
                        else:
                            _, pdb_summary = summarize_structure_region(
                                uploaded_structure.getvalue(),
                                filename=uploaded_structure.name,
                                structure_id=uploaded_structure.name,
                                selected_chain_id=selected_chain,
                            )
                            st.caption("Pocket: Full chain（リガンドなし）")

                        if pdb_summary is not None:
                            pocket_charge = pdb_summary["pocket_charge_guess"]
                            pocket_hydrophobicity = pdb_summary["pocket_hydrophobicity_guess"]
                            st.info(
                                f"Auto-detected pocket:  \n"
                                f"Charge → **{pocket_charge}**  \n"
                                f"Hydrophobicity → **{pocket_hydrophobicity}**"
                            )

                    else:  # Expert mode
                        selected_chain = st.selectbox("Select chain", available_chains, index=0)

                        pocket_mode = st.radio(
                            "Pocket source",
                            ["Manual region", "Ligand neighborhood"],
                            index=0,
                        )

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
                                selected_ligand = next(
                                    x for x in ligand_options if x["label"] == selected_ligand_label
                                )
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

        # ── Expert mode パラメータ ────────────────────
        if mode == "Expert":
            use_structure_features = st.checkbox(
                "Use selected structure to auto-fill pocket bias",
                value=True,
            )

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

            st.markdown("### Generation")
            num_candidates = st.slider("Number of candidates", 100, 5000, 1000, 100)
            min_len = st.slider("Min peptide length", 4, 20, 7)
            max_len = st.slider("Max peptide length", 4, 30, 12)

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
                min_value=0.0, max_value=1.0, value=0.30, step=0.05,
            )
            max_diverse_candidates = st.slider(
                "Max diverse candidates to keep",
                min_value=10, max_value=500, value=100, step=10,
            )

            st.markdown("### Known ligand / motif comparison")
            preset_dict = get_default_motif_presets()
            preset_labels = list(preset_dict.keys())
            selected_presets = st.multiselect(
                "Preset motifs",
                options=preset_labels,
                default=[],
                help="Reference motifs for comparison. / PoC向け参照モチーフ。厳密な既知 binder 保証ではなく説明性向上用です。",
            )
            custom_known_sequences_text = st.text_area(
                "Custom known sequences (one per line / comma / semicolon / space separated)",
                value="",
                height=120,
            )
            preset_sequences = [preset_dict[label] for label in selected_presets]
            custom_sequences = parse_known_sequences(custom_known_sequences_text)
            known_sequences = list(dict.fromkeys(preset_sequences + custom_sequences))

            st.markdown("### Docking (Phase B-1)")
            from core.docking import is_vina_available, compute_pocket_center
            vina_ok = is_vina_available()
            if not vina_ok:
                st.caption("Docking unavailable: vina binary not found (bin/vina)")
                enable_docking = False
                docking_params = None
            else:
                enable_docking = st.checkbox("Enable AutoDock Vina docking", value=False)
                if enable_docking:
                    auto_center = None
                    if uploaded_structure is not None and pdb_summary is not None:
                        auto_center = compute_pocket_center(
                            uploaded_structure.getvalue(),
                            uploaded_structure.name,
                            pdb_summary,
                        )
                    default_cx = auto_center[0] if auto_center else 0.0
                    default_cy = auto_center[1] if auto_center else 0.0
                    default_cz = auto_center[2] if auto_center else 0.0
                    if auto_center:
                        st.caption(f"Auto-computed pocket center: ({default_cx}, {default_cy}, {default_cz})")
                    cx = st.number_input("Box center X [Å]", value=default_cx, step=1.0, format="%.2f")
                    cy = st.number_input("Box center Y [Å]", value=default_cy, step=1.0, format="%.2f")
                    cz = st.number_input("Box center Z [Å]", value=default_cz, step=1.0, format="%.2f")
                    box_sx = st.slider("Box size X [Å]", 10, 40, 20)
                    box_sy = st.slider("Box size Y [Å]", 10, 40, 20)
                    box_sz = st.slider("Box size Z [Å]", 10, 40, 20)
                    docking_top_n = st.slider("Top N candidates to dock", 5, 30, 10)
                    docking_exhaustiveness = st.select_slider(
                        "Exhaustiveness", options=[4, 8, 16, 32], value=8,
                        help="Higher = more accurate but slower.",
                    )
                    docking_params = {
                        "box_center": (cx, cy, cz),
                        "box_size": (float(box_sx), float(box_sy), float(box_sz)),
                        "top_n": docking_top_n,
                        "exhaustiveness": docking_exhaustiveness,
                    }
                else:
                    docking_params = None

        # ── Selectivity（Expert mode 専用）─────────────
        selectivity_params = None
        if mode == "Expert":
            st.markdown("### Selectivity")
            enable_selectivity = st.checkbox(
                "Enable Selectivity Mode",
                value=False,
                help=(
                    "Selectivity score = target fit − off-target fit.\n"
                    "選択性スコア = ターゲット適合 − オフターゲット適合\n\n"
                    "事前に知られたOff-targetに対し、結合選択性を定量評価します。"
                ),
            )

            if enable_selectivity:
                from core.offtarget_db import OFFTARGET_DB, list_target_keys

                offtarget_structure_bytes = None
                offtarget_structure_filename = None
                pdb_summary_offtarget = None
                offtarget_charge = "neutral"
                offtarget_hydrophobicity = "medium"
                offtarget_label = "Off-target"

                # Off-target 指定方法の選択
                ot_source = st.radio(
                    "Off-target source",
                    ["Known off-target DB から選択", "構造ファイルをアップロード"],
                    key="ot_source_mode",
                )

                if ot_source == "Known off-target DB から選択":
                    target_keys = list_target_keys()
                    selected_target_key = st.selectbox(
                        "Target を選択",
                        target_keys,
                        help="登録済みターゲットに対応する既知 Off-target が表示されます。",
                        key="ot_db_target_key",
                    )
                    ot_candidates = OFFTARGET_DB[selected_target_key]["offtargets"]
                    ot_labels = [x["label"] for x in ot_candidates]
                    selected_ot_label = st.selectbox(
                        "Off-target を選択",
                        ot_labels,
                        key="ot_db_offtarget_label",
                    )
                    selected_ot = next(x for x in ot_candidates if x["label"] == selected_ot_label)
                    st.caption(selected_ot["description"])

                    offtarget_label = selected_ot["label"]

                    if st.button("Off-target 構造を RCSB から取得", key="ot_db_download"):
                        try:
                            file_bytes, filename = download_rcsb_mmcif(selected_ot["pdb_id"])
                            st.session_state["ot_db_bytes"] = file_bytes
                            st.session_state["ot_db_filename"] = filename
                            st.success(f"Downloaded: {selected_ot['pdb_id']}")
                        except Exception as e:
                            st.error(f"Download failed: {e}")

                    # ダウンロード済みの場合は解析
                    if st.session_state.get("ot_db_bytes") is not None:
                        offtarget_structure_bytes = st.session_state["ot_db_bytes"]
                        offtarget_structure_filename = st.session_state["ot_db_filename"]
                        try:
                            ot_struct, ot_init_summary = load_structure_and_summary(
                                offtarget_structure_bytes,
                                filename=offtarget_structure_filename,
                                structure_id=offtarget_structure_filename,
                            )
                            ot_chains = ot_init_summary["chains"]
                            if ot_chains:
                                ot_recommended = get_recommended_chain(ot_struct)
                                ot_chain_idx = (
                                    ot_chains.index(ot_recommended)
                                    if ot_recommended in ot_chains else 0
                                )
                                ot_chain = st.selectbox(
                                    "Off-target chain",
                                    ot_chains,
                                    index=ot_chain_idx,
                                    key="ot_db_chain",
                                )
                                st.caption(f"推奨: Chain {ot_recommended}（最多残基）")

                                ot_ligands = get_ligand_options(ot_struct, selected_chain_id=ot_chain)
                                if ot_ligands:
                                    ot_lig_labels = [x["label"] for x in ot_ligands]
                                    sel_ot_lig_label = st.selectbox(
                                        "Off-target ligand", ot_lig_labels, index=0,
                                        key="ot_db_ligand",
                                    )
                                    sel_ot_lig = next(
                                        x for x in ot_ligands if x["label"] == sel_ot_lig_label
                                    )
                                    _, pdb_summary_offtarget = summarize_structure_ligand_pocket(
                                        offtarget_structure_bytes,
                                        filename=offtarget_structure_filename,
                                        structure_id=offtarget_structure_filename,
                                        ligand_chain_id=sel_ot_lig["chain_id"],
                                        ligand_residue_id=sel_ot_lig["residue_id"],
                                        radius=6.0,
                                    )
                                else:
                                    _, pdb_summary_offtarget = summarize_structure_region(
                                        offtarget_structure_bytes,
                                        filename=offtarget_structure_filename,
                                        structure_id=offtarget_structure_filename,
                                        selected_chain_id=ot_chain,
                                    )

                                if pdb_summary_offtarget is not None:
                                    offtarget_charge = pdb_summary_offtarget["pocket_charge_guess"]
                                    offtarget_hydrophobicity = pdb_summary_offtarget["pocket_hydrophobicity_guess"]
                                    st.info(
                                        f"Off-target auto-detected:  \n"
                                        f"Charge → **{offtarget_charge}**  \n"
                                        f"Hydrophobicity → **{offtarget_hydrophobicity}**"
                                    )
                        except Exception as e:
                            st.error(f"Off-target parse error: {e}")

                else:  # 構造ファイルをアップロード
                    offtarget_label = st.text_input("Off-target label", value="Off-target", key="ot_upload_label")
                    offtarget_uploaded = st.file_uploader(
                        "Off-target PDB or mmCIF",
                        type=["pdb", "cif", "mmcif"],
                        key="offtarget_structure_file",
                    )
                    if offtarget_uploaded is not None:
                        offtarget_structure_bytes = offtarget_uploaded.getvalue()
                        offtarget_structure_filename = offtarget_uploaded.name
                        try:
                            offtarget_struct, offtarget_initial_summary = load_structure_and_summary(
                                offtarget_structure_bytes,
                                filename=offtarget_uploaded.name,
                                structure_id=offtarget_uploaded.name,
                            )
                            offtarget_chains = offtarget_initial_summary["chains"]
                            if offtarget_chains:
                                ot_recommended = get_recommended_chain(offtarget_struct)
                                ot_chain_idx = (
                                    offtarget_chains.index(ot_recommended)
                                    if ot_recommended in offtarget_chains else 0
                                )
                                offtarget_chain = st.selectbox(
                                    "Off-target chain", offtarget_chains,
                                    index=ot_chain_idx, key="offtarget_chain",
                                )
                                st.caption(f"推奨: Chain {ot_recommended}（最多残基）")
                                offtarget_pocket_mode = st.radio(
                                    "Off-target pocket source",
                                    ["Ligand neighborhood", "Manual region"],
                                    index=0, key="offtarget_pocket_mode",
                                )
                                if offtarget_pocket_mode == "Manual region":
                                    offtarget_res_nums = get_chain_residue_numbers(offtarget_struct, offtarget_chain)
                                    if offtarget_res_nums:
                                        min_res_ot = min(offtarget_res_nums)
                                        max_res_ot = max(offtarget_res_nums)
                                        offtarget_res_range = st.slider(
                                            "Off-target residue range",
                                            min_value=min_res_ot, max_value=max_res_ot,
                                            value=(min_res_ot, max_res_ot),
                                            key="offtarget_residue_range",
                                        )
                                        _, pdb_summary_offtarget = summarize_structure_region(
                                            offtarget_structure_bytes,
                                            filename=offtarget_uploaded.name,
                                            structure_id=offtarget_uploaded.name,
                                            selected_chain_id=offtarget_chain,
                                            residue_start=offtarget_res_range[0],
                                            residue_end=offtarget_res_range[1],
                                        )
                                    else:
                                        _, pdb_summary_offtarget = summarize_structure_region(
                                            offtarget_structure_bytes,
                                            filename=offtarget_uploaded.name,
                                            structure_id=offtarget_uploaded.name,
                                            selected_chain_id=offtarget_chain,
                                        )
                                else:  # Ligand neighborhood
                                    offtarget_ligand_opts = get_ligand_options(
                                        offtarget_struct, selected_chain_id=offtarget_chain
                                    )
                                    if offtarget_ligand_opts:
                                        ot_lig_labels = [x["label"] for x in offtarget_ligand_opts]
                                        sel_ot_lig_label = st.selectbox(
                                            "Off-target ligand", ot_lig_labels, index=0,
                                            key="offtarget_ligand",
                                        )
                                        sel_ot_lig = next(
                                            x for x in offtarget_ligand_opts if x["label"] == sel_ot_lig_label
                                        )
                                        ot_pocket_radius = st.slider(
                                            "Off-target pocket radius (Å)",
                                            min_value=3.0, max_value=12.0, value=6.0, step=0.5,
                                            key="offtarget_pocket_radius",
                                        )
                                        _, pdb_summary_offtarget = summarize_structure_ligand_pocket(
                                            offtarget_structure_bytes,
                                            filename=offtarget_uploaded.name,
                                            structure_id=offtarget_uploaded.name,
                                            ligand_chain_id=sel_ot_lig["chain_id"],
                                            ligand_residue_id=sel_ot_lig["residue_id"],
                                            radius=ot_pocket_radius,
                                        )
                                    else:
                                        st.warning("No ligand-like HETATM residues found.")
                                        pdb_summary_offtarget = offtarget_initial_summary

                                if pdb_summary_offtarget is not None:
                                    offtarget_charge = pdb_summary_offtarget["pocket_charge_guess"]
                                    offtarget_hydrophobicity = pdb_summary_offtarget["pocket_hydrophobicity_guess"]
                                    st.info(
                                        f"Off-target auto-detected:  \n"
                                        f"Charge → **{offtarget_charge}**  \n"
                                        f"Hydrophobicity → **{offtarget_hydrophobicity}**"
                                    )
                        except Exception as e:
                            st.error(f"Off-target structure parse error: {e}")

                # pocket bias 表示（自動入力済み、手動上書き可）
                st.caption("Selectivity score = target rescoring − off-target rescoring  ｜  手動上書き可")
                offtarget_charge = st.selectbox(
                    "Off-target pocket charge", charge_options,
                    index=charge_options.index(offtarget_charge),
                    key="offtarget_charge",
                )
                offtarget_hydrophobicity = st.selectbox(
                    "Off-target pocket hydrophobicity", hydrophobicity_options,
                    index=hydrophobicity_options.index(offtarget_hydrophobicity),
                    key="offtarget_hydrophobicity",
                )

                st.markdown("**Selectivity weight λ** — Recommended: 0.3")
                selectivity_lambda = st.slider(
                    "λ", min_value=0.0, max_value=1.0, value=0.3, step=0.05,
                    help=(
                        "Ranking score = final_score + λ × selectivity_score\n\n"
                        "λ = 0.0: rank by binding affinity only (off-target risk ignored).\n"
                        "Higher λ: promotes candidates with low off-target affinity.\n"
                        "λ = 1.0: binding affinity and selectivity weighted equally.\n"
                        "Recommended 0.3: balanced — influences ranking without large reversals.\n\n"
                        "ランキングスコア = final_score + λ × selectivity_score"
                    ),
                )
                lambda_label = (
                    "Binding affinity only (selectivity ignored)"
                    if selectivity_lambda == 0.0
                    else f"Selectivity applied (λ={selectivity_lambda:.2f})"
                )
                st.caption(f"Current: {lambda_label}")

                selectivity_params = {
                    "offtarget_label": offtarget_label,
                    "pocket_charge_offtarget": offtarget_charge,
                    "pocket_hydrophobicity_offtarget": offtarget_hydrophobicity,
                    "selectivity_lambda": selectivity_lambda,
                    "offtarget_structure_bytes": offtarget_structure_bytes,
                    "offtarget_structure_filename": offtarget_structure_filename,
                    "pdb_summary_offtarget": pdb_summary_offtarget,
                }

        # ── Seed sequence（バリアント生成）────────────
        st.markdown("### Seed Sequence")
        seed_sequence_input = st.text_input(
            "Known binding sequence (optional)",
            value="",
            placeholder="例: KLAKLAK",
            key="seed_sequence_input",
            help=(
                "既知の結合ペプチド配列を入力すると、ランダム生成の代わりにバリアントを網羅生成します。\n\n"
                "If a known binding sequence is provided, variants are exhaustively generated "
                "instead of random sequences."
            ),
        )
        seed_sequence = None
        variant_strategies = None
        if seed_sequence_input.strip():
            from core.variant_generator import validate_seed_sequence, estimate_variant_count
            is_valid, err_msg = validate_seed_sequence(seed_sequence_input)
            if is_valid:
                seed_sequence = seed_sequence_input.strip().upper()
                if mode == "Expert":
                    _strategy_labels = {
                        "single_mutant": "Single mutant scan — 各位置を全アミノ酸に変換",
                        "alanine_scan":  "Alanine scan — 各位置をAlaに置換",
                        "truncation":    "Truncation — N/C末端から最大3残基削る",
                        "scramble":      "Scramble — ランダムシャッフル（対照用）",
                    }
                    selected_strategies = st.multiselect(
                        "Variant strategies",
                        options=list(_strategy_labels.keys()),
                        default=["single_mutant", "truncation"],
                        format_func=lambda x: _strategy_labels[x],
                    )
                    variant_strategies = selected_strategies if selected_strategies else ["single_mutant"]
                else:
                    # Simple mode: 自動選択
                    variant_strategies = ["single_mutant", "truncation"]

                n_est = estimate_variant_count(
                    seed_sequence, variant_strategies,
                    avoid_residues=["C"] if avoid_cysteine else [],
                )
                st.caption(
                    f"シード: `{seed_sequence}`（{len(seed_sequence)} AA）  "
                    f"予測バリアント数: 約 {n_est} 件"
                )
            else:
                st.warning(f"配列エラー: {err_msg}")

        # ── Run button ───────────────────────────────
        _no_structure = uploaded_structure is None
        if mode == "Simple" and _no_structure:
            st.info("構造ファイルをアップロードまたは RCSB から検索してください。")
        run_button = st.button(
            "▶ Run" if mode == "Simple" else "Generate and Filter",
            type="primary",
            disabled=(_no_structure and mode == "Simple"),
        )

        # ── Structure summary ─────────────────────────
        if pdb_parse_error is not None:
            st.error(f"Failed to parse structure: {pdb_parse_error}")
        elif pdb_summary is not None:
            with st.expander(
                "Parsed structure summary",
                expanded=(mode == "Expert"),
            ):
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
                st.write(f"Residues: {pdb_summary['residue_count']}")
                st.write(f"Atoms: {pdb_summary['atom_count']}")
                st.write(f"Auto charge: {pdb_summary['pocket_charge_guess']}")
                st.write(f"Auto hydrophobicity: {pdb_summary['pocket_hydrophobicity_guess']}")
        else:
            st.caption("Load a structure to see the parsed summary here.")

    # Structure bytes for 3D viewer
    structure_bytes = uploaded_structure.getvalue() if uploaded_structure is not None else None
    structure_filename = uploaded_structure.name if uploaded_structure is not None else None

    return {
        "target_name": target_name if target_name else "peptide",
        "enable_docking": enable_docking,
        "docking_params": docking_params,
        "num_candidates": num_candidates,
        "min_len": min_len,
        "max_len": max_len,
        "pocket_charge": pocket_charge,
        "pocket_hydrophobicity": pocket_hydrophobicity,
        "avoid_cysteine": avoid_cysteine,
        "max_abs_charge": max_abs_charge,
        "max_hydrophobicity": max_hydrophobicity,
        "max_repeat_residue": max_repeat_residue,
        "preferred_len_min": preferred_len_min,
        "preferred_len_max": preferred_len_max,
        "use_diversity_filter": use_diversity_filter,
        "min_diversity_distance": min_diversity_distance,
        "max_diverse_candidates": max_diverse_candidates,
        "known_sequences": known_sequences,
        "pdb_summary": pdb_summary,
        "structure_bytes": structure_bytes,
        "structure_filename": structure_filename,
        "run_button": run_button,
        "selectivity_params": selectivity_params,
        "seed_sequence": seed_sequence,
        "variant_strategies": variant_strategies,
    }

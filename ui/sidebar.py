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
                        st.info("No results found on RCSB. / RCSBで候補が見つかりませんでした。検索語を変えて試してください。")
                except error.HTTPError as e:
                    st.error(f"RCSB search HTTP error / RCSB検索でHTTPエラー: {e}")
                except Exception as e:
                    st.error(f"RCSB search error / RCSB検索でエラー: {type(e).__name__}: {e}")

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
                        st.success(f"Loaded structure: {selected_record['pdb_id']}")
                    except error.HTTPError as e:
                        st.error(f"Structure download HTTP error / 構造ダウンロードでHTTPエラー: {e}")
                    except Exception as e:
                        st.error(f"Structure download error / 構造ダウンロードでエラー: {type(e).__name__}: {e}")

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
            help="Reference motifs for comparison. For explainability, not guaranteed known binders. / PoC向け参照モチーフ。厳密な既知 binder 保証ではなく説明性向上用です。",
        )

        custom_known_sequences_text = st.text_area(
            "Custom known sequences (one per line / comma / semicolon / space separated)",
            value="",
            height=120,
            help="例: KLAKLAKKLAKLAK",
        )

        st.markdown("### Docking (Phase B-1)")
        from core.docking import is_vina_available, compute_pocket_center
        vina_ok = is_vina_available()

        if not vina_ok:
            st.caption("Docking unavailable: vina binary not found (bin/vina) / vina バイナリが見つかりません（bin/vina）")
            enable_docking = False
            docking_params = None
        else:
            enable_docking = st.checkbox("Enable AutoDock Vina docking", value=False)
            if enable_docking:
                # ポケット中心を自動計算（構造が読み込まれている場合）
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
                    "Exhaustiveness",
                    options=[4, 8, 16, 32],
                    value=8,
                    help="Higher = more accurate but slower. / 高いほど精度が上がるが時間がかかる",
                )

                docking_params = {
                    "box_center": (cx, cy, cz),
                    "box_size": (float(box_sx), float(box_sy), float(box_sz)),
                    "top_n": docking_top_n,
                    "exhaustiveness": docking_exhaustiveness,
                }
            else:
                docking_params = None

        st.markdown("### Selectivity (Phase C-1)")
        enable_selectivity = st.checkbox(
            "Enable Selectivity Mode",
            value=False,
            help="Computes selectivity score = target fit − off-target fit. / ターゲットとオフターゲットのポケット適合差（選択性スコア）を計算します",
        )
        selectivity_params = None
        if enable_selectivity:
            offtarget_label = st.text_input("Off-target label", value="Off-target")
            st.caption(
                "Specify off-target pocket properties. "
                "Selectivity score = target rescoring − off-target rescoring  "
                "｜  オフターゲットのポケット特性を指定してください。"
                "選択性スコア = ターゲット rescoring − オフターゲット rescoring"
            )
            offtarget_charge = st.selectbox(
                "Off-target pocket charge",
                charge_options,
                index=0,
                key="offtarget_charge",
            )
            offtarget_hydrophobicity = st.selectbox(
                "Off-target pocket hydrophobicity",
                hydrophobicity_options,
                index=1,
                key="offtarget_hydrophobicity",
            )

            st.markdown("**Selectivity weight λ (ranking influence) / 選択性重み λ（ランキングへの反映度）** — Recommended / 推奨: 0.3")
            selectivity_lambda = st.slider(
                "λ",
                min_value=0.0,
                max_value=1.0,
                value=0.3,
                step=0.05,
                help=(
                    "Ranking score = final_score + λ × selectivity_score\n\n"
                    "λ = 0.0: rank by binding affinity only (off-target risk ignored).\n"
                    "Higher λ: promotes candidates with low off-target affinity.\n"
                    "λ = 1.0: binding affinity and selectivity weighted equally.\n"
                    "Recommended 0.3: balanced — influences ranking without large reversals.\n\n"
                    "ランキングスコア = final_score + λ × selectivity_score\n"
                    "λ = 0.0 のとき、ターゲット結合力のみでランキングします。\n"
                    "λ を大きくするほど、オフターゲットへの親和性が低い候補を上位に評価します。\n"
                    "推奨値 0.3: ランキングに影響は出るが final_score を大きく逆転させない適度な強さ。"
                ),
            )
            lambda_label = (
                "Binding affinity only (selectivity ignored) / ターゲット結合力のみ（選択性を考慮しない）"
                if selectivity_lambda == 0.0
                else f"Selectivity applied (λ={selectivity_lambda:.2f}) / 選択性を反映（λ={selectivity_lambda:.2f}）"
            )
            st.caption(f"Current setting / 現在の設定: {lambda_label}")

            selectivity_params = {
                "offtarget_label": offtarget_label,
                "pocket_charge_offtarget": offtarget_charge,
                "pocket_hydrophobicity_offtarget": offtarget_hydrophobicity,
                "selectivity_lambda": selectivity_lambda,
            }

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
            st.caption("Load a structure to see the parsed summary here. / 構造を選択すると parsed structure summary が表示されます。")

    # Assemble known_sequences from presets + custom input
    preset_sequences = [preset_dict[label] for label in selected_presets]
    custom_sequences = parse_known_sequences(custom_known_sequences_text)

    known_sequences = []
    seen = set()
    for seq in preset_sequences + custom_sequences:
        if seq not in seen:
            seen.add(seq)
            known_sequences.append(seq)

    # Structure bytes for 3D viewer
    structure_bytes = uploaded_structure.getvalue() if uploaded_structure is not None else None
    structure_filename = uploaded_structure.name if uploaded_structure is not None else None

    return {
        "target_name": target_name,
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
    }

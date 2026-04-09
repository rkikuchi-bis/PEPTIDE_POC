import streamlit as st


# カラーパレット
_COLOR_OTHER_CHAINS = "#bbbbbb"   # 非選択チェーン：薄いグレー
_COLOR_SELECTED_CHAIN = "#5b9bd5"  # 選択チェーン：ブルー
_COLOR_POCKET_STICK = "orangeCarbon"  # ポケット残基：オレンジ
_COLOR_LIGAND = "greenCarbon"         # リガンド：グリーン
_POCKET_SURFACE_OPACITY = 0.25


def _file_format(filename: str) -> str:
    lower = filename.lower()
    if lower.endswith(".pdb"):
        return "pdb"
    return "cif"


def _resi_range(pdb_summary: dict) -> tuple[int | None, int | None]:
    """ポケット残基の開始・終了番号を返す（なければ None）。"""
    start = pdb_summary.get("residue_start") or pdb_summary.get("used_residue_min")
    end = pdb_summary.get("residue_end") or pdb_summary.get("used_residue_max")
    return start, end


def render_structure_viewer(
    structure_bytes: bytes,
    structure_filename: str,
    pdb_summary: dict,
    width: int = 720,
    height: int = 460,
) -> None:
    """
    py3Dmol を使ってタンパク質構造を 3D 表示する。

    - 全チェーン：薄いグレーのカートゥーン
    - 選択チェーン：ブルーのカートゥーン
    - ポケット領域（manual_region）：オレンジのスティック＋半透明サーフェス
    - リガンド（ligand_neighborhood）：グリーンのスティック
    """
    try:
        import py3Dmol
    except ImportError:
        st.warning("py3Dmol がインストールされていません。`uv add py3Dmol` を実行してください。")
        return

    structure_text = structure_bytes.decode("utf-8", errors="ignore")
    fmt = _file_format(structure_filename)

    view = py3Dmol.view(width=width, height=height)
    view.addModel(structure_text, fmt)

    # ── ベーススタイル：全体をグレーのカートゥーン ──
    view.setStyle({}, {"cartoon": {"color": _COLOR_OTHER_CHAINS, "opacity": 0.6}})

    selected_chain = pdb_summary.get("selected_chain")
    source_mode = pdb_summary.get("source_mode", "manual_region")

    # ── 選択チェーン：ブルーのカートゥーン ──
    if selected_chain:
        view.setStyle(
            {"chain": selected_chain},
            {"cartoon": {"color": _COLOR_SELECTED_CHAIN}},
        )

    # ── ポケット強調 ──
    if source_mode == "manual_region" and selected_chain:
        res_start, res_end = _resi_range(pdb_summary)
        if res_start is not None and res_end is not None:
            resi_str = f"{int(res_start)}-{int(res_end)}"
            pocket_sel = {"chain": selected_chain, "resi": resi_str}

            # スティック表示
            view.addStyle(pocket_sel, {"stick": {"colorscheme": _COLOR_POCKET_STICK, "radius": 0.25}})

            # 半透明サーフェス
            view.addSurface(
                py3Dmol.VDW,
                {"opacity": _POCKET_SURFACE_OPACITY, "color": "orange"},
                pocket_sel,
            )

    elif source_mode == "ligand_neighborhood":
        ligand_names = pdb_summary.get("ligand_names", [])

        # リガンドをグリーンのスティックで表示
        for lig in ligand_names:
            view.addStyle({"resn": lig}, {"stick": {"colorscheme": _COLOR_LIGAND, "radius": 0.4}})

        # ポケット残基（おおよその範囲）をオレンジのスティックで表示
        res_start, res_end = _resi_range(pdb_summary)
        if selected_chain and res_start is not None and res_end is not None:
            resi_str = f"{int(res_start)}-{int(res_end)}"
            view.addStyle(
                {"chain": selected_chain, "resi": resi_str},
                {"stick": {"colorscheme": _COLOR_POCKET_STICK, "radius": 0.2, "opacity": 0.6}},
            )

    view.zoomTo()
    view.spin(False)

    html = view._make_html()
    st.iframe(html, height=height + 15)


def render_viewer_section(
    structure_bytes: bytes | None,
    structure_filename: str | None,
    pdb_summary: dict | None,
) -> None:
    """
    app.py から呼ぶエントリポイント。
    構造が読み込まれているときだけビューアを表示する。
    """
    if structure_bytes is None or pdb_summary is None:
        return

    with st.expander("3D Structure Viewer", expanded=True):
        chain = pdb_summary.get("selected_chain", "-")
        mode = pdb_summary.get("source_mode", "-")
        res_start, res_end = _resi_range(pdb_summary)
        ligands = pdb_summary.get("ligand_names", [])

        # 凡例
        legend_parts = [
            "🔵 選択チェーン",
            "🟠 ポケット領域",
        ]
        if ligands:
            legend_parts.append("🟢 リガンド")
        st.caption(
            f"Chain: **{chain}** | Mode: **{mode}**"
            + (f" | Residues: {res_start}–{res_end}" if res_start and res_end else "")
            + (f" | Ligand: {', '.join(ligands)}" if ligands else "")
            + f"　　{'  '.join(legend_parts)}"
        )

        render_structure_viewer(structure_bytes, structure_filename, pdb_summary)

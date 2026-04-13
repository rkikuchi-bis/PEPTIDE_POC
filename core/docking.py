"""
Phase B-1: AutoDock Vina バイナリによるペプチドドッキングスコア計算

フロー:
  1. ペプチド配列 → RDKit 3D構造 → meeko PDBQT (リガンド)
  2. 受容体PDB → meeko PDBQT (受容体、キャッシュ付き)
  3. vina バイナリ実行 (subprocess)
  4. スコア [kcal/mol] を返却

使用方法:
    from core.docking import dock_top_candidates
    result_df = dock_top_candidates(result_df, pdb_path, box_center, box_size, top_n=10)
"""
from __future__ import annotations

import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import pandas as pd

# vina バイナリのデフォルトパス（プロジェクトルートの bin/ 以下）
_PROJECT_ROOT = Path(__file__).parent.parent
_VINA_BINARY = _PROJECT_ROOT / "bin" / "vina"


def compute_pocket_center(
    structure_bytes: bytes,
    filename: str,
    pdb_summary: dict,
) -> Optional[tuple[float, float, float]]:
    """
    pdb_summary の情報を元に構造ファイルからポケット中心座標を計算する。

    - Manual region モード: 指定チェーン・残基範囲の重心
    - Ligand モード: リガンド原子の重心

    Returns
    -------
    (x, y, z) or None
    """
    try:
        from Bio.PDB import PDBParser, MMCIFParser
        import io
        import numpy as np

        # パーサー選択
        fmt = pdb_summary.get("file_format", "")
        if "cif" in fmt.lower() or filename.lower().endswith((".cif", ".mmcif")):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        with tempfile.NamedTemporaryFile(
            suffix="." + filename.rsplit(".", 1)[-1], delete=False
        ) as tmp:
            tmp.write(structure_bytes)
            tmp_path = tmp.name

        try:
            structure = parser.get_structure("s", tmp_path)
        finally:
            os.unlink(tmp_path)

        chain_id = pdb_summary.get("selected_chain")
        res_start = pdb_summary.get("residue_start")
        res_end = pdb_summary.get("residue_end")
        ligand_names = pdb_summary.get("ligand_names", [])
        source_mode = pdb_summary.get("source_mode", "")

        coords = []
        model = next(structure.get_models())

        if "ligand" in source_mode.lower() and ligand_names:
            # リガンド原子の重心
            for chain in model.get_chains():
                if chain_id and chain.get_id() != chain_id:
                    continue
                for residue in chain.get_residues():
                    if residue.get_resname().strip() in [ln.strip() for ln in ligand_names]:
                        for atom in residue.get_atoms():
                            coords.append(atom.get_vector().get_array())
        else:
            # 指定チェーン・残基範囲の重心
            for chain in model.get_chains():
                if chain_id and chain.get_id() != chain_id:
                    continue
                for residue in chain.get_residues():
                    seq_id = residue.get_id()[1]
                    if res_start is not None and seq_id < res_start:
                        continue
                    if res_end is not None and seq_id > res_end:
                        continue
                    for atom in residue.get_atoms():
                        coords.append(atom.get_vector().get_array())

        if not coords:
            # フォールバック：全原子の重心
            for atom in model.get_atoms():
                coords.append(atom.get_vector().get_array())

        if not coords:
            return None

        center = np.mean(coords, axis=0)
        return (round(float(center[0]), 2), round(float(center[1]), 2), round(float(center[2]), 2))

    except Exception:
        return None


def is_vina_available() -> bool:
    """vina バイナリが実行可能かチェック。"""
    if not _VINA_BINARY.exists():
        return False
    try:
        result = subprocess.run(
            [str(_VINA_BINARY), "--version"],
            capture_output=True, text=True, timeout=5,
        )
        return result.returncode == 0
    except Exception:
        return False


# ─────────────────────────────────────────────
# リガンド（ペプチド）PDBQT 準備
# ─────────────────────────────────────────────

def _sequence_to_pdbqt(sequence: str, out_path: str, flexible: bool = False) -> bool:
    """
    ペプチド配列 → RDKit 3D構造 → obabel PDBQT ファイル。

    flexible=False (デフォルト):
        obabel -xr で全結合を剛体 PDBQT 化。torsion 数を抑えて高速処理。
    flexible=True:
        obabel の標準変換（-xr なし）でtorsion を保持した柔軟 PDBQT を生成。
        短鎖ペプチド（≤5残基）に適用することでスコアの精度が向上する。

    成功時 True、失敗時 False を返す。
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        # ペプチド配列 → RDKit 分子
        mol = Chem.MolFromSequence(sequence)
        if mol is None:
            return False

        # 3D 座標を生成（複数の設定でフォールバック）
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        result = AllChem.EmbedMolecule(mol, params)

        if result != 0:
            params2 = AllChem.ETKDGv3()
            params2.randomSeed = 0
            params2.useRandomCoords = True
            result = AllChem.EmbedMolecule(mol, params2)

        if result != 0:
            return False

        AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)

        with tempfile.TemporaryDirectory() as tmpdir:
            sdf_path = os.path.join(tmpdir, "lig.sdf")
            writer = Chem.SDWriter(sdf_path)
            writer.write(mol)
            writer.close()

            if flexible:
                # 柔軟ドッキング: obabel の標準変換で torsion を保持した PDBQT を生成
                r = subprocess.run(
                    ["obabel", sdf_path, "-O", out_path,
                     "--partialcharge", "gasteiger"],
                    capture_output=True, text=True, timeout=30,
                )
                return r.returncode == 0 and os.path.exists(out_path)
            else:
                # 剛体ドッキング: -xr で ATOM 行のみ出力 → ROOT/ENDROOT/TORSDOF 0 を手動付加
                rigid_path = os.path.join(tmpdir, "lig_rigid.pdbqt")
                r = subprocess.run(
                    ["obabel", sdf_path, "-O", rigid_path, "-xr",
                     "--partialcharge", "gasteiger"],
                    capture_output=True, text=True, timeout=30,
                )
                if r.returncode != 0 or not os.path.exists(rigid_path):
                    return False

                with open(rigid_path) as f:
                    lines = f.readlines()
                remark_lines = [l for l in lines if l.startswith("REMARK")]
                atom_lines = [l for l in lines if l.startswith("ATOM") or l.startswith("HETATM")]
                fixed = remark_lines + ["ROOT\n"] + atom_lines + ["ENDROOT\n", "TORSDOF 0\n"]
                with open(out_path, "w") as f:
                    f.writelines(fixed)
                return True

    except Exception:
        return False


# 柔軟ドッキングを適用する最大ペプチド長（この残基数以下が対象）
FLEXIBLE_DOCKING_MAX_LENGTH = 5


# ─────────────────────────────────────────────
# 受容体 PDBQT 準備
# ─────────────────────────────────────────────

def _prepare_receptor_pdbqt(pdb_path: str, out_path: str) -> bool:
    """
    受容体 PDB / mmCIF → obabel 剛体 PDBQT ファイル。

    mmCIF (.cif) の場合は BioPython で PDB 形式に変換してから obabel に渡す。
    PDB の場合は ATOM 行のみ抽出して obabel に渡す。
    成功時 True、失敗時 False を返す。
    """
    try:
        is_cif = pdb_path.lower().endswith((".cif", ".mmcif"))

        if is_cif:
            # mmCIF → BioPython で PDB 文字列に変換（固定幅 PDB 形式）
            import io
            from Bio.PDB import MMCIFParser, PDBIO, Select

            class ProteinOnly(Select):
                def accept_residue(self, residue):
                    return residue.get_id()[0] == " "  # HETATM を除外

            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("receptor", pdb_path)
            pdbio = PDBIO()
            pdbio.set_structure(structure)
            buf = io.StringIO()
            pdbio.save(buf, select=ProteinOnly())
            pdb_text = buf.getvalue()
            lines = pdb_text.splitlines(keepends=True)
        else:
            with open(pdb_path) as f:
                lines = f.readlines()

        # ATOM / TER / END のみ残す（mmCIF 変換後は既に HETATM 除去済みだが念のため）
        protein_lines = [
            l for l in lines
            if l.startswith("ATOM") or l.startswith("TER") or l.startswith("END")
        ]
        with tempfile.NamedTemporaryFile(
            suffix=".pdb", mode="w", delete=False
        ) as tmp:
            tmp.writelines(protein_lines)
            tmp_path = tmp.name

        try:
            result = subprocess.run(
                ["obabel", tmp_path, "-O", out_path, "-xr",
                 "--partialcharge", "gasteiger"],
                capture_output=True, text=True, timeout=60,
            )
            return result.returncode == 0 and os.path.exists(out_path)
        finally:
            os.unlink(tmp_path)

    except Exception:
        return False


# ─────────────────────────────────────────────
# ドッキング実行
# ─────────────────────────────────────────────

def _run_vina(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    out_pdbqt: str,
    center: tuple[float, float, float],
    size: tuple[float, float, float],
    exhaustiveness: int = 8,
    num_modes: int = 3,
    timeout: int = 120,
) -> Optional[float]:
    """
    vina バイナリを subprocess で実行し、最良スコア [kcal/mol] を返す。
    失敗時は None。
    """
    cmd = [
        str(_VINA_BINARY),
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--out", out_pdbqt,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
    ]
    try:
        result = subprocess.run(
            cmd,
            capture_output=True, text=True, timeout=timeout,
        )
        if result.returncode != 0:
            return None

        # vina 出力の affinity 表をパース
        # 例: "   1         -7.2      0.000      0.000"  (負が正常、正はクラッシュ)
        for line in result.stdout.splitlines():
            m = re.match(r"\s+1\s+([-+]?\d+\.?\d*)", line)
            if m:
                return float(m.group(1))
        return None

    except subprocess.TimeoutExpired:
        return None
    except Exception:
        return None


# ─────────────────────────────────────────────
# 公開 API
# ─────────────────────────────────────────────

def dock_peptide(
    sequence: str,
    receptor_pdbqt_path: str,
    box_center: tuple[float, float, float],
    box_size: tuple[float, float, float] = (20.0, 20.0, 20.0),
    exhaustiveness: int = 8,
    flexible: Optional[bool] = None,
) -> Optional[float]:
    """
    単一ペプチド配列をドッキングしてスコア [kcal/mol] を返す。

    Parameters
    ----------
    sequence : str
        1文字アミノ酸配列
    receptor_pdbqt_path : str
        受容体 PDBQT ファイルパス（事前に prepare_receptor_pdbqt() で準備）
    box_center : tuple[float, float, float]
        ポケット中心座標 (x, y, z) [Å]
    box_size : tuple[float, float, float]
        サーチボックスサイズ (x, y, z) [Å]（デフォルト 20x20x20）
    exhaustiveness : int
        探索の徹底度（デフォルト 8、速度優先なら 4）
    flexible : bool or None
        True=柔軟、False=剛体、None=配列長で自動判定
        (≤FLEXIBLE_DOCKING_MAX_LENGTH 残基 → 柔軟、それ以上 → 剛体)

    Returns
    -------
    float or None
        ドッキングスコア [kcal/mol]（小さいほど良い）、失敗時は None
    """
    if not is_vina_available():
        return None

    # 自動判定: 短鎖は柔軟ドッキング、長鎖は剛体
    use_flexible = flexible if flexible is not None else (len(sequence) <= FLEXIBLE_DOCKING_MAX_LENGTH)

    with tempfile.TemporaryDirectory() as tmpdir:
        ligand_pdbqt = os.path.join(tmpdir, "ligand.pdbqt")
        out_pdbqt = os.path.join(tmpdir, "out.pdbqt")

        if not _sequence_to_pdbqt(sequence, ligand_pdbqt, flexible=use_flexible):
            return None

        # 柔軟ドッキングはサーチ空間が広いため timeout を延長
        timeout = 180 if use_flexible else 120
        return _run_vina(
            receptor_pdbqt=receptor_pdbqt_path,
            ligand_pdbqt=ligand_pdbqt,
            out_pdbqt=out_pdbqt,
            center=box_center,
            size=box_size,
            exhaustiveness=exhaustiveness,
            timeout=timeout,
        )


def prepare_receptor_pdbqt(pdb_path: str, cache_dir: Optional[str] = None) -> Optional[str]:
    """
    受容体 PDB ファイルを PDBQT に変換して返す（パス文字列）。

    Parameters
    ----------
    pdb_path : str
        入力 PDB ファイルパス
    cache_dir : str, optional
        PDBQT キャッシュ保存先ディレクトリ（省略時は inputs/ 直下）

    Returns
    -------
    str or None
        生成された PDBQT ファイルパス、失敗時は None
    """
    pdb_path = str(pdb_path)
    stem = Path(pdb_path).stem

    if cache_dir is None:
        cache_dir = str(_PROJECT_ROOT / "inputs")
    os.makedirs(cache_dir, exist_ok=True)

    out_path = os.path.join(cache_dir, f"{stem}_receptor.pdbqt")

    # キャッシュが存在すればそのまま返す
    if os.path.exists(out_path):
        return out_path

    success = _prepare_receptor_pdbqt(pdb_path, out_path)
    return out_path if success else None


def dock_top_candidates(
    result_df: pd.DataFrame,
    receptor_pdbqt_path: str,
    box_center: tuple[float, float, float],
    box_size: tuple[float, float, float] = (20.0, 20.0, 20.0),
    top_n: int = 15,
    exhaustiveness: int = 8,
    flexible: Optional[bool] = None,
) -> pd.DataFrame:
    """
    result_df の上位 top_n 候補に対してドッキングを実行し、
    docking_score カラムを追加して返す。

    ドッキングを実行しない候補（top_n 以降）は NaN。
    vina が利用不可の場合はカラムを追加せずそのまま返す。

    Parameters
    ----------
    result_df : pd.DataFrame
        run_pipeline() の出力 DataFrame（rank カラムを含む）
    receptor_pdbqt_path : str
        受容体 PDBQT ファイルパス
    box_center, box_size : tuple
        ドッキングボックス設定
    top_n : int
        ドッキングを実行する上位候補数（デフォルト 15）
    exhaustiveness : int
        vina の探索徹底度
    flexible : bool or None
        True=柔軟、False=剛体、None=配列長で自動判定（各配列ごとに適用）
    """
    if not is_vina_available():
        return result_df

    out = result_df.copy()
    out["docking_score"] = float("nan")
    out["docking_mode"] = ""

    for idx in out.index:
        if "rank" in out.columns and out.at[idx, "rank"] > top_n:
            continue

        seq = out.at[idx, "sequence"]
        use_flexible = flexible if flexible is not None else (len(seq) <= FLEXIBLE_DOCKING_MAX_LENGTH)

        score = dock_peptide(
            sequence=seq,
            receptor_pdbqt_path=receptor_pdbqt_path,
            box_center=box_center,
            box_size=box_size,
            exhaustiveness=exhaustiveness,
            flexible=use_flexible,
        )
        if score is not None:
            out.at[idx, "docking_score"] = score
        out.at[idx, "docking_mode"] = "flexible" if use_flexible else "rigid"

    return out

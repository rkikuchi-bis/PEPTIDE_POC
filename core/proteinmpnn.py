"""
Phase B-2: ProteinMPNN による配列スコアリング

LightGBM との OpenMP 競合を回避するため、ProteinMPNN の推論は
scripts/mpnn_scorer.py をサブプロセスとして起動して実行する。

複数配列を一括スコアリングする場合は score_sequences_batch() を使うこと
（サブプロセス起動コストを一度にまとめられる）。

受容体構造を条件付けとしてスコアリングする場合は
score_sequences_with_receptor() を使うこと（Phase B-2+/B-2++）。

Streamlit に依存しない純粋な Python モジュール。
"""
from __future__ import annotations

import pickle
import subprocess
import sys
import tempfile
import os
from pathlib import Path
from typing import Sequence

# ProteinMPNN ウェイトのパス（存在チェックに使用）
_ROOT = Path(__file__).resolve().parent.parent
_WEIGHTS_PATH = _ROOT / "models" / "proteinmpnn" / "v_48_020.pt"
_SCORER_SCRIPT = _ROOT / "scripts" / "mpnn_scorer.py"
_SCORER_RECEPTOR_SCRIPT = _ROOT / "scripts" / "mpnn_scorer_receptor.py"


def is_proteinmpnn_available() -> bool:
    """ProteinMPNN ウェイトが存在するか確認する。"""
    return _WEIGHTS_PATH.exists() and _SCORER_SCRIPT.exists()


def score_sequences_batch(sequences: Sequence[str]) -> list[float]:
    """
    複数ペプチド配列の ProteinMPNN スコアを一括計算する。

    サブプロセスを一度だけ起動して全配列をまとめて処理するため、
    score_with_proteinmpnn() を繰り返し呼ぶより高速。

    Args:
        sequences: ペプチド配列のリスト

    Returns:
        各配列に対応するスコア [0, 1] のリスト。
        エラー時は 0.5（中立）を返す。
    """
    if not sequences:
        return []

    if not is_proteinmpnn_available():
        return [0.5] * len(sequences)

    try:
        result = subprocess.run(
            [sys.executable, str(_SCORER_SCRIPT)] + list(sequences),
            capture_output=True,
            text=True,
            timeout=300,  # 最大5分
        )
        lines = result.stdout.strip().split("\n")
        scores = []
        for i, line in enumerate(lines):
            try:
                scores.append(float(line.strip()))
            except ValueError:
                scores.append(0.5)
        # 配列数と出力行数が合わない場合は残りを 0.5 で埋める
        while len(scores) < len(sequences):
            scores.append(0.5)
        return scores[:len(sequences)]
    except Exception:
        return [0.5] * len(sequences)


def score_with_proteinmpnn(sequence: str) -> float:
    """
    1配列の ProteinMPNN スコアを計算する。

    複数配列をまとめて処理する場合は score_sequences_batch() を使うこと。

    Returns:
        0.0〜1.0 のスコア（高いほど設計可能性が高い）。
        エラー時は 0.5（中立）を返す。
    """
    scores = score_sequences_batch([sequence])
    return scores[0] if scores else 0.5


def score_sequences_with_receptor(
    sequences: Sequence[str],
    structure_text: str,
    file_format: str,
    centroid: tuple[float, float, float],
) -> tuple[list[float], list[bool]]:
    """
    受容体構造を条件付けとして、複数ペプチド配列を一括スコアリングする（Phase B-2+/B-2++）。

    受容体バックボーン (chain_M=0) + ペプチド骨格 (chain_M=1) を結合し、
    ペプチド残基の対数尤度のみをスコアとして返す。

    ペプチド骨格は ESMFold API で予測した座標を使用する（Phase B-2++）。
    API 呼び出しが失敗した配列は理想 αヘリックスで代替する（Phase B-2+ 互換）。

    Args:
        sequences: ペプチド配列のリスト
        structure_text: 受容体構造のテキスト（PDB または mmCIF 形式）
        file_format: "pdb" または "cif"
        centroid: ポケット重心座標 (x, y, z)

    Returns:
        (scores, esmfold_used_flags) のタプル
        - scores: 各配列に対応するスコア [0, 1] のリスト
        - esmfold_used_flags: 各配列で ESMFold 骨格を使用したか否かのフラグリスト
        エラー時は scores=0.5（中立）、flags=False を返す。
    """
    if not sequences:
        return [], []

    if not _WEIGHTS_PATH.exists() or not _SCORER_RECEPTOR_SCRIPT.exists():
        return [0.5] * len(sequences), [False] * len(sequences)

    # ESMFold による骨格予測（Phase B-2++）
    try:
        from core.pepfold import predict_backbones_batch
        peptide_coords_list = predict_backbones_batch(list(sequences))
    except Exception:
        peptide_coords_list = [None] * len(sequences)

    esmfold_used = [c is not None for c in peptide_coords_list]

    # 受容体 PDB を一時ファイルに書き出す
    suffix = ".pdb" if file_format == "pdb" else ".cif"
    tmp_receptor = None
    tmp_coords = None
    try:
        with tempfile.NamedTemporaryFile(
            suffix=suffix, delete=False, mode="w", encoding="utf-8"
        ) as tmp:
            tmp.write(structure_text)
            tmp_receptor = tmp.name

        # mmCIF の場合は BioPython で PDB に変換して渡す
        if file_format == "cif":
            tmp_receptor = _convert_cif_to_pdb_temp(structure_text)
            if tmp_receptor is None:
                return [0.5] * len(sequences), [False] * len(sequences)

        # ESMFold 骨格座標を pickle で一時ファイルに保存してサブプロセスに渡す
        with tempfile.NamedTemporaryFile(
            suffix=".pkl", delete=False, mode="wb"
        ) as tmp:
            pickle.dump(peptide_coords_list, tmp)
            tmp_coords = tmp.name

        cx, cy, cz = centroid
        result = subprocess.run(
            [
                sys.executable, str(_SCORER_RECEPTOR_SCRIPT),
                "--receptor", tmp_receptor,
                "--centroid", str(cx), str(cy), str(cz),
                "--peptide-coords", tmp_coords,
            ] + list(sequences),
            capture_output=True,
            text=True,
            timeout=600,  # 最大10分（受容体解析込み）
        )
        lines = result.stdout.strip().split("\n")
        scores = []
        for line in lines:
            try:
                scores.append(float(line.strip()))
            except ValueError:
                scores.append(0.5)
        while len(scores) < len(sequences):
            scores.append(0.5)
        return scores[:len(sequences)], esmfold_used

    except Exception:
        return [0.5] * len(sequences), [False] * len(sequences)
    finally:
        for path in (tmp_receptor, tmp_coords):
            if path and os.path.exists(path):
                os.unlink(path)


def _convert_cif_to_pdb_temp(cif_text: str) -> str | None:
    """mmCIF テキストを PDB 形式の一時ファイルに変換して返す。失敗時は None。"""
    try:
        from io import StringIO
        from Bio.PDB import MMCIFParser, PDBIO

        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("rec", StringIO(cif_text))

        pdb_io = PDBIO()
        pdb_io.set_structure(structure)

        with tempfile.NamedTemporaryFile(
            suffix=".pdb", delete=False, mode="w", encoding="utf-8"
        ) as tmp:
            pdb_io.save(tmp)
            return tmp.name
    except Exception:
        return None

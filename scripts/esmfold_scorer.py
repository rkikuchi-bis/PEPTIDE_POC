"""
ESMFold ローカル推論スクリプト（Phase B-2++）

core/pepfold.py からサブプロセスとして呼び出される。
LightGBM との OpenMP 競合を避けるため、このプロセスには LightGBM は存在しない。

使い方:
    python scripts/esmfold_scorer.py \
        --output /tmp/coords.pkl \
        ACDEFGHIK WFYWFYWFY ...

出力:
    --output に list[np.ndarray | None] を pickle で保存する。
    各要素は (L, 4, 3) の骨格座標（N, CA, C, O）。予測失敗時は None。

初回実行時:
    fair-esm が ESMFold モデル重み（約 700MB）を
    ~/.cache/torch/hub/ に自動ダウンロードする。
"""
import sys
import os
import pickle
import argparse

import numpy as np
import torch

# OpenMP 競合を念のため抑制
os.environ.setdefault("OMP_NUM_THREADS", "1")

AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def _parse_pdb_backbone(pdb_text: str) -> np.ndarray | None:
    """PDB テキストから N, CA, C, O 座標を抽出して (L, 4, 3) 配列を返す。"""
    residues: dict[tuple, dict[str, np.ndarray]] = {}
    residue_order: list[tuple] = []

    for line in pdb_text.splitlines():
        record = line[:6].strip()
        if record not in ("ATOM", "HETATM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name not in ("N", "CA", "C", "O"):
            continue
        resname = line[17:20].strip()
        if resname not in AA3_TO_AA1:
            continue
        chain = line[21:22].strip()
        try:
            res_seq = int(line[22:26].strip())
        except ValueError:
            continue

        key = (chain, res_seq, resname)
        if key not in residues:
            residues[key] = {}
            residue_order.append(key)

        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            residues[key][atom_name] = np.array([x, y, z], dtype=np.float32)
        except ValueError:
            continue

    if not residue_order:
        return None

    coords = np.zeros((len(residue_order), 4, 3), dtype=np.float32)
    for i, key in enumerate(residue_order):
        atoms = residues[key]
        ca = atoms.get("CA", np.zeros(3, dtype=np.float32))
        coords[i, 0] = atoms.get("N", ca)
        coords[i, 1] = ca
        coords[i, 2] = atoms.get("C", ca)
        coords[i, 3] = atoms.get("O", ca)

    return coords


def _load_model():
    """ESMFold モデルをロードする。MPS → CPU の順でフォールバック。"""
    import esm  # type: ignore
    model = esm.pretrained.esmfold_v1()
    model = model.eval()

    if torch.backends.mps.is_available():
        try:
            model = model.to("mps")
            return model, torch.device("mps")
        except Exception:
            pass

    return model.to("cpu"), torch.device("cpu")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, help="出力 pickle ファイルパス")
    parser.add_argument("sequences", nargs="*", help="予測するペプチド配列")
    args = parser.parse_args()

    sequences = args.sequences
    if not sequences:
        with open(args.output, "wb") as f:
            pickle.dump([], f)
        sys.exit(0)

    try:
        model, device = _load_model()
    except Exception as e:
        print(f"[esmfold_scorer] モデルロード失敗: {e}", file=sys.stderr)
        with open(args.output, "wb") as f:
            pickle.dump([None] * len(sequences), f)
        sys.exit(0)

    results = []
    for seq in sequences:
        try:
            with torch.no_grad():
                pdb_text = model.infer_pdb(seq)
            coords = _parse_pdb_backbone(pdb_text)
            results.append(coords)
        except Exception as e:
            print(f"[esmfold_scorer] {seq[:10]}... 予測失敗: {e}", file=sys.stderr)
            results.append(None)

    with open(args.output, "wb") as f:
        pickle.dump(results, f)

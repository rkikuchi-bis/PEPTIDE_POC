"""
受容体構造を条件付けとした ProteinMPNN スコアリング（Phase B-2+）

core/proteinmpnn.py からサブプロセスとして呼び出される。
LightGBM との OpenMP 競合を避けるため、このプロセスには LightGBM は存在しない。

使い方:
    python scripts/mpnn_scorer_receptor.py \
        --receptor /tmp/receptor.pdb \
        --centroid 10.5 20.3 -5.1 \
        ACDEFGHIK WFYWFYWFY ...
    -> 1行1スコアで標準出力

設計:
    受容体バックボーン (chain_M=0) + ペプチド理想ヘリックス (chain_M=1) を結合し、
    ペプチド残基の対数尤度のみをスコアとして返す。
    ペプチドは pocket centroid を中心に配置（グラフ上で受容体残基と隣接させる）。
"""
import sys
import os
import math
import argparse

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_MPNN_DIR = os.path.join(_ROOT, "models", "proteinmpnn")
if _MPNN_DIR not in sys.path:
    sys.path.insert(0, _MPNN_DIR)

import numpy as np
import torch
from Bio.PDB import PDBParser

_WEIGHTS_PATH = os.path.join(_MPNN_DIR, "v_48_020.pt")
_ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
_AA_TO_IDX = {aa: i for i, aa in enumerate(_ALPHABET)}
_STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# 受容体が大きすぎる場合は centroid 周辺のみ使用する（高速化のため）
_MAX_RECEPTOR_RESIDUES = 200


def _extract_receptor_backbone(
    pdb_path: str, centroid: np.ndarray
) -> tuple[np.ndarray, str]:
    """
    PDB から受容体バックボーン座標 (N, CA, C, O) を抽出する。

    受容体残基が多すぎる場合は centroid に近い上位 _MAX_RECEPTOR_RESIDUES 残基のみ使用。

    Returns:
        coords: (L, 4, 3) - N, CA, C, O 座標（欠損原子は CA 座標で代替）
        sequence: str
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("rec", pdb_path)

    residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":  # HETATM skip
                    continue
                resname = residue.get_resname().strip()
                if resname not in AA3_TO_AA1:
                    continue
                if "CA" not in residue:
                    continue
                residues.append(residue)
        break  # first model only

    if not residues:
        return np.zeros((0, 4, 3), dtype=np.float32), ""

    # centroid に近い残基を優先（受容体コンテキストの効率化）
    if len(residues) > _MAX_RECEPTOR_RESIDUES:
        ca_coords = np.array([r["CA"].coord for r in residues])
        dists = np.linalg.norm(ca_coords - centroid, axis=1)
        idx = np.argsort(dists)[:_MAX_RECEPTOR_RESIDUES]
        # 元の順序を保持してから取り出す
        idx_sorted = np.sort(idx)
        residues = [residues[i] for i in idx_sorted]

    coords = np.zeros((len(residues), 4, 3), dtype=np.float32)
    seq = ""

    for i, residue in enumerate(residues):
        ca = residue["CA"].coord
        n = residue["N"].coord if "N" in residue else ca
        c = residue["C"].coord if "C" in residue else ca
        o = residue["O"].coord if "O" in residue else ca

        coords[i, 0] = n
        coords[i, 1] = ca
        coords[i, 2] = c
        coords[i, 3] = o
        seq += AA3_TO_AA1[residue.get_resname().strip()]

    return coords, seq


def _ideal_helix_at_centroid(length: int, centroid: np.ndarray) -> np.ndarray:
    """
    pocket centroid を中心に配置した理想 αヘリックス座標を生成する。

    ヘリックス軸は z 方向。centroid にヘリックス中央を合わせる。
    """
    RISE = 1.5
    RADIUS = 2.3
    DELTA = math.radians(100.0)

    coords = np.zeros((length, 4, 3), dtype=np.float32)
    z_center = centroid[2]
    z_offset = z_center - (length - 1) * RISE / 2.0

    for i in range(length):
        angle = i * DELTA
        ca = np.array([
            centroid[0] + RADIUS * math.cos(angle),
            centroid[1] + RADIUS * math.sin(angle),
            z_offset + i * RISE,
        ], dtype=np.float32)

        tangent = np.array([
            -math.sin(angle), math.cos(angle),
            RISE / math.sqrt(RADIUS ** 2 + RISE ** 2),
        ], dtype=np.float32)
        tangent /= np.linalg.norm(tangent)
        radial = np.array([math.cos(angle), math.sin(angle), 0.0], dtype=np.float32)

        coords[i, 0] = ca - tangent * 1.46   # N
        coords[i, 1] = ca                     # CA
        coords[i, 2] = ca + tangent * 1.52   # C
        coords[i, 3] = (ca + tangent * 1.52) + radial * 1.24  # O

    return coords


def _nll_to_score(nll: float, scale: float = 3.3) -> float:
    """
    NLL をシグモイド変換して [0, 1] スコアに変換する。

    scale = スコア 0.5 の中点となる NLL 値。
    受容体条件付きモードの実測 NLL 平均が約 3.3 のため、
    この値を中点として設定している（構造フリー版の 2.0 より高い）。
    """
    return 1.0 / (1.0 + math.exp(nll - scale))


def load_model(device):
    from protein_mpnn_utils import ProteinMPNN  # type: ignore
    ckpt = torch.load(_WEIGHTS_PATH, map_location="cpu", weights_only=False)
    model = ProteinMPNN(
        num_letters=21, node_features=128, edge_features=128,
        hidden_dim=128, num_encoder_layers=3, num_decoder_layers=3,
        vocab=21, k_neighbors=ckpt["num_edges"],
        augment_eps=0.0, dropout=0.0, ca_only=False,
    )
    model.load_state_dict(ckpt["model_state_dict"])
    model.to(device)
    model.eval()
    return model


@torch.no_grad()
def score_sequence_with_receptor(
    model,
    receptor_coords: np.ndarray,
    receptor_seq: str,
    peptide_sequence: str,
    centroid: np.ndarray,
    device,
) -> float:
    """
    受容体座標をコンテキストとして、ペプチド配列をスコアリングする。

    - 受容体残基: chain_M = 0（固定コンテキスト、スコアリング対象外）
    - ペプチド残基: chain_M = 1（デザイン対象、スコアリング対象）

    NLL はペプチド残基のみで計算する。
    """
    clean_pep = "".join(aa for aa in peptide_sequence.upper() if aa in _STANDARD_AA)
    if len(clean_pep) < 3:
        return 0.5

    L_rec = len(receptor_coords)
    L_pep = len(clean_pep)
    L_total = L_rec + L_pep

    # ペプチド座標を centroid に配置
    pep_coords = _ideal_helix_at_centroid(L_pep, centroid)

    # 座標結合: (L_total, 4, 3)
    if L_rec > 0:
        all_coords = np.concatenate([receptor_coords, pep_coords], axis=0)
    else:
        all_coords = pep_coords

    X = torch.tensor(all_coords, device=device).unsqueeze(0)  # (1, L, 4, 3)

    # 配列インデックス: 受容体 + ペプチド
    rec_idx = [_AA_TO_IDX.get(aa, _AA_TO_IDX["X"]) for aa in receptor_seq]
    pep_idx = [_AA_TO_IDX.get(aa, _AA_TO_IDX["X"]) for aa in clean_pep]
    S = torch.tensor([rec_idx + pep_idx], dtype=torch.long, device=device)

    # マスク: すべて有効
    mask = torch.ones(1, L_total, device=device)

    # chain_M: 受容体=0（固定）、ペプチド=1（スコアリング対象）
    chain_M = torch.zeros(1, L_total, device=device)
    if L_rec > 0:
        chain_M[0, L_rec:] = 1.0
    else:
        chain_M[0, :] = 1.0

    # residue_idx: sequential
    residue_idx = torch.arange(L_total, device=device).unsqueeze(0)

    # chain_encoding: 受容体=1, ペプチド=2
    chain_encoding = torch.ones(1, L_total, dtype=torch.long, device=device)
    chain_encoding[0, L_rec:] = 2

    randn = torch.zeros(1, L_total, device=device)

    log_probs = model(X, S, mask, chain_M, residue_idx, chain_encoding, randn)

    # ペプチド部分のみの NLL を計算
    pep_start = L_rec
    pep_log_probs = log_probs[0, pep_start:pep_start + L_pep, :]  # (L_pep, 21)
    pep_S = S[0, pep_start:pep_start + L_pep]  # (L_pep,)
    nll = float(-pep_log_probs[torch.arange(L_pep), pep_S].mean().item())

    return round(_nll_to_score(nll), 4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptor", required=True, help="受容体 PDB ファイルパス")
    parser.add_argument("--centroid", type=float, nargs=3, required=True,
                        metavar=("X", "Y", "Z"), help="ポケット重心座標")
    parser.add_argument("sequences", nargs="*", help="スコアリングするペプチド配列")
    args = parser.parse_args()

    sequences = args.sequences
    if not sequences:
        sys.exit(0)

    centroid = np.array(args.centroid, dtype=np.float32)
    device = torch.device("mps") if torch.backends.mps.is_available() else torch.device("cpu")

    # 受容体バックボーン抽出
    try:
        receptor_coords, receptor_seq = _extract_receptor_backbone(args.receptor, centroid)
    except Exception:
        receptor_coords = np.zeros((0, 4, 3), dtype=np.float32)
        receptor_seq = ""

    # モデルロード
    try:
        model = load_model(device)
    except Exception:
        for _ in sequences:
            print(0.5)
        sys.exit(0)

    for seq in sequences:
        try:
            score = score_sequence_with_receptor(
                model, receptor_coords, receptor_seq, seq, centroid, device
            )
        except Exception:
            score = 0.5
        print(score)

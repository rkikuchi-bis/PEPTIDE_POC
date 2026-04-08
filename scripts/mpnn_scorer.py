"""
ProteinMPNN スタンドアロンスコアラー

core/proteinmpnn.py からサブプロセスとして呼び出される。
LightGBM との OpenMP 競合を避けるため、このプロセスには LightGBM は存在しない。

使い方:
    python scripts/mpnn_scorer.py ACDEFGHIKL WFYWFYWFYW ...
    -> 1行1スコアで標準出力

sys.path は呼び出し元が正しく設定すること（または PYTHONPATH を使う）。
"""
import sys
import os
import math

# プロジェクトルートを sys.path に追加
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_MPNN_DIR = os.path.join(_ROOT, "models", "proteinmpnn")
if _MPNN_DIR not in sys.path:
    sys.path.insert(0, _MPNN_DIR)

import numpy as np
import torch

_WEIGHTS_PATH = os.path.join(_MPNN_DIR, "v_48_020.pt")
_ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
_AA_TO_IDX = {aa: i for i, aa in enumerate(_ALPHABET)}
_STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def _ideal_helix_coords(length: int) -> np.ndarray:
    RISE = 1.5
    RADIUS = 2.3
    DELTA = math.radians(100.0)
    coords = np.zeros((length, 4, 3), dtype=np.float32)
    for i in range(length):
        angle = i * DELTA
        ca = np.array([RADIUS * math.cos(angle),
                        RADIUS * math.sin(angle),
                        i * RISE], dtype=np.float32)
        tangent = np.array([-math.sin(angle), math.cos(angle),
                             RISE / math.sqrt(RADIUS**2 + RISE**2)], dtype=np.float32)
        tangent /= np.linalg.norm(tangent)
        radial = np.array([math.cos(angle), math.sin(angle), 0.0], dtype=np.float32)
        coords[i, 0] = ca - tangent * 1.46
        coords[i, 1] = ca
        coords[i, 2] = ca + tangent * 1.52
        coords[i, 3] = (ca + tangent * 1.52) + radial * 1.24
    return coords


def _nll_to_score(nll: float, scale: float = 2.0) -> float:
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
def score_sequence(model, sequence: str, device) -> float:
    clean = "".join(aa for aa in sequence.upper() if aa in _STANDARD_AA)
    if len(clean) < 3:
        return 0.5
    L = len(clean)
    coords = _ideal_helix_coords(L)
    X = torch.tensor(coords, device=device).unsqueeze(0)
    S = torch.tensor([[_AA_TO_IDX.get(aa, _AA_TO_IDX["X"]) for aa in clean]],
                      dtype=torch.long, device=device)
    mask = torch.ones(1, L, device=device)
    chain_M = torch.ones(1, L, device=device)
    residue_idx = torch.arange(L, device=device).unsqueeze(0)
    chain_encoding = torch.ones(1, L, dtype=torch.long, device=device)
    randn = torch.zeros(1, L, device=device)
    log_probs = model(X, S, mask, chain_M, residue_idx, chain_encoding, randn)
    nll = float(-log_probs[0, torch.arange(L), S[0]].mean().item())
    return round(_nll_to_score(nll), 4)


if __name__ == "__main__":
    sequences = sys.argv[1:]
    if not sequences:
        sys.exit(0)

    device = torch.device("mps") if torch.backends.mps.is_available() else torch.device("cpu")
    try:
        model = load_model(device)
    except Exception:
        # モデルロード失敗 → 全て 0.5 を返す
        for _ in sequences:
            print(0.5)
        sys.exit(0)

    for seq in sequences:
        try:
            score = score_sequence(model, seq, device)
        except Exception:
            score = 0.5
        print(score)

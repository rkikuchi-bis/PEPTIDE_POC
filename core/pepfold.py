"""
Phase B-2++: ESMFold ローカル推論によるペプチド骨格構造予測

fair-esm パッケージの ESMFold（esmfold_v1）を使ってペプチド配列から
3D 骨格座標（N, CA, C, O）を取得する。

初回実行時にモデル重み（約 700MB）が ~/.cache/torch/hub/ に自動ダウンロードされる。

LightGBM との OpenMP 競合を避けるため、推論は
scripts/esmfold_scorer.py をサブプロセスとして実行する。

上限:
    環境変数 PEPFOLD_MAX_SEQS で一度に予測する最大配列数を設定できる（デフォルト 30）。
    超過分は None を返し、呼び出し元が理想ヘリックスで代替する。
"""
from __future__ import annotations

import os
import pickle
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Sequence

import numpy as np

_ROOT = Path(__file__).resolve().parent.parent
_SCORER_SCRIPT = _ROOT / "scripts" / "esmfold_scorer.py"
_DEFAULT_MAX_SEQS = int(os.environ.get("PEPFOLD_MAX_SEQS", "30"))


def is_esmfold_available() -> bool:
    """fair-esm がインストールされているか確認する。"""
    try:
        import esm  # noqa: F401
        return _SCORER_SCRIPT.exists()
    except ImportError:
        return False


def predict_backbones_batch(
    sequences: Sequence[str],
    max_seqs: int | None = None,
) -> list[np.ndarray | None]:
    """
    複数ペプチド配列の骨格座標を一括予測する（サブプロセス経由）。

    max_seqs を超えた配列は None（理想ヘリックスで代替）を返す。

    Args:
        sequences: ペプチド配列のリスト
        max_seqs:  最大予測数（None の場合は環境変数 PEPFOLD_MAX_SEQS のデフォルト値）

    Returns:
        各配列に対応する (L, 4, 3) numpy 配列または None のリスト
    """
    if not sequences:
        return []

    if max_seqs is None:
        # import 時ではなく呼び出し時に環境変数を読む（実行時の PEPFOLD_MAX_SEQS=0 に対応）
        max_seqs = int(os.environ.get("PEPFOLD_MAX_SEQS", "30"))

    if not is_esmfold_available():
        return [None] * len(sequences)

    # 最大数以内の配列だけサブプロセスに渡す
    seqs_to_predict = list(sequences[:max_seqs])
    n_total = len(sequences)

    tmp_output = None
    try:
        with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as f:
            tmp_output = f.name

        subprocess.run(
            [sys.executable, str(_SCORER_SCRIPT), "--output", tmp_output]
            + seqs_to_predict,
            capture_output=True,
            timeout=600,  # 最大10分（初回モデルDL込み）
        )

        with open(tmp_output, "rb") as f:
            predicted: list[np.ndarray | None] = pickle.load(f)


        # max_seqs を超えた分は None で埋める
        while len(predicted) < len(seqs_to_predict):
            predicted.append(None)
        predicted = predicted[:len(seqs_to_predict)]
        predicted += [None] * (n_total - len(seqs_to_predict))
        return predicted

    except Exception:
        return [None] * n_total
    finally:
        if tmp_output and os.path.exists(tmp_output):
            os.unlink(tmp_output)

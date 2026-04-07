"""
Phase A-2: データセット準備スクリプト

RCSB PDB REST API からペプチド複合体エントリを取得し、
Phase A-1 の特徴量で特徴ベクトル化したデータセットを作成する。

使用方法:
    uv run python scripts/prepare_dataset.py

出力:
    data/peptide_dataset.csv  (陽性例 + 陰性例の特徴ベクトル)
"""
from __future__ import annotations

import random
import sys
import time
from pathlib import Path

import pandas as pd
import requests

# プロジェクトルートを sys.path に追加（core/ をインポートするため）
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from core.rescorer import _protparam_features  # noqa: E402

# ─────────────────────────────────────────────
# 定数
# ─────────────────────────────────────────────
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity"

PEPTIDE_MIN_LEN = 3
PEPTIDE_MAX_LEN = 30
N_POSITIVE = 500   # 陽性例の上限（API取得数）
N_NEGATIVE = 500   # 陰性例の数（ランダム生成）

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
OUTPUT_PATH = ROOT / "data" / "peptide_dataset.csv"


# ─────────────────────────────────────────────
# RCSB Search API: ペプチドエンティティを検索
# ─────────────────────────────────────────────

def search_peptide_entities(n: int = N_POSITIVE) -> list[str]:
    """
    RCSB Search API でペプチド（3〜30残基）を含む複合体エントリを検索する。

    Returns:
        ["4HHB_1", "1ABC_2", ...] 形式のエンティティIDリスト
    """
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_entity_polymer_type",
                        "operator": "exact_match",
                        "value": "Protein"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_sample_sequence_length",
                        "operator": "range",
                        "value": {
                            "from": PEPTIDE_MIN_LEN,
                            "to": PEPTIDE_MAX_LEN,
                            "include_lower": True,
                            "include_upper": True
                        }
                    }
                }
            ]
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": n
            }
        }
    }

    resp = requests.post(RCSB_SEARCH_URL, json=query, timeout=30)
    resp.raise_for_status()

    data = resp.json()
    results = data.get("result_set", [])
    return [r["identifier"] for r in results]


# ─────────────────────────────────────────────
# RCSB Data API: 配列を取得
# ─────────────────────────────────────────────

def fetch_sequence(entity_id: str) -> str | None:
    """
    "4HHB_1" 形式のエンティティIDから配列を取得する。

    Returns:
        ワンレター配列文字列 or None
    """
    parts = entity_id.split("_")
    if len(parts) != 2:
        return None
    pdb_id, ent_id = parts

    url = f"{RCSB_DATA_URL}/{pdb_id}/{ent_id}"
    try:
        resp = requests.get(url, timeout=10)
        if resp.status_code != 200:
            return None
        data = resp.json()
        seq = data.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", "")
        return seq.replace("\n", "").replace(" ", "") if seq else None
    except Exception:
        return None


# ─────────────────────────────────────────────
# 特徴量計算
# ─────────────────────────────────────────────

def featurize(sequence: str, label: int) -> dict | None:
    """
    配列を特徴ベクトルに変換する。
    非標準アミノ酸が多く有効長が不足する場合は None を返す。
    """
    clean = "".join(aa for aa in sequence.upper() if aa in STANDARD_AA)
    if len(clean) < PEPTIDE_MIN_LEN:
        return None

    feats = _protparam_features(clean)
    length = len(clean)

    net_charge = sum(1 for aa in clean if aa in "KRH") - sum(1 for aa in clean if aa in "DE")
    avg_hydrophobicity = sum(1 for aa in clean if aa in "AILMFWVY") / length

    return {
        "sequence": clean,
        "length": length,
        "label": label,
        "net_charge": net_charge,
        "avg_hydrophobicity": round(avg_hydrophobicity, 4),
        **{k: round(v, 4) for k, v in feats.items()},
    }


# ─────────────────────────────────────────────
# 陰性例生成（ランダム）
# ─────────────────────────────────────────────

def generate_negative_examples(n: int) -> list[dict]:
    """
    完全ランダム配列を陰性例（非結合ペプチド）として生成する。
    実際の結合ペプチドとは統計的特性が異なることを前提とした簡易陰性例。
    """
    aa_list = list(STANDARD_AA)
    negatives = []
    while len(negatives) < n:
        length = random.randint(PEPTIDE_MIN_LEN, PEPTIDE_MAX_LEN)
        seq = "".join(random.choices(aa_list, k=length))
        feat = featurize(seq, label=0)
        if feat:
            negatives.append(feat)
    return negatives


# ─────────────────────────────────────────────
# メイン
# ─────────────────────────────────────────────

def main():
    print("=== Phase A-2: データセット準備 ===")

    # 陽性例の検索
    print(f"\n[1/4] RCSB Search API でペプチドエンティティを検索中 (上限 {N_POSITIVE} 件)...")
    entity_ids = search_peptide_entities(N_POSITIVE)
    print(f"  → {len(entity_ids)} 件のエンティティを取得")

    # 配列の取得
    print(f"\n[2/4] 配列を取得中 ({len(entity_ids)} 件)...")
    positives = []
    for i, eid in enumerate(entity_ids):
        if i % 50 == 0:
            print(f"  {i}/{len(entity_ids)}...")
        seq = fetch_sequence(eid)
        if seq:
            feat = featurize(seq, label=1)
            if feat:
                positives.append(feat)
        time.sleep(0.05)  # API レートリミット対策
    print(f"  → 有効な陽性例: {len(positives)} 件")

    # 陰性例の生成
    print(f"\n[3/4] ランダム陰性例を生成中 ({N_NEGATIVE} 件)...")
    negatives = generate_negative_examples(N_NEGATIVE)
    print(f"  → 陰性例: {len(negatives)} 件")

    # データセット保存
    print("\n[4/4] データセットを保存中...")
    df = pd.DataFrame(positives + negatives)
    df = df.sample(frac=1, random_state=42).reset_index(drop=True)  # シャッフル

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_PATH, index=False)

    print(f"\n完了! 保存先: {OUTPUT_PATH}")
    print(f"  陽性例 (label=1): {(df['label'] == 1).sum()} 件")
    print(f"  陰性例 (label=0): {(df['label'] == 0).sum()} 件")
    print(f"  カラム: {list(df.columns)}")
    print(f"\n先頭5行:")
    print(df.head())


if __name__ == "__main__":
    main()

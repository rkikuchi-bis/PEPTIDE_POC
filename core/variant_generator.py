"""
variant_generator.py
既知結合配列（シード）からバリアントを網羅生成するモジュール。
ランダム生成の代わりにリード最適化ワークフローで使用する。
"""
import random
from typing import Dict, List, Optional, Tuple

from core.generator import AMINO_ACIDS, _simple_gen_score

VALID_AAS = set("ACDEFGHIKLMNPQRSTVWY")


def validate_seed_sequence(sequence: str) -> Tuple[bool, str]:
    """
    シード配列を検証する。
    戻り値: (is_valid, error_message)
    """
    seq = sequence.strip().upper()
    if not seq:
        return False, "配列が空です。"
    invalid = set(seq) - VALID_AAS
    if invalid:
        return False, f"無効なアミノ酸文字: {', '.join(sorted(invalid))}"
    if len(seq) < 4:
        return False, "最低4残基必要です。"
    if len(seq) > 30:
        return False, "最大30残基まで対応しています。"
    return True, ""


def estimate_variant_count(
    seed_sequence: str,
    strategies: List[str],
    avoid_residues: Optional[List[str]] = None,
    n_scramble: int = 50,
) -> int:
    """バリアント生成数の概算を返す（UI表示用）。"""
    avoid = set(avoid_residues or [])
    n = 1  # seed 自体
    mutable_aas = [aa for aa in AMINO_ACIDS if aa not in avoid]

    if "single_mutant" in strategies:
        for aa in seed_sequence:
            n += sum(1 for m in mutable_aas if m != aa)

    if "alanine_scan" in strategies:
        n += sum(1 for aa in seed_sequence if aa != "A")

    if "truncation" in strategies:
        min_keep = max(4, len(seed_sequence) - 3)
        for trunc_len in range(len(seed_sequence) - 1, min_keep - 1, -1):
            n_trim = len(seed_sequence) - trunc_len
            if n_trim > 0:
                n += 1  # N末端トランケーション
            n += 1      # C末端トランケーション

    if "scramble" in strategies:
        n += n_scramble

    return n


def generate_variants(
    seed_sequence: str,
    strategies: List[str],
    pocket_charge: str = "neutral",
    pocket_hydrophobicity: str = "medium",
    avoid_residues: Optional[List[str]] = None,
    n_scramble: int = 50,
) -> List[Dict]:
    """
    シード配列からバリアントを生成する。

    strategies に含められる値:
        'single_mutant'  — 各位置を他の19種類のアミノ酸に1残基ずつ変換
        'alanine_scan'   — 各位置をAlaに置換（重要残基の同定用）
        'truncation'     — N/C末端から最大3残基ずつ削る（最低4残基を保持）
        'scramble'       — ランダムシャッフルをn_scramble個生成（対照用）

    戻り値:
        [{"sequence": str, "gen_score": float, "source": str}, ...]
    """
    avoid = set(avoid_residues or [])
    seed = seed_sequence.strip().upper()
    candidates: List[Dict] = []
    seen: set = set()

    def add(seq: str, source: str) -> None:
        if seq and seq not in seen:
            seen.add(seq)
            gen_score = _simple_gen_score(seq, pocket_charge, pocket_hydrophobicity)
            candidates.append({
                "sequence": seq,
                "gen_score": round(gen_score, 4),
                "source": source,
            })

    # シード配列自体を追加
    add(seed, "seed")

    if "single_mutant" in strategies:
        for i, original_aa in enumerate(seed):
            for aa in AMINO_ACIDS:
                if aa == original_aa or aa in avoid:
                    continue
                variant = seed[:i] + aa + seed[i + 1:]
                add(variant, "single_mutant")

    if "alanine_scan" in strategies:
        for i, original_aa in enumerate(seed):
            if original_aa == "A":
                continue
            variant = seed[:i] + "A" + seed[i + 1:]
            add(variant, "alanine_scan")

    if "truncation" in strategies:
        min_keep = max(4, len(seed) - 3)
        for trunc_len in range(len(seed) - 1, min_keep - 1, -1):
            n_trim = len(seed) - trunc_len
            if n_trim > 0:
                # N末端を削る（C末端側を保持）
                add(seed[n_trim:], f"n_truncation_{n_trim}")
            # C末端を削る（N末端側を保持）
            add(seed[:trunc_len], f"c_truncation_{len(seed) - trunc_len}")

    if "scramble" in strategies:
        seq_list = list(seed)
        for _ in range(n_scramble):
            random.shuffle(seq_list)
            variant = "".join(seq_list)
            add(variant, "scramble")

    return candidates

import random
from typing import Dict, List

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

POSITIVE = set("KRH")
NEGATIVE = set("DE")
HYDROPHOBIC = set("AILMFWVY")
POLAR = set("STNQ")

def _build_weighted_pool(
    pocket_charge: str,
    pocket_hydrophobicity: str,
    avoid_residues: List[str] | None = None,
) -> List[str]:
    avoid_residues = set(avoid_residues or [])
    pool: List[str] = []

    for aa in AMINO_ACIDS:
        if aa in avoid_residues:
            continue

        weight = 10

        if pocket_charge == "negative" and aa in POSITIVE:
            weight += 8
        elif pocket_charge == "positive" and aa in NEGATIVE:
            weight += 8

        if pocket_hydrophobicity == "high" and aa in HYDROPHOBIC:
            weight += 8
        elif pocket_hydrophobicity == "low" and aa in POLAR:
            weight += 5

        pool.extend([aa] * max(weight, 1))

    return pool

def _simple_gen_score(sequence: str, pocket_charge: str, pocket_hydrophobicity: str) -> float:
    score = 0.5

    pos_count = sum(1 for x in sequence if x in POSITIVE)
    neg_count = sum(1 for x in sequence if x in NEGATIVE)
    hydrophobic_count = sum(1 for x in sequence if x in HYDROPHOBIC)

    if pocket_charge == "negative":
        score += min(pos_count * 0.03, 0.2)
    elif pocket_charge == "positive":
        score += min(neg_count * 0.03, 0.2)

    if pocket_hydrophobicity == "high":
        score += min(hydrophobic_count * 0.02, 0.2)

    if len(sequence) >= 7:
        score += 0.05

    return min(score, 0.99)

def generate_candidates(
    n: int,
    min_len: int,
    max_len: int,
    pocket_charge: str = "neutral",
    pocket_hydrophobicity: str = "medium",
    avoid_residues: List[str] | None = None,
) -> List[Dict]:
    if min_len > max_len:
        raise ValueError("min_len must be <= max_len")

    pool = _build_weighted_pool(
        pocket_charge=pocket_charge,
        pocket_hydrophobicity=pocket_hydrophobicity,
        avoid_residues=avoid_residues,
    )

    candidates: List[Dict] = []

    for _ in range(n):
        length = random.randint(min_len, max_len)
        seq = "".join(random.choice(pool) for _ in range(length))
        gen_score = _simple_gen_score(seq, pocket_charge, pocket_hydrophobicity)

        candidates.append(
            {
                "sequence": seq,
                "gen_score": round(gen_score, 4),
                "source": "rule_based_v1",
            }
        )

    return candidates
from __future__ import annotations

from typing import List, Sequence

import pandas as pd

AROMATIC_RESIDUES = set("FWYH")
POSITIVE_RESIDUES = set("KRH")
NEGATIVE_RESIDUES = set("DE")


def sanitize_sequence(seq: str) -> str:
    if seq is None:
        return ""
    seq = str(seq).strip().upper()
    return "".join(ch for ch in seq if "A" <= ch <= "Z")


def parse_known_sequences(text: str) -> List[str]:
    if not text:
        return []

    normalized = (
        text.replace(",", "\n")
        .replace(";", "\n")
        .replace("\t", "\n")
        .replace(" ", "\n")
    )
    items = [sanitize_sequence(x) for x in normalized.splitlines()]
    items = [x for x in items if x]

    seen = set()
    out = []
    for x in items:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def get_default_motif_presets() -> dict[str, str]:
    return {
        "Basic amphipathic short helix": "KLAKLAKKLAKLAK",
        "Cationic aromatic-rich": "RWRWRW",
        "Lys-rich binder-like": "KKKLFKKILKYL",
        "Hydrophobic short motif": "FLIVIG",
        "Mixed charged motif": "KRFDEKRF",
    }


def aromatic_ratio(seq: str) -> float:
    seq = sanitize_sequence(seq)
    if not seq:
        return 0.0
    return sum(1 for aa in seq if aa in AROMATIC_RESIDUES) / len(seq)


def charge_pattern(seq: str) -> str:
    seq = sanitize_sequence(seq)
    out = []
    for aa in seq:
        if aa in POSITIVE_RESIDUES:
            out.append("+")
        elif aa in NEGATIVE_RESIDUES:
            out.append("-")
        else:
            out.append("0")
    return "".join(out)


def positional_identity(seq_a: str, seq_b: str) -> float:
    a = sanitize_sequence(seq_a)
    b = sanitize_sequence(seq_b)
    if not a or not b:
        return 0.0
    matches = sum(1 for x, y in zip(a, b) if x == y)
    return matches / max(len(a), len(b))


def charge_pattern_similarity(seq_a: str, seq_b: str) -> float:
    a = charge_pattern(seq_a)
    b = charge_pattern(seq_b)
    if not a or not b:
        return 0.0
    matches = sum(1 for x, y in zip(a, b) if x == y)
    return matches / max(len(a), len(b))


def kmer_set(seq: str, k: int = 2) -> set[str]:
    seq = sanitize_sequence(seq)
    if not seq:
        return set()
    if len(seq) < k:
        return {seq}
    return {seq[i : i + k] for i in range(len(seq) - k + 1)}


def kmer_jaccard(seq_a: str, seq_b: str, k: int = 2) -> float:
    a = kmer_set(seq_a, k=k)
    b = kmer_set(seq_b, k=k)
    if not a and not b:
        return 0.0
    union = a | b
    if not union:
        return 0.0
    return len(a & b) / len(union)


def compare_one_sequence_to_known(seq: str, known_sequences: Sequence[str]) -> dict:
    seq = sanitize_sequence(seq)

    result = {
        "aromatic_ratio": round(aromatic_ratio(seq), 4),
        "charge_pattern": charge_pattern(seq),
        "best_known_motif": "",
        "known_aromatic_gap_min": None,
        "known_charge_pattern_similarity_max": None,
        "known_sequence_identity_max": None,
        "known_kmer_jaccard_max": None,
        "motif_compare_score": None,
    }

    if not seq or not known_sequences:
        return result

    seq_aromatic = aromatic_ratio(seq)

    best_motif = ""
    best_score = -1.0

    aromatic_gaps = []
    charge_sims = []
    identity_sims = []
    kmer_sims = []

    for motif in known_sequences:
        motif = sanitize_sequence(motif)
        if not motif:
            continue

        motif_aromatic = aromatic_ratio(motif)
        aromatic_gap = abs(seq_aromatic - motif_aromatic)
        charge_sim = charge_pattern_similarity(seq, motif)
        identity_sim = positional_identity(seq, motif)
        kmer_sim = kmer_jaccard(seq, motif, k=2)

        compare_score = (
            0.35 * charge_sim
            + 0.35 * identity_sim
            + 0.20 * kmer_sim
            + 0.10 * (1.0 - min(aromatic_gap, 1.0))
        )

        aromatic_gaps.append(aromatic_gap)
        charge_sims.append(charge_sim)
        identity_sims.append(identity_sim)
        kmer_sims.append(kmer_sim)

        if compare_score > best_score:
            best_score = compare_score
            best_motif = motif

    result.update(
        {
            "best_known_motif": best_motif,
            "known_aromatic_gap_min": round(min(aromatic_gaps), 4) if aromatic_gaps else None,
            "known_charge_pattern_similarity_max": round(max(charge_sims), 4) if charge_sims else None,
            "known_sequence_identity_max": round(max(identity_sims), 4) if identity_sims else None,
            "known_kmer_jaccard_max": round(max(kmer_sims), 4) if kmer_sims else None,
            "motif_compare_score": round(best_score, 4) if best_score >= 0 else None,
        }
    )
    return result


def compare_candidates_to_known(df: pd.DataFrame, known_sequences: Sequence[str]) -> pd.DataFrame:
    if df is None:
        return pd.DataFrame()

    out = df.copy()

    if "sequence" not in out.columns:
        return out

    known_sequences = [sanitize_sequence(x) for x in known_sequences if sanitize_sequence(x)]

    if not known_sequences:
        out["aromatic_ratio"] = out["sequence"].fillna("").map(lambda x: round(aromatic_ratio(str(x)), 4))
        out["charge_pattern"] = out["sequence"].fillna("").map(lambda x: charge_pattern(str(x)))
        out["best_known_motif"] = ""
        out["known_aromatic_gap_min"] = None
        out["known_charge_pattern_similarity_max"] = None
        out["known_sequence_identity_max"] = None
        out["known_kmer_jaccard_max"] = None
        out["motif_compare_score"] = None
        return out

    metrics = out["sequence"].fillna("").map(
        lambda x: compare_one_sequence_to_known(str(x), known_sequences)
    )
    metrics_df = pd.DataFrame(list(metrics))

    out = pd.concat([out.reset_index(drop=True), metrics_df.reset_index(drop=True)], axis=1)
    return out
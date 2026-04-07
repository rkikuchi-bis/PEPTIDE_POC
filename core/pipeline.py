from pathlib import Path
from datetime import datetime

import pandas as pd

from core.generator import generate_candidates
from core.filters import apply_filters, add_basic_properties
from core.rescorer import rescore_candidates
from core.utils import save_run_csv
from core.diversity import diversify_candidates
from core.motif_compare import compare_candidates_to_known


def run_pipeline(
    num_candidates: int,
    min_len: int,
    max_len: int,
    pocket_charge: str,
    pocket_hydrophobicity: str,
    avoid_cysteine: bool,
    max_abs_charge: int,
    max_hydrophobicity: float,
    max_repeat_residue: int,
    preferred_len_min: int,
    preferred_len_max: int,
    use_diversity_filter: bool,
    min_diversity_distance: float,
    max_diverse_candidates: int,
    known_sequences: list,
    target_name: str,
) -> tuple[pd.DataFrame, str]:
    candidates = generate_candidates(
        n=num_candidates,
        min_len=min_len,
        max_len=max_len,
        pocket_charge=pocket_charge,
        pocket_hydrophobicity=pocket_hydrophobicity,
        avoid_residues=["C"] if avoid_cysteine else [],
    )

    df = pd.DataFrame(candidates)
    df = add_basic_properties(df)

    filtered_df = apply_filters(
        df,
        min_len=min_len,
        max_len=max_len,
        max_abs_charge=max_abs_charge,
        max_hydrophobicity=max_hydrophobicity,
        max_repeat_residue=max_repeat_residue,
        remove_near_duplicates=True,
    )

    rescored_df = rescore_candidates(
        filtered_df,
        pocket_charge=pocket_charge,
        pocket_hydrophobicity=pocket_hydrophobicity,
        preferred_len_min=preferred_len_min,
        preferred_len_max=preferred_len_max,
    )

    rescored_df = rescored_df.sort_values(
        by=["final_score", "rescoring_score", "property_score", "gen_score"],
        ascending=False,
    ).reset_index(drop=True)

    if use_diversity_filter:
        rescored_df = diversify_candidates(
            rescored_df,
            sequence_col="sequence",
            score_col="final_score",
            min_normalized_distance=min_diversity_distance,
            max_candidates=max_diverse_candidates,
        )
    else:
        rescored_df = rescored_df.copy()
        rescored_df["diversity_kept"] = True
        if "diversity_min_distance" not in rescored_df.columns:
            rescored_df["diversity_min_distance"] = 1.0

    rescored_df = compare_candidates_to_known(
        rescored_df,
        known_sequences=known_sequences,
    )

    rescored_df = rescored_df.reset_index(drop=True)
    rescored_df.insert(0, "rank", range(1, len(rescored_df) + 1))

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path("outputs")
    output_file = output_dir / f"{target_name}_{timestamp}.csv"
    saved_path = save_run_csv(rescored_df, str(output_file))

    return rescored_df, saved_path

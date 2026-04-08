"""ML + ProteinMPNN + rescore_candidates 統合テスト"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd
from core.rescorer import rescore_candidates

df = pd.DataFrame({
    "sequence": ["ACDEFGHIKL", "WFYWFYWFYW", "KRKRKRKRKR"],
    "length": [10, 10, 10],
    "gen_score": [0.7, 0.6, 0.5],
    "property_score": [0.8, 0.5, 0.6],
    "net_charge": [0, 0, 4],
    "avg_hydrophobicity": [0.4, 0.8, 0.1],
})

result = rescore_candidates(df, pocket_charge="negative", pocket_hydrophobicity="high")
cols = ["sequence", "final_score", "ml_score", "proteinmpnn_score", "rescoring_score"]
print(result[cols].to_string(index=False))
print()
print("ml_score_available:", result["ml_score_available"].iloc[0])
print("proteinmpnn_score_available:", result["proteinmpnn_score_available"].iloc[0])
print("ALL OK")

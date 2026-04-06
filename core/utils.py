from pathlib import Path
import pandas as pd

def save_run_csv(df: pd.DataFrame, path: str) -> str:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    return str(output_path.resolve())
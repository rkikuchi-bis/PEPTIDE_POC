# Project: Peptide Discovery PoC

## Goal
タンパク質ポケット構造（PDB/mmCIF）を入力として、結合候補ペプチド配列を生成・フィルタ・スコアリング・可視化するStreamlitアプリのPoC開発。
創薬初期の仮説生成ツールとして、段階的にスコアリングの科学的精度を高めることを目標とする。

---

## Architecture

```
PEPTIDE_POC/
├── app.py                  # 薄いオーケストレーター（約60行）
├── core/
│   ├── generator.py        # バイアス付きランダム配列生成
│   ├── filters.py          # 物性フィルタ（電荷・疎水性・繰り返し）
│   ├── rescorer.py         # BioPython ProtParam スコアリング（Phase A-1）
│   ├── diversity.py        # Hamming距離による多様性制御
│   ├── pdb_utils.py        # PDB/mmCIF パース・ポケット解析
│   ├── motif_compare.py    # 既知モチーフとの比較
│   ├── pipeline.py         # 生成→フィルタ→スコア→多様性→比較のオーケストレーター
│   ├── rcsb_client.py      # RCSB API クライアント
│   ├── structure_scorer.py # 構造優先度スコアリング
│   └── utils.py            # CSV保存
├── ui/
│   ├── sidebar.py          # render_sidebar() → dict を返す
│   ├── results.py          # render_results(result_df, pdb_summary)
│   └── structure_viewer.py # py3Dmol 3D可視化
├── data/
│   └── example_targets/
├── prompts/                # Claude再起動用プロンプト
├── CLAUDE.md               # このファイル（Claude用）
├── CONTEXT.md              # 現在の開発状態
└── TODO.md                 # 短期タスク
```

### データフロー
```
User Input (Sidebar)
  → Structure Upload / RCSB Search
  → PDB Parse & Pocket Analysis (pdb_utils)
  → run_pipeline():
      generate_candidates()
      → add_basic_properties()
      → apply_filters()
      → rescore_candidates()       ← Phase A-1: ProtParam
      → diversify_candidates()
      → compare_candidates_to_known()
  → render_results() + render_viewer_section()
```

---

## Current Features

- **構造入力**: ローカルPDB/mmCIFアップロード または RCSB検索・ダウンロード
- **ポケット解析**: Manual region（残基範囲指定）または Ligand neighborhood（半径指定）
- **候補生成**: ポケット電荷・疎水性バイアス付きランダム生成
- **フィルタ**: 電荷・疎水性・繰り返し残基・重複除去
- **スコアリング（Phase A-1）**: BioPython ProtParam
  - 等電点 pI → 電荷マッチスコア
  - GRAVY → 疎水性マッチスコア
  - 不安定性指数 → 安定性スコア
  - 芳香族性・ヘリックス傾向・配列多様性 → 複雑性スコア
- **多様性制御**: Hamming距離ベースのgreedy選択
- **モチーフ比較**: 既知配列との同一性・k-mer Jaccard・電荷パターン類似度
- **3D可視化**: py3Dmol（チェーン・ポケット・リガンドを色分け表示）
- **CSV出力**: `outputs/{target}_{timestamp}.csv`

---

## Scoring Roadmap

| Phase | 内容 | 状態 |
|---|---|---|
| A-1 | BioPython ProtParam（pI, GRAVY, 不安定性指数） | ✅ 完了 |
| A-2 | PepBDB学習モデル（RandomForest/LightGBM） | 未着手 |
| B-1 | AutoDock Vina ドッキングスコア | 未着手 |
| B-2 | ProteinMPNN 配列逆設計（M4 MPS対応） | 未着手 |

---

## Constraints

- Python 3.11〜3.12（`pyproject.toml` で指定）
- パッケージ管理: `uv`（`pip` は使わない）
- 実行: `uv run streamlit run app.py`
- `core/` 内のモジュールは Streamlit に依存しない（純粋なPythonのみ）
- `ui/` 内のモジュールは Streamlit に依存してよい
- `app.py` は薄いオーケストレーターに保つ（UI・ロジックを混在させない）
- 開発環境: MacBook Pro M4 MAX 128GB（MPS利用可能、NVIDIA不要）

---

## Coding Rules

- 新しいビジネスロジックは必ず `core/` に追加する
- 新しいUI要素は必ず `ui/` に追加する
- `app.py` には imports・session_state初期化・3つの関数呼び出しのみ
- 依存追加は `uv add <package>`（pyproject.tomlが自動更新される）
- 関数は単一責任。副作用（st.write等）をコアロジックに混ぜない
- コメントは日本語でもよい

---

## Output Format

- 結果DataFrameは `outputs/{target_name}_{timestamp}.csv` に保存
- 主要カラム:
  `rank, sequence, final_score, rescoring_score, gen_score, property_score,`
  `length, net_charge, avg_hydrophobicity,`
  `isoelectric_point, gravy, instability_index, aromaticity,`
  `molecular_weight, helix_fraction, turn_fraction, sheet_fraction,`
  `charge_match_score, hydrophobic_match_score, stability_score, complexity_score,`
  `diversity_kept, diversity_min_distance,`
  `best_known_motif, known_sequence_identity_max, motif_compare_score,`
  `rescoring_notes`

# 🧬 Peptide Discovery PoC

**AI-assisted peptide candidate generation and scoring for early drug discovery**  
**創薬初期の仮説生成を支援する、AIペプチド候補設計・スコアリングツール**

---

## Overview / 概要

**[English]**  
Peptide Discovery PoC is a Streamlit-based web application that takes a protein pocket structure (PDB/mmCIF) as input and generates, filters, scores, and visualizes binding peptide candidates. Designed as a hypothesis-generation tool for early-stage drug discovery, it integrates multiple scoring layers — from physicochemical properties to machine learning and structural design — and provides Japanese-language explanations of each candidate's rationale.

**[日本語]**  
Peptide Discovery PoC は、タンパク質ポケット構造（PDB/mmCIF）を入力として、結合候補ペプチド配列を生成・フィルタ・スコアリング・可視化する Streamlit 製 Web アプリです。創薬初期の仮説生成ツールとして、物性スコアから機械学習・構造解析まで複数のスコアリング層を統合し、各候補の推薦理由を日本語で出力します。

---

## Key Features / 主な機能

| Feature | Description |
|---|---|
| **Simple / Expert mode** | Simple: 構造をアップロードして ▶ Run だけ。Expert: 全パラメータ手動制御 |
| **Structure Input** | ローカル PDB/mmCIF アップロード または RCSB PDB から直接検索・ダウンロード |
| **Pocket Analysis** | Manual residue range または Ligand neighborhood（リガンド自動検出） |
| **Candidate Generation** | ポケット電荷・疎水性バイアス付きランダム配列生成 |
| **Multi-layer Scoring** | ProtParam (A-1) → LightGBM ML (A-2) → ProteinMPNN (B-2/B-2+/B-2++) → Vina Docking (B-1) |
| **Selectivity Score** | ターゲット vs. オフターゲットの結合選択性スコア（Phase C-1）|
| **Off-target DB** | 14ターゲット対応の既知 Off-target データベース内蔵（EGFR, hERG, CYP3A4 等）|
| **Explainability** | 各候補に対する日本語推薦理由（pI / GRAVY / 安定性 / ML / MPNN / ドッキング統合）|
| **3D Visualization** | py3Dmol インタラクティブビューア（ポケット・リガンド・ペプチド候補を色分け表示）|
| **CSV Export** | 全スコアを含む結果テーブルを `outputs/` に保存 |

---

## Scoring Pipeline / スコアリングパイプライン

```
Candidate Generation  (pocket charge/hydrophobicity biased random)
  ↓
Filter  (charge / hydrophobicity / repeat residue / deduplication)
  ↓
Rescoring
  ├── Phase A-1: BioPython ProtParam  (pI, GRAVY, instability index, aromaticity)
  ├── Phase A-2: LightGBM ML model   (RCSB peptide complex data, CV AUC = 0.891)
  └── Phase B-2/B-2+/B-2++: ProteinMPNN
        B-2:   structure-free (ideal α-helix surrogate)
        B-2+:  receptor-conditioned (pocket centroid + ideal helix)
        B-2++: receptor-conditioned + ESMFold predicted backbone
  ↓
Diversity Filter  (Hamming distance-based greedy selection)
  ↓
ESMFold Re-scoring  (B-2++: top diverse candidates re-scored with predicted backbone)
  ↓
Motif Comparison  (identity / k-mer Jaccard / charge pattern similarity)
  ↓
Selectivity Score (Phase C-1)  [optional, Expert mode]
  selective_final_score = final_score + λ × selectivity_score  (default λ = 0.3)
  ↓
Explainability (Direction A)  — 日本語推薦理由
  ↓
Optional Docking (Phase B-1): AutoDock Vina  [local only, Expert mode]
```

---

## Quick Start / デモ

> Try with the classic HIV protease structure **1HSG** from RCSB PDB.  
> RCSB から **1HSG**（HIVプロテアーゼ）を検索してそのまま試せます。

### Simple mode（推奨）

1. サイドバーで **Simple** モードを選択
2. **"Search RCSB by target"** → `HIV protease` と入力 → **Search RCSB**
3. `1HSG` を選択 → **Use selected structure**
4. Chain と Ligand が自動推奨される（Chain A / Ligand MK1）
5. **▶ Run** を押す
6. ランク付き候補・スコア・3D ビューア・推薦理由を確認

### Expert mode（全機能）

- Pocket bias / Generation 数 / Diversity filter / Docking / Selectivity を手動制御
- Known off-target DB から選択し、選択性スコアを計算

---

## Architecture / アーキテクチャ

```
peptide_poc/
├── app.py                    # 薄いオーケストレーター（~60行）
├── core/
│   ├── generator.py          # バイアス付きランダム配列生成
│   ├── filters.py            # 物性フィルタ
│   ├── rescorer.py           # マルチレイヤースコアリング統合
│   ├── ml_scorer.py          # LightGBM 推論（Phase A-2）
│   ├── proteinmpnn.py        # ProteinMPNN スコアリング（Phase B-2/B-2+/B-2++）
│   ├── pepfold.py            # ESMFold ローカル推論（Phase B-2++）
│   ├── docking.py            # AutoDock Vina ドッキング（Phase B-1、ローカルのみ）
│   ├── selectivity.py        # 選択性スコア（Phase C-1）
│   ├── explainer.py          # 日本語推薦理由生成（Direction A）
│   ├── helix_utils.py        # 理想αヘリックス座標生成（3D Viewer 用）
│   ├── offtarget_db.py       # 既知 Off-target DB（14ターゲット）
│   ├── diversity.py          # Hamming 距離多様性制御
│   ├── pdb_utils.py          # PDB/mmCIF パース・ポケット解析
│   ├── motif_compare.py      # 既知モチーフ比較
│   ├── pipeline.py           # パイプラインオーケストレーター
│   ├── structure_scorer.py   # 構造優先度スコアリング
│   ├── rcsb_client.py        # RCSB API クライアント
│   └── utils.py              # CSV 保存ユーティリティ
├── ui/
│   ├── sidebar.py            # render_sidebar() → dict
│   ├── results.py            # render_results()
│   └── structure_viewer.py   # py3Dmol 3D ビューア + ペプチド重畳
├── models/
│   ├── peptide_classifier.joblib    # 学習済み LightGBM モデル
│   └── proteinmpnn/v_48_020.pt      # ProteinMPNN 重み (6.4 MB)
├── scripts/
│   ├── prepare_dataset.py    # 学習データ生成（RCSB API）
│   ├── train_classifier.py   # LightGBM 学習スクリプト
│   ├── mpnn_scorer.py        # ProteinMPNN スコアラー（サブプロセス用）
│   ├── mpnn_scorer_receptor.py  # 受容体条件付き ProteinMPNN
│   └── esmfold_scorer.py     # ESMFold 推論スクリプト
└── data/
    └── peptide_dataset.csv   # 学習用データセット（陽性 485 件 + 陰性 500 件）
```

---

## Installation / セットアップ

### Requirements / 必要環境

- Python 3.11 または 3.12
- [uv](https://github.com/astral-sh/uv)（パッケージマネージャー）

### Local setup / ローカル実行

```bash
git clone https://github.com/rkikuchi-bis/PEPTIDE_POC.git
cd PEPTIDE_POC

uv python pin 3.11
uv sync
uv run streamlit run app.py
```

### Optional: AutoDock Vina（Phase B-1）

```bash
# macOS (arm64)
brew install openbabel

# Vina バイナリを bin/vina に配置
# https://github.com/ccsb-scripps/AutoDock-Vina/releases
```

### Optional: ESMFold（Phase B-2++）

```bash
# fair-esm と依存パッケージのインストール（uv 経由）
# モデルは初回実行時に自動ダウンロード（~2GB）
# ~/.cache/torch/hub/checkpoints/esmfold_3B_v1.pt
```

---

## Scoring Details / スコアリング詳細

### Phase A-1: BioPython ProtParam

| Feature | Usage |
|---|---|
| Isoelectric point (pI) | ポケット電荷との補完性スコア |
| GRAVY | ポケット疎水性とのマッチスコア |
| Instability index | in vivo 安定性推定スコア |
| Aromaticity, helix fraction | 配列複雑性スコア |

### Phase A-2: ML Scoring (LightGBM)

- 学習データ: RCSB PDB ペプチド複合体 陽性 485 件 + ランダム陰性 500 件
- 5-fold CV AUC: **0.891** / Holdout AUC: **0.881**
- 特徴量: Phase A-1 全特徴量 + 長さ・電荷・疎水性

### Phase B-2 / B-2+ / B-2++: ProteinMPNN

| Mode | Description |
|---|---|
| B-2 | 構造フリー（理想αヘリックス骨格でスコア） |
| B-2+ | 受容体条件付き（ポケット centroid に配置した理想ヘリックス） |
| B-2++ | 受容体条件付き + ESMFold 予測骨格（diversity 後の少数候補に適用） |

### Phase B-1: AutoDock Vina（ローカル限定）

- ≤5 残基: 柔軟ドッキング / ≥6 残基: 剛体ドッキング（自動切り替え）
- obabel 経由で PDB → PDBQT 変換
- スコア単位: kcal/mol（負値が強い結合を示す）

### Phase C-1: Selectivity Score

```
selectivity_score        = rescoring_score(target) − rescoring_score(off-target)
selective_final_score    = final_score + λ × selectivity_score   (default λ = 0.3)
```

- 🟢 selectivity ≥ 0.10: ターゲット選択的
- 🟡 −0.05 < selectivity < 0.10: 中程度
- 🔴 selectivity ≤ −0.05: オフターゲット親和性リスク

---

## Cloud Deployment / クラウドデプロイ

Streamlit Cloud に `requirements.txt` でデプロイ可能です。  
ESMFold と AutoDock Vina はリソース制約により自動無効化されます。

| Feature | Streamlit Cloud | Local |
|---|---|---|
| Structure input & 3D viewer | ✅ | ✅ |
| ProtParam + ML scoring | ✅ | ✅ |
| ProteinMPNN B-2 / B-2+ | ✅ (CPU) | ✅ (MPS) |
| Selectivity + Explainability | ✅ | ✅ |
| AutoDock Vina docking (B-1) | ❌ | ✅ |
| ESMFold backbone (B-2++) | ❌ | ✅ |

---

## Limitations / 制約事項

**[English]**  
This is a proof-of-concept tool for hypothesis generation. Candidate sequences are generated by biased random sampling, not structural design. Scores are for relative ranking purposes and do not guarantee actual binding activity. Selectivity scores based on physicochemical differences have limited resolution for closely related protein families (Phase C-2 docking-based selectivity is planned for higher resolution). All results should be validated by experimental assay.

**[日本語]**  
本ツールは仮説生成のための PoC（概念実証）です。候補配列はバイアス付きランダム生成であり、構造ベース設計ではありません。スコアは相対的なランキング指標であり、実際の結合活性を保証するものではありません。物性差分ベースの選択性スコアは、同一ファミリー内タンパク質の識別精度に限界があります（ドッキングスコア差分ベースの Phase C-2 を将来計画）。実験による検証が必要です。

---

## Tech Stack / 使用技術

`Python 3.11` · `Streamlit` · `BioPython` · `LightGBM` · `PyTorch (MPS)` · `ProteinMPNN` · `ESMFold (fair-esm)` · `AutoDock Vina` · `RDKit` · `py3Dmol` · `RCSB PDB API`

---

## Author / 著者

Developed as a portfolio project demonstrating multi-layer AI scoring in early-stage drug discovery.  
創薬初期における多段階 AI スコアリングのポートフォリオとして開発。

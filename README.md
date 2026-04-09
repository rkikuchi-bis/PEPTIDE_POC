# Peptide Discovery PoC

**AI-assisted peptide candidate generation and scoring for early drug discovery**  
**創薬初期の仮説生成を支援する、AIペプチド候補設計・スコアリングツール**

---

## Overview / 概要

**[English]**  
Peptide Discovery PoC is a Streamlit-based web application that takes a protein pocket structure (PDB/mmCIF) as input and generates, filters, scores, and visualizes binding peptide candidates. Designed as a hypothesis-generation tool for early-stage drug discovery, it integrates multiple scoring layers — from physicochemical properties to machine learning and structural analysis — and provides Japanese-language explanations of each candidate's rationale.

**[日本語]**  
Peptide Discovery PoC は、タンパク質ポケット構造（PDB/mmCIF）を入力として、結合候補ペプチド配列を生成・フィルタ・スコアリング・可視化する Streamlit 製 Web アプリです。創薬初期の仮説生成ツールとして、物性スコアから機械学習・構造解析まで複数のスコアリング層を統合し、各候補の推薦理由を日本語で出力します。

---

## Key Features / 主な機能

| Feature | Description |
|---|---|
| **Structure Input** | Upload local PDB/mmCIF or search & download directly from RCSB PDB |
| **Pocket Analysis** | Manual residue range or ligand neighborhood (auto-detected) |
| **Candidate Generation** | Charge- and hydrophobicity-biased random sequence generation |
| **Multi-layer Scoring** | ProtParam (A-1) → ML/LightGBM (A-2) → ProteinMPNN (B-2) → Docking/Vina (B-1) |
| **Selectivity Score** | Quantifies on-target vs. off-target preference (Phase C-1) |
| **Explainability** | Japanese-language rationale for each top candidate (Direction A) |
| **3D Visualization** | Interactive structure viewer with pocket and ligand highlighting |
| **CSV Export** | Full results table saved to `outputs/` |

---

## Scoring Pipeline / スコアリングパイプライン

```
Candidate Generation (biased random)
  ↓
Filter  (charge / hydrophobicity / repeat residue / deduplication)
  ↓
Rescoring
  ├── Phase A-1: BioPython ProtParam  (pI, GRAVY, instability index, aromaticity)
  ├── Phase A-2: LightGBM ML model   (trained on RCSB peptide complex data, AUC = 0.891)
  └── Phase B-2: ProteinMPNN         (sequence designability score, structure-free)
  ↓
Diversity Filter  (Hamming distance-based greedy selection)
  ↓
Motif Comparison  (identity / k-mer Jaccard / charge pattern similarity)
  ↓
Selectivity Score (Phase C-1)  [optional]
  ↓
Explainability    (Direction A) — Japanese rationale per candidate
```

---

## Demo / デモ

> Try with the classic HIV protease structure **1HSG** from RCSB PDB.  
> RCSB から **1HSG**（HIVプロテアーゼ）を検索してそのまま試せます。

1. Select **"Search RCSB by target"** → type `HIV protease` → click **Search RCSB**
2. Select `1HSG` → click **Use selected structure**
3. Choose pocket mode: **Ligand neighborhood** (ligand: MK1)
4. Click **Generate and Filter**
5. Review ranked candidates, scores, 3D viewer, and Japanese explanations

---

## Architecture / アーキテクチャ

```
peptide_poc/
├── app.py                  # Thin orchestrator (~60 lines)
├── core/
│   ├── generator.py        # Biased random sequence generation
│   ├── filters.py          # Physicochemical filters
│   ├── rescorer.py         # Multi-layer rescoring integration
│   ├── ml_scorer.py        # LightGBM inference (Phase A-2)
│   ├── proteinmpnn.py      # ProteinMPNN scoring via subprocess (Phase B-2)
│   ├── docking.py          # AutoDock Vina docking (Phase B-1, local only)
│   ├── selectivity.py      # Selectivity score (Phase C-1)
│   ├── explainer.py        # Japanese rationale generation (Direction A)
│   ├── diversity.py        # Hamming-distance diversity control
│   ├── pdb_utils.py        # PDB/mmCIF parsing & pocket analysis
│   ├── motif_compare.py    # Known motif comparison
│   ├── pipeline.py         # Pipeline orchestration
│   └── rcsb_client.py      # RCSB API client
├── ui/
│   ├── sidebar.py          # render_sidebar() → dict
│   ├── results.py          # render_results()
│   └── structure_viewer.py # py3Dmol 3D viewer
├── models/
│   ├── peptide_classifier.joblib   # Trained LightGBM model
│   └── proteinmpnn/v_48_020.pt     # ProteinMPNN weights (6.4 MB)
└── requirements.txt
```

---

## Installation / セットアップ

### Requirements / 必要環境
- Python 3.11 or 3.12
- [uv](https://github.com/astral-sh/uv) (package manager)

### Local setup / ローカル実行

```bash
git clone https://github.com/rkikuchi-bis/PEPTIDE_POC.git
cd PEPTIDE_POC

uv python pin 3.11
uv sync
uv run streamlit run app.py
```

### Optional local features / ローカル限定機能

AutoDock Vina docking (Phase B-1) requires:

```bash
# macOS (arm64)
brew install openbabel
# Download vina binary to bin/vina (mac_aarch64)
# https://github.com/ccsb-scripps/AutoDock-Vina/releases
```

---

## Scoring Details / スコアリング詳細

### Phase A-1: BioPython ProtParam

| Feature | Use |
|---|---|
| Isoelectric point (pI) | Charge complementarity with pocket |
| GRAVY | Hydrophobic match with pocket |
| Instability index | In vivo stability estimate |
| Aromaticity, helix fraction | Structural complexity score |

### Phase A-2: ML Scoring (LightGBM)

- Training data: 485 peptide-positive complexes from RCSB PDB + 500 random negatives  
- 5-fold CV AUC: **0.891** / Holdout AUC: **0.881**  
- Features: all Phase A-1 physicochemical descriptors

### Phase B-2: ProteinMPNN

- Sequence designability score using vanilla ProteinMPNN weights (v_48_020)  
- Structure-free mode: ideal α-helix backbone as surrogate

### Phase C-1: Selectivity Score

```
selectivity_score = rescoring_score(target) − rescoring_score(off-target)
selective_final_score = final_score + λ × selectivity_score  (default λ = 0.3)
```

Positive selectivity score indicates target-selective candidates.  
選択性スコアが正値の候補はターゲット選択的と評価されます。

---

## Cloud Deployment / クラウドデプロイ

The app is deployable to [Streamlit Cloud](https://share.streamlit.io) using `requirements.txt`.  
Streamlit Cloud では ESMFold と AutoDock Vina は自動的に無効化されます。

| Feature | Streamlit Cloud | Local |
|---|---|---|
| Structure input & 3D viewer | ✅ | ✅ |
| ProtParam + ML scoring | ✅ | ✅ |
| ProteinMPNN scoring | ✅ (CPU) | ✅ |
| Selectivity + Explainability | ✅ | ✅ |
| AutoDock Vina docking | ❌ | ✅ |
| ESMFold backbone (B-2++) | ❌ | ✅ |

---

## Limitations / 制約事項

**[English]**  
This is a proof-of-concept tool for hypothesis generation. Candidate sequences are generated by biased random sampling, not structural design. Scores are for relative ranking and do not guarantee actual binding activity. Results should be validated by experimental assay.

**[日本語]**  
本ツールは仮説生成のためのPoC（概念実証）です。候補配列はバイアス付きランダム生成であり、構造ベース設計ではありません。スコアは相対的なランキング指標であり、実際の結合活性を保証するものではありません。実験による検証が必要です。

---

## Tech Stack / 使用技術

`Python 3.11` · `Streamlit` · `BioPython` · `LightGBM` · `PyTorch` · `ProteinMPNN` · `py3Dmol` · `RCSB PDB API`

---

## Author / 著者

Developed as a portfolio project demonstrating AI application in early-stage drug discovery.  
創薬初期AIアプリケーションのポートフォリオとして開発。

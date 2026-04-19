# 🧬 Peptide Discovery PoC

**AI-assisted peptide candidate generation and scoring for early drug discovery**

![Python](https://img.shields.io/badge/python-3.11%20%7C%203.12-blue?logo=python&logoColor=white)
![Streamlit](https://img.shields.io/badge/Streamlit-1.45%2B-FF4B4B?logo=streamlit&logoColor=white)
![License](https://img.shields.io/badge/license-MIT-green)
![Platform](https://img.shields.io/badge/platform-macOS%20%7C%20Linux-lightgrey)

---

## What This Is

A Streamlit web application that takes a protein pocket structure (PDB/mmCIF) as input and generates, filters, scores, and visualizes binding peptide candidates.

**Intended use:** Hypothesis generation at the start of a drug discovery campaign — before running expensive wet-lab assays.

**Key design decisions:**
- **Layered scoring**: physics → ML → generative structure models → docking. Each layer refines the ranking without blocking earlier results.
- **Graceful degradation**: every optional module (ML, ProteinMPNN, ESMFold, Vina) has a fallback so the app works with zero local setup.
- **Selectivity-first**: built-in off-target database lets you penalize promiscuous binders at design time.

---

## Demo

> Try it with the HIV protease structure **1HSG** from RCSB PDB — no local files needed.

**Simple mode (recommended first run):**
1. Select **Simple** mode in the sidebar
2. Search RCSB: type `HIV protease` → click **Search RCSB**
3. Select `1HSG` → **Use selected structure**
4. Chain A and Ligand MK1 are auto-recommended
5. Click **▶ Run** — ranked candidates, scores, 3D viewer, and Japanese-language explanations appear

**Expert mode** unlocks manual control over generation count, diversity filter, docking, and selectivity scoring.

---

## Scoring Pipeline

```
Input: protein pocket structure (PDB / mmCIF)
  │
  ▼
Candidate Generation
  Biased random sampling (pocket charge + hydrophobicity weighted)
  — or — Variant generation (single-mutant / alanine scan / truncation)
  │
  ▼
Property Filters
  Length · Net charge · Avg hydrophobicity · Repeat-residue · Deduplication
  │
  ▼
Multi-layer Rescoring
  Phase A-1  BioPython ProtParam    pI · GRAVY · instability index · aromaticity
  Phase A-2  LightGBM ML model      CV AUC = 0.891 (RCSB peptide complex dataset)
  Phase B-2  ProteinMPNN            structure-free (ideal α-helix backbone)
  Phase B-2+ ProteinMPNN            receptor-conditioned (pocket centroid + ideal helix)
  Phase B-2++ ProteinMPNN           receptor-conditioned + ESMFold predicted backbone
  │
  ▼
Diversity Filter   Hamming distance-based greedy selection
  │
  ▼
ESMFold Re-scoring (B-2++: top diverse candidates only)
  │
  ▼
Motif Comparison   Sequence identity · k-mer Jaccard · charge pattern similarity
  │
  ▼
Optional post-processing
  Phase C-1  Selectivity score      rescoring diff (target − off-target)
  Phase B-1  AutoDock Vina          rigid/flexible docking  [local only]
  Phase C-2  Docking selectivity    Vina score diff (off-target − target)  [local only]
  │
  ▼
Explainability (Direction A)
  Japanese + English natural-language rationale for each top candidate
```

---

## Architecture

```
peptide_poc/
├── app.py                        # Thin orchestrator (~70 lines)
├── core/                         # Pure Python — no Streamlit dependency
│   ├── pipeline.py               # End-to-end orchestration
│   ├── generator.py              # Pocket-biased sequence generation
│   ├── variant_generator.py      # Mutation / truncation variants
│   ├── filters.py                # Property filters
│   ├── rescorer.py               # Multi-layer scoring (A-1 → B-2++)
│   ├── ml_scorer.py              # LightGBM inference (Phase A-2)
│   ├── proteinmpnn.py            # ProteinMPNN subprocess wrapper (B-2/B-2+/B-2++)
│   ├── pepfold.py                # ESMFold backbone prediction (B-2++)
│   ├── docking.py                # AutoDock Vina wrapper (B-1)
│   ├── selectivity.py            # Selectivity scores (C-1 / C-2)
│   ├── explainer.py              # Japanese / English explanations (Direction A)
│   ├── diversity.py              # Hamming distance diversity filter
│   ├── motif_compare.py          # Known-motif comparison
│   ├── pdb_utils.py              # PDB / mmCIF parsing & pocket analysis
│   ├── rcsb_client.py            # RCSB PDB REST API client
│   ├── offtarget_db.py           # Built-in off-target DB (14 targets)
│   ├── admet_scorer.py           # Heuristic bioactivity scoring
│   ├── helix_utils.py            # Ideal α-helix coordinate generation
│   ├── structure_scorer.py       # RCSB result ranking
│   └── utils.py                  # CSV export
├── ui/                           # Streamlit UI modules
│   ├── sidebar.py                # render_sidebar() → params dict
│   ├── actions.py                # Post-pipeline action panels
│   ├── results.py                # Results table & candidate detail panel
│   └── structure_viewer.py       # py3Dmol 3D visualization
├── models/
│   ├── peptide_classifier.joblib # Trained LightGBM weights (Phase A-2)
│   └── proteinmpnn/v_48_020.pt   # ProteinMPNN weights (6.4 MB)
├── scripts/
│   ├── prepare_dataset.py        # RCSB training data pipeline
│   ├── train_classifier.py       # LightGBM 5-fold CV training
│   ├── mpnn_scorer.py            # ProteinMPNN subprocess runner
│   ├── mpnn_scorer_receptor.py   # Receptor-conditioned MPNN runner
│   └── esmfold_scorer.py         # ESMFold inference runner
└── data/
    └── peptide_dataset.csv       # 485 positive + 500 negative examples
```

**Design note:** ProteinMPNN and ESMFold run as subprocesses to avoid an OpenMP thread conflict between PyTorch and LightGBM on Apple Silicon (MPS).

---

## Installation

### Requirements

- Python 3.11 or 3.12
- [uv](https://github.com/astral-sh/uv) package manager

### Quick start

```bash
git clone https://github.com/rkikuchi-bis/PEPTIDE_POC.git
cd PEPTIDE_POC

uv python pin 3.11
uv sync
uv run streamlit run app.py
```

### Optional: AutoDock Vina (Phase B-1)

```bash
# macOS
brew install openbabel

# Download the Vina binary for your platform and place it at bin/vina
# https://github.com/ccsb-scripps/AutoDock-Vina/releases
```

### Optional: ESMFold (Phase B-2++)

```bash
# Dependencies are included in pyproject.toml.
# The model checkpoint (~2 GB) is downloaded automatically on first use to:
# ~/.cache/torch/hub/checkpoints/esmfold_3B_v1.pt
#
# Control the number of sequences sent to ESMFold:
export PEPFOLD_MAX_SEQS=30  # default; set to 0 to disable
```

### Streamlit Cloud deployment

The `requirements.txt` targets Streamlit Cloud (CPU-only, no Vina or ESMFold):

```bash
# Streamlit Cloud reads requirements.txt automatically.
# ESMFold and AutoDock Vina are disabled in cloud environments.
```

| Feature | Streamlit Cloud | Local |
|---|---|---|
| Structure input & 3D viewer | ✅ | ✅ |
| ProtParam + ML scoring | ✅ | ✅ |
| ProteinMPNN B-2 / B-2+ | ✅ (CPU) | ✅ (MPS) |
| Selectivity + Explainability | ✅ | ✅ |
| AutoDock Vina (B-1) | ❌ | ✅ |
| ESMFold backbone (B-2++) | ❌ | ✅ |

---

## Scoring Details

### Phase A-1: BioPython ProtParam

| Feature | Scoring rationale |
|---|---|
| Isoelectric point (pI) | Complementarity with pocket charge at physiological pH 7.4 |
| GRAVY (Kyte-Doolittle) | Match to pocket hydrophobicity; sigmoid-normalized |
| Instability index | In vivo stability estimate; threshold at 40 |
| Aromaticity | π-π stacking potential (F/W/Y fraction) |
| Helix fraction | α-helix propensity; peak reward near 35% |

### Phase A-2: LightGBM ML Model

- **Training data:** RCSB PDB peptide complexes — 485 positive + 500 negative examples
- **Validation:** 5-fold CV AUC **0.891** / holdout AUC **0.881**
- **Features:** all Phase A-1 physicochemical descriptors + length, charge, hydrophobicity

### Phase B-2 / B-2+ / B-2++: ProteinMPNN

| Mode | Backbone | Receptor context |
|---|---|---|
| B-2 | Ideal α-helix | None |
| B-2+ | Ideal α-helix | Pocket centroid |
| B-2++ | ESMFold predicted | Pocket centroid |

NLL scores from ProteinMPNN are converted to a [0, 1] range via sigmoid normalization (scale = 3.3).

### Phase B-1: AutoDock Vina

- **≤ 5 residues:** flexible docking (side-chain rotamers sampled)
- **≥ 6 residues:** rigid docking (backbone locked)
- **Units:** kcal/mol — more negative = stronger predicted binding

### Phase C-1 / C-2: Selectivity

```
C-1  selectivity_score     = rescoring(target)  − rescoring(off-target)
C-2  docking_selectivity   = docking(off-target) − docking(target)   [more positive = more selective]

selective_final_score = final_score + λ × selectivity_score   (default λ = 0.3)
```

- 🟢 selectivity ≥ 0.10 — target-selective
- 🟡 −0.05 < selectivity < 0.10 — moderate
- 🔴 selectivity ≤ −0.05 — off-target affinity risk

---

## Limitations

This is a proof-of-concept tool for hypothesis generation.

- Sequences are generated by biased random sampling, not structure-based design.
- Scores are relative rankings; they do not predict binding affinity in absolute terms.
- Physicochemical-based selectivity (C-1) has limited resolution for closely related protein families.
- All results require experimental validation.

---

## Tech Stack

`Python 3.11` · `Streamlit` · `BioPython` · `LightGBM` · `PyTorch (MPS)` · `ProteinMPNN` · `ESMFold (fair-esm)` · `AutoDock Vina` · `RDKit` · `py3Dmol` · `RCSB PDB API`

---

## License

MIT — see [LICENSE](LICENSE).

---

## Author

Developed as a portfolio project demonstrating multi-layer AI scoring for early-stage drug discovery.

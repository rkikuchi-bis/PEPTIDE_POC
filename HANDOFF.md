# Peptide Discovery App v2 — 引き継ぎドキュメント

作成日: 2026-04-13  
目的: peptide_poc から peptide_v2 への再設計に際し、新セッションの Claude に背景・方針・具体的な次のステップを伝える。

---

## 1. このプロジェクトの目的

社内研究者・製薬会社が「使ってみたい」と感じるレベルの、ペプチド創薬初期仮説生成ツールを作る。  
デモ品質から脱却し、**科学的精度を優先した設計**に転換する。

---

## 2. なぜ v2 を 1 から作るか

旧版（peptide_poc）の最大の問題は「**ランダム生成 → フィルタ**」という設計思想。  
研究者はこの部分を見た瞬間に信頼を失う。設計思想ごと変える必要があるため新規作成を選択した。

### 旧版で捨てる部分
- `core/generator.py` — ランダム配列生成
- `core/rescorer.py` — BioPython ProtParam ベースのスコアリング
- `core/pipeline.py` — 生成→フィルタのパイプライン全体
- `core/docking.py` — AutoDock Vina 剛体ドッキング
- `core/ml_scorer.py` — LightGBM バイナリ分類器
- `core/proteinmpnn.py` — スコアリング専用の ProteinMPNN 使用
- `core/pepfold.py` — ESMFold ローカル推論
- `core/selectivity.py` — Phase C-1/C-2 選択性スコア

### 旧版から移植する部分（そのままコピー可）
- `core/pdb_utils.py` — PDB/mmCIF パース・ポケット解析・centroid 計算
- `core/rcsb_client.py` — RCSB PDB REST API クライアント
- `ui/structure_viewer.py` — py3Dmol 3D可視化
- `ui/sidebar.py` の**構造入力・ポケット設定部分のみ**（スコアリング UI は不要）
- `data/example_targets/` — サンプル構造ファイル

---

## 3. 新しいアーキテクチャ（v2 の設計）

```
Input: PDB / mmCIF（ローカルアップロード または RCSB検索）
  ↓
ポケット解析（pdb_utils 流用）
  ↓
ProteinMPNN で配列生成（構造ベース逆設計）← ランダム生成を完全廃止
  ↓
Boltz-2 で複合体構造予測（ペプチド＋受容体）
  ↓
【二軸スコアリング】
  - iPSAE スコア（界面構造信頼度、Boltz-2 出力）
  - ΔG / Kd（PRODIGY、オープンソース熱力学スコアリング）
  ↓
二軸で候補を提示・可視化
```

### 科学的根拠
この設計は以下の論文にインスパイアされている:  
**"Single-Pass Discrete Diffusion Predicts High-Affinity Peptide Binders at >1,000 Sequences per Second across 150 Receptor Targets"**  
Andre Watson, Ligandal Inc., bioRxiv 2026.03.14.711748

論文の主要な主張:
- 構造ベースのペプチド設計において、iPSAE（構造信頼度）と ΔG（熱力学的有利性）は**直交する指標**であり、両方が必要
- iBSAE ≥ 0.5 を足切り閾値として、ΔG でランキングする戦略が最も有効
- BoltzGen（Boltz-2 + ProteinMPNN）は LigandForge（クローズドソース）とマッチしたサンプル数では同等の品質

つまり **Boltz-2（オープンソース）+ ProteinMPNN（オープンソース）+ PRODIGY（オープンソース）** で、論文の BoltzGen 相当のパイプラインが構築できる。

---

## 4. 技術スタック（v2）

### 開発環境
- MacBook Pro M4 MAX 128GB（MPS 利用可能）
- Python 3.11〜3.12
- パッケージ管理: `uv`（`pip` は使わない）
- 実行: `uv run streamlit run app.py`

### 主要パッケージ（インストール予定）
```
boltz          # Boltz-2 複合体構造予測
prodigy-prot   # ΔG / Kd 熱力学スコアリング
torch==2.11.0  # MPS 動作確認済み（旧環境より）
```

### 旧環境で確認済みのパッケージ（uv add 済み）
```
biopython, py3Dmol, streamlit, pandas, numpy, requests
rdkit, gemmi（必要であれば）
```

### 旧環境で未確認
```
boltz          # 要インストール・要動作確認
prodigy-prot   # 要インストール・要動作確認
```

---

## 5. v2 のコーディングルール

旧版から継承するルール:
- 新しいビジネスロジックは必ず `core/` に追加
- 新しいUI要素は必ず `ui/` に追加
- `app.py` は薄いオーケストレーター（imports・session_state・関数呼び出しのみ）
- `core/` 内のモジュールは Streamlit に依存しない
- 依存追加は `uv add <package>`

追加ルール（v2から）:
- スコアリングは必ず iPSAE + ΔG の二軸で行う（単一スコアでのランキングは避ける）
- 計算に時間がかかる処理（Boltz-2）には `st.spinner` と進捗表示を必ず付ける
- 科学的な制約・限界は UI に明示する（「これは計算上の予測です」）

---

## 6. 新セッションで最初にやること

### Step 1: 環境構築確認
```bash
uv run python -c "import torch; print(torch.__version__); print('MPS:', torch.backends.mps.is_available())"
```

### Step 2: Boltz-2 インストールと動作確認
```bash
uv add boltz
uv run python -c "import boltz; print(boltz.__version__)"
```

### Step 3: PRODIGY インストールと動作確認
```bash
uv add prodigy-prot
uv run python -c "from prodigy_prot import predict_IC; print('PRODIGY OK')"
```

### Step 4: 旧版から移植
以下のファイルを peptide_poc/ からコピーする:
```
peptide_poc/core/pdb_utils.py       → peptide_v2/core/pdb_utils.py
peptide_poc/core/rcsb_client.py     → peptide_v2/core/rcsb_client.py
peptide_poc/ui/structure_viewer.py  → peptide_v2/ui/structure_viewer.py
peptide_poc/data/example_targets/   → peptide_v2/data/example_targets/
```

### Step 5: 新規作成するファイル
```
peptide_v2/core/mpnn_generator.py   # ProteinMPNN で配列生成（構造ベース）
peptide_v2/core/boltz_predictor.py  # Boltz-2 複合体構造予測・iPSAE 計算
peptide_v2/core/prodigy_scorer.py   # PRODIGY ΔG / Kd スコアリング
peptide_v2/core/pipeline.py         # 新パイプライン
peptide_v2/ui/sidebar.py            # 構造入力・ポケット設定 UI
peptide_v2/ui/results.py            # 二軸スコア表示 UI
peptide_v2/app.py                   # 薄いオーケストレーター
```

---

## 7. 優先度の判断基準

機能を追加する際の優先基準（この順番で判断する）:

1. **科学的妥当性** — 研究者が「根拠がある」と思えるか
2. **説明可能性** — なぜこのペプチドが良いか説明できるか
3. **使いやすさ** — 構造を入れて Run するだけで動くか
4. **スループット** — M4 Mac で現実的な待ち時間か

---

## 8. ユーザーについて

- ライフサイエンス企業でのAI/ML業務委託を目指している
- 控えめな表現を好む
- 製薬会社・社内研究者への提供を視野に入れているが、すぐ売り込む意図はない
- 科学的精度を最優先することを明確に決断している

---

## 9. 参照すべき旧版ファイル

詳細な実装の参考として手元に残しておく旧版ファイル:
```
peptide_poc/CLAUDE.md       # 旧版アーキテクチャ・コーディングルール
peptide_poc/CONTEXT.md      # 旧版の実装状態（Phase A〜C の詳細）
peptide_poc/core/pdb_utils.py
peptide_poc/core/rcsb_client.py
peptide_poc/ui/structure_viewer.py
```

---

以上。新セッションでは Step 1 の環境確認から始めてください。

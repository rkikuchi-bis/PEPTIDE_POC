# Current State
最終更新: 2026-04-08（Phase B-2+ 完了: 受容体条件付き ProteinMPNN）

---

## Completed

### モジュール分割（リファクタリング）
- `app.py` を847行のモノリスから約60行の薄いオーケストレーターに分割
- `core/pipeline.py` — 候補生成パイプラインのオーケストレーター
- `core/rcsb_client.py` — HTTP通信・RCSB API・`InMemoryUploadedStructure`
- `core/structure_scorer.py` — 構造優先度スコアリング関数
- `ui/sidebar.py` — `render_sidebar() -> dict`
- `ui/results.py` — `render_results(result_df, pdb_summary)`
- `ui/structure_viewer.py` — py3Dmol 3D可視化

### 3D構造可視化（py3Dmol）
- 依存追加: `py3Dmol>=2.0.0`（`uv add` 済み）
- 表示内容: チェーン（ブルー）、ポケット（オレンジ）、リガンド（グリーン）
- `ui/structure_viewer.py` → `render_viewer_section()` として `app.py` から呼び出し
- 注意: `_repr_html_()` はJupyter専用のため `_make_html()` を使用

### Phase A-1: BioPython ProtParam スコアリング
- `core/rescorer.py` を全面刷新
- 旧: 独自アミノ酸分率によるヒューリスティック
- 新: BioPython ProteinAnalysis による科学的特徴量
  - `isoelectric_point` (pI) → 電荷マッチスコア（生理的pH 7.4基準）
  - `gravy` (GRAVY) → 疎水性マッチスコア（シグモイド正規化）
  - `instability_index` → 安定性スコア（Guruprasad et al., 1990）
  - `aromaticity`, `helix_fraction`, 配列多様性 → 複雑性スコア
  - 新規出力カラム: `gravy`, `instability_index`, `isoelectric_point`,
    `aromaticity`, `helix_fraction`, `turn_fraction`, `sheet_fraction`, `molecular_weight`

### 環境移行
- 開発環境: Windows → MacBook Pro M4 MAX 128GB に完全移行
- GitHub push 済み: `https://github.com/rkikuchi-bis/PEPTIDE_POC`
- `.gitignore` に `outputs/`, `.claude/` を追加

---

### Phase A-2: ML スコアリング
- `scripts/prepare_dataset.py` — RCSB PDB REST API（v2）からペプチド陽性例485件取得、ランダム陰性例500件生成
  - `data/peptide_dataset.csv` に保存（特徴量: length, net_charge, avg_hydrophobicity + Phase A-1 全特徴量）
- `scripts/train_classifier.py` — RandomForest / LightGBM の 5-fold CV 比較
  - LightGBM が最良（CV AUC=0.891、Holdout AUC=0.881、Accuracy=81%）
  - `models/peptide_classifier.joblib` に保存
- `core/ml_scorer.py` — `score_with_ml(sequence) -> float`（joblib ロード・キャッシュ付き）
- `core/rescorer.py` — `ml_score` を統合、`final_score` の重みを自動切替
  - ML あり: gen×0.20 + property×0.20 + rescoring×0.30 + ml×0.30
  - ML なし: gen×0.30 + property×0.30 + rescoring×0.40（後方互換）
- `ui/results.py` — テーブルと詳細パネルに `ml_score` / `ML binding probability` を追加

---

## Issues

- なし（現時点で既知のバグなし）

---

### Phase B-1: AutoDock Vina ドッキング（実装中）
- `bin/vina` — mac_aarch64 バイナリをダウンロード済み（Python binding は M4 非対応）
- 追加依存: `rdkit`, `meeko>=0.7.1`, `gemmi`
- `core/docking.py` — ペプチド → PDBQT → Vina 実行 → スコア返却
  - `_sequence_to_pdbqt()` — RDKit `MolFromSequence` + meeko `PDBQTWriterLegacy`
  - `prepare_receptor_pdbqt()` — 受容体 PDB → PDBQT（キャッシュ付き）
  - `dock_top_candidates()` — 上位 N 候補に一括ドッキング
- `ui/results.py` — `docking_score` カラム・詳細パネルに追加
- 依存ツール: `obabel`（Homebrew）、`rdkit`、`gemmi`
- 受容体変換: HETATM 除去 → `obabel -xr`（剛体 PDBQT）
- リガンド変換: RDKit 3D生成 → `obabel -xr` → ROOT/ENDROOT/TORSDOF 0 付加
- ドッキング: Vina バイナリ（subprocess）→ スコア [kcal/mol] パース
- **制約**: 7残基以上は剛体ドッキング（スコアが正値になることあり）。相対ランキングには利用可能。
- エンドツーエンドテスト完了（1HSG + 5候補）

---

### Phase B-2: ProteinMPNN スコアリング（実装完了）
- `torch==2.11.0` インストール済み（MPS 動作確認済み）
- `models/proteinmpnn/v_48_020.pt` — vanilla_model_weights (6.4MB) ダウンロード済み
- `models/proteinmpnn/protein_mpnn_utils.py` — モデルクラス
- `scripts/mpnn_scorer.py` — スタンドアロンスコアラー（LightGBM との OpenMP 競合を回避するためサブプロセスとして実行）
- `core/proteinmpnn.py` — `score_sequences_batch()`, `score_with_proteinmpnn()` API
- `core/rescorer.py` — `proteinmpnn_score` を統合、`final_score` の重み自動切替
  - ML + MPNN あり: gen×0.15 + property×0.15 + rescoring×0.25 + ml×0.25 + mpnn×0.20
  - ML あり: gen×0.20 + property×0.20 + rescoring×0.30 + ml×0.30
  - なし: gen×0.30 + property×0.30 + rescoring×0.40
- `ui/results.py` — `proteinmpnn_score` カラム・詳細パネルに追加
- **設計**: 理想αヘリックス骨格（rise=1.5Å, radius=2.3Å, rotation=100°/residue）の対数尤度スコア
- **制約**: 構造フリーのヒューリスティック。相対ランキング用途（絶対値の科学的解釈は限定的）
- **OpenMP 競合**: LightGBM + PyTorch の競合は mpnn_scorer.py をサブプロセス実行することで解決

---

### Phase B-2+: 受容体条件付き ProteinMPNN スコアリング（実装完了）
- `scripts/mpnn_scorer_receptor.py` — 受容体バックボーン + ペプチド理想ヘリックスを結合
  - 受容体 (chain_M=0): グラフ上の固定コンテキスト（スコア対象外）
  - ペプチド (chain_M=1): pocket centroid に配置、NLL をスコアに変換
  - 受容体が大きい場合（>200残基）は centroid 近傍の上位 200 残基のみ使用
- `core/pdb_utils.py` — `get_pocket_ca_centroid(structure, pdb_summary)` を追加
  - Manual region / Ligand neighborhood 両モードに対応
- `core/proteinmpnn.py` — `score_sequences_with_receptor()` を追加
  - mmCIF は PDB 変換して渡す（BioPython PDBIO 経由）
- `core/rescorer.py` — `rescore_candidates()` に `structure_text`, `file_format`, `pocket_centroid` 引数を追加
  - 3 引数が揃えば受容体条件付きモード、揃わなければ構造フリーモード（B-2）へ自動フォールバック
  - `proteinmpnn_receptor_conditioned` カラムでどちらのモードかを記録
- `core/pipeline.py` — 上記 3 引数を `rescore_candidates()` に伝搬
- `app.py` — Run ボタン時に構造パース + `get_pocket_ca_centroid()` を呼び出し、`run_pipeline()` に渡す

## Next Steps

### Phase B-1+ 改善案（任意）
- 短鎖（≤5残基）は柔軟ドッキング・長鎖は剛体の自動切り替え
- obabel 依存削減（meeko が arm64 対応次第）

### Phase B-2++ 改善案（任意）
- PepFold による初期骨格生成 → ProteinMPNN への入力に使用（ダミーヘリックスより精度向上）
- AlphaFold2-Multimer による結合構造再予測

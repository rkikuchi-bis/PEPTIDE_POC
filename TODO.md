# TODO
最終更新: 2026-04-09（Phase B-2++ 動作確認完了）

---

## 完了済み短期タスク（Phase A-2）

- [x] RCSB PDB REST API からペプチド複合体データを取得・前処理する（Option B 実装済み）
      → `scripts/prepare_dataset.py` で陽性485件・陰性500件取得
      → `data/peptide_dataset.csv` に保存

- [ ] 【将来検討】PepBDB 全件ダウンロードによる高精度学習（Option A）
      URL: http://huanglab.phys.hust.edu.cn/pepbdb/db/download/pepbdb-20200318.tgz
      サイズ: 659MB、PDB構造ファイルから配列抽出が必要
      → Option B のモデル精度が不十分な場合に採用を検討する

- [x] `scripts/train_classifier.py` を作成（LightGBM、CV AUC=0.891）
- [x] `core/ml_scorer.py` を新設（`score_with_ml` 実装）
- [x] `core/rescorer.py` に ml_score を統合
- [x] `ui/results.py` に ml_score カラムを追加表示

---

## 中期タスク（Phase B-1）

- [x] `vina` バイナリ (mac_aarch64) を `bin/vina` にダウンロード（Python bindingはM4非対応）
- [x] `gemmi`, `rdkit`, `meeko` のインストール
- [x] `core/docking.py` を新設
      - `is_vina_available()` — バイナリ存在チェック
      - `_sequence_to_pdbqt()` — RDKit + meeko でペプチド → PDBQT
      - `prepare_receptor_pdbqt()` — PDB → PDBQT（キャッシュ付き）
      - `dock_peptide()` — 単一ペプチドドッキング
      - `dock_top_candidates()` — 上位 N 候補に一括実行
- [x] `ui/results.py` に `docking_score` カラムを追加
- [x] `ui/sidebar.py` にドッキング設定UI（ボックス中心自動計算・サイズ・top_n・exhaustiveness）を追加
- [x] `app.py` に "Run Docking" ボタンを追加（pipeline とは独立したオプション実行）
- [x] 受容体 PDBQT 準備の動作検証（1HSG で確認済み）
      - obabel による受容体 PDBQT 変換（HETATM除去 → -xr オプション）
      - obabel によるリガンド PDBQT 変換（剛体: ROOT/ENDROOT/TORSDOF 0 付加）
      - dock_top_candidates のエンドツーエンドテスト完了
      - 注意: 7残基以上は剛体ドッキング（相対ランキング用途、絶対値は不正確）
             4残基以下は柔軟ドッキングで負のスコアが得られる

## Phase B-1 完了 ✅

## Phase B-1+ 完了 ✅

- [x] 短鎖ペプチド（≤5残基）のみ柔軟ドッキング、長鎖は剛体の自動切り替え
      → `core/docking.py`: `FLEXIBLE_DOCKING_MAX_LENGTH=5`、`dock_peptide(flexible=None)` で自動判定
      → `docking_mode` カラム追加、詳細パネルに flexible/rigid を表示
- [ ] obabel 依存を減らす（meeko が将来 arm64 対応した場合に置き換え）

---

## Phase B-2 完了 ✅

## Phase B-2+ 完了 ✅

- [x] 受容体構造を条件付けとした ProteinMPNN スコアリング
      （受容体 PDB + ペプチドダミー座標を結合 → より科学的なスコア）
      → `scripts/mpnn_scorer_receptor.py`, `core/pdb_utils.get_pocket_ca_centroid()`,
         `core/proteinmpnn.score_sequences_with_receptor()` 実装済み
      → 構造未指定時は B-2（構造フリー）に自動フォールバック

## Phase B-2++ 完了 ✅

- [x] ESMFold ローカル推論（fair-esm）によるペプチド骨格生成
      → `core/pepfold.py`、`scripts/esmfold_scorer.py` 実装済み
      → 設計: 全候補はヘリックスで高速処理 → diversity 後の少数候補に ESMFold で再スコアリング
      → 予測失敗時は理想ヘリックスへ自動フォールバック（後方互換）
      → 環境変数: PEPFOLD_MAX_SEQS（デフォルト30）
      → 動作確認済み: Rescoring notes に MPNN(ESMFold+receptor) が表示される

## Direction A + B 完了 ✅（2026-04-09）

- [x] Direction A: Explainability — `core/explainer.py`, `ui/results.py` 更新
      → 各候補の推薦理由を日本語自然言語で説明（pI/GRAVY/安定性/ML/MPNN/ドッキング/選択性を統合）
- [x] Direction B: Selectivity Phase C-1 — `core/selectivity.py`, `ui/sidebar.py`, `app.py` 更新
      → ターゲットとオフターゲットの rescoring 差分で選択性スコアを定量化
      → テーブルの選択性ソート、詳細パネルの 🟢/🟡/🔴 表示

## λ デフォルト化 完了 ✅（2026-04-09）

- [x] λ=0.3 をデフォルト値に設定（`ui/sidebar.py`）
- [x] サイドバーに「推奨: 0.3」ラベルと根拠をツールチップに記載

## UI リファクタリング（完了 2026-04-10）

- [x] Simple / Expert mode 分割（`ui/sidebar.py` 全面刷新）
- [x] Target label → Project name（空白可、CSV ファイル名用）
- [x] Chain 自動推奨（`get_recommended_chain()` 追加）
- [x] Selectivity を Expert mode 専用に移動
- [x] Off-target DB 実装（`core/offtarget_db.py`、14ターゲット収録）
- [x] Known off-target DB から選択 → RCSB 自動DL フロー
- [x] オフターゲット構造アップロード → pocket centroid 自動計算

## デモ品質改善（完了 2026-04-10）

- [x] Bug fix: `default_offtarget_hydrophobicity` 未定義変数（Expert Selectivity 使用時にクラッシュ）
- [x] Simple mode: 構造未ロード時に Run ボタンを disabled + ガイドメッセージ
- [x] Simple mode: Run ボタンラベルを「▶ Run」に変更
- [x] app.py: タイトルに emoji 追加、キャプション日英化
- [x] app.py: 生成完了メッセージを「完了 ✅ N件」形式に簡潔化
- [x] results.py: 初期状態メッセージを日英化・改善
- [x] results.py: 表示カラムをデモ向けに再整理（主要スコアを先頭、詳細を後方）
- [x] results.py: 推薦理由ラベルを短縮

## 将来タスク（優先度低）

- [ ] Phase C-2: ドッキングスコア差分による選択性（選択性の根本的解決）
- [ ] 3D Viewer ペプチド重畳: Option 2（ESMFold骨格）・Option 3（Vinaドッキングポーズ保存）
- [ ] AlphaFold2-Multimer による結合構造再予測（B-3: 長期・高コスト）
- [ ] obabel 依存削減（meeko が arm64 対応次第）

---

## 完了済み ✅

- [x] app.py のモジュール分割（ui/ core/ に整理）
- [x] py3Dmol による3D構造可視化
- [x] Phase A-1: BioPython ProtParam スコアリング
- [x] GitHub push・Mac環境への移行
- [x] CLAUDE.md / CONTEXT.md / TODO.md / prompts/ の整備
- [x] Phase A-2: ML スコアリング（LightGBM、AUC=0.891）

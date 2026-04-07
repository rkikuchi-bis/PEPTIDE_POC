# TODO
最終更新: 2026-04-08（Phase B-1 完了）

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

## 今後の改善案（Phase B-1+）

- [ ] 短鎖ペプチド（≤5残基）のみ柔軟ドッキング、長鎖は剛体の自動切り替え
- [ ] obabel 依存を減らす（meeko が将来 arm64 対応した場合に置き換え）

---

## 長期タスク（Phase B-2）

- [ ] ProteinMPNN の M4 MPS 対応確認
      → `device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")`
- [ ] PepFold による初期骨格生成パイプライン
- [ ] ProteinMPNN → 配列確率分布 → スコア変換
- [ ] AlphaFold2-Multimer による結合構造再予測

---

## 完了済み ✅

- [x] app.py のモジュール分割（ui/ core/ に整理）
- [x] py3Dmol による3D構造可視化
- [x] Phase A-1: BioPython ProtParam スコアリング
- [x] GitHub push・Mac環境への移行
- [x] CLAUDE.md / CONTEXT.md / TODO.md / prompts/ の整備
- [x] Phase A-2: ML スコアリング（LightGBM、AUC=0.891）

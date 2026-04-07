# TODO
最終更新: 2026-04-08（Phase A-2 完了）

---

## 短期タスク（Phase A-2）

- [ ] RCSB PDB REST API からペプチド複合体データを取得・前処理する（Option B）
      → ペプチド含む複合体エントリを検索（配列のみ取得、構造不要）
      → 結合ペプチド配列（陽性例）を抽出
      → ランダム生成配列（陰性例）を作成
      → Phase A-1 の特徴量で特徴ベクトル化

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

- [ ] `vina` (Python バインディング) の動作確認（Mac M4対応）
- [ ] `meeko` によるペプチド→PDBQT変換の実装
- [ ] `core/docking.py` を新設
- [ ] 上位10〜20候補のみドッキング実行するロジック

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

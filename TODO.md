# TODO
最終更新: 2026-04-08

---

## 短期タスク（Phase A-2）

- [ ] PepBDB データセットをダウンロード・前処理する
      URL: http://huanglab.phys.hust.edu.cn/pepbdb/
      → 結合ペプチド配列（陽性例）を抽出
      → ランダム生成配列（陰性例）を作成
      → Phase A-1 の特徴量で特徴ベクトル化

- [ ] `scripts/train_classifier.py` を作成（学習スクリプト）
      → RandomForest または LightGBM で二値分類
      → 評価: AUC-ROC, Precision/Recall
      → モデルを `models/peptide_classifier.joblib` に保存

- [ ] `core/ml_scorer.py` を新設
      → `score_with_ml(sequence: str) -> float` を実装
      → モデルをロードしてスコアを返す

- [ ] `core/rescorer.py` に ml_score を統合
      → `final_score` の重みを調整（ml_score追加分）

- [ ] `ui/results.py` に ml_score カラムを追加表示

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

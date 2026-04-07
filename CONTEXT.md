# Current State
最終更新: 2026-04-08

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

## Issues

- なし（現時点で既知のバグなし）

---

## Next Steps

### 直近（Phase A-2）
1. PepBDB からタンパク質結合ペプチドデータを取得
2. Phase A-1の特徴量で特徴ベクトル化
3. RandomForest/LightGBM で二値分類モデルを学習
4. `core/ml_scorer.py` を新設、`core/pipeline.py` に組み込み
5. 学習済みモデルを `models/peptide_classifier.joblib` として保存

### その後（Phase B-1）
- AutoDock Vina Python バインディングによるドッキングスコア追加
- 上位候補のみ（10〜20件）に実行

### 将来（Phase B-2）
- ProteinMPNN による配列逆設計（M4 MPS対応）
- `device = torch.device("mps")` で動作確認済みのコードを書く予定

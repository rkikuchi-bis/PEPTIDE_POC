"""
既知ターゲット / Off-target ペアのキュレーションDB。

構造:
    OFFTARGET_DB[target_key] = {
        "target_label": str,
        "target_pdb": str,          # 代表的なターゲット PDB ID（参考表示のみ）
        "offtargets": [
            {
                "label":       str,   # 表示名
                "pdb_id":      str,   # RCSB からダウンロードする PDB ID
                "description": str,   # 交差反応の科学的背景（日英）
            },
            ...
        ]
    }

チェーン選択は get_recommended_chain()（最多残基チェーン自動選択）に委ねる。
"""

OFFTARGET_DB: dict[str, dict] = {

    # ── Kinases ─────────────────────────────────────────────────────
    "EGFR": {
        "target_label": "EGFR",
        "target_pdb": "2J6M",
        "offtargets": [
            {
                "label": "HER2 (3PP0)",
                "pdb_id": "3PP0",
                "description": (
                    "EGFRファミリーメンバー。キナーゼドメインが高相同性のため、"
                    "EGFR阻害薬が HER2 に交差結合しやすい。"
                    "/ ErbB2 kinase — high homology to EGFR; common EGFR inhibitors show cross-reactivity."
                ),
            },
            {
                "label": "SRC (2SRC)",
                "pdb_id": "2SRC",
                "description": (
                    "非受容体型チロシンキナーゼ。ATP結合部位の形状が類似しており、"
                    "広域スペクトルキナーゼ阻害薬で問題となる。"
                    "/ Non-receptor tyrosine kinase; similar ATP-binding pocket shape causes cross-inhibition."
                ),
            },
            {
                "label": "MET (1R0P)",
                "pdb_id": "1R0P",
                "description": (
                    "受容体型チロシンキナーゼ。EGFR耐性後に MET 増幅が起こることが多く、"
                    "選択性評価が重要。"
                    "/ RTK amplified as EGFR-resistance mechanism; selectivity is clinically important."
                ),
            },
        ],
    },

    "BCR-ABL / ABL1": {
        "target_label": "BCR-ABL / ABL1",
        "target_pdb": "1IEP",
        "offtargets": [
            {
                "label": "KIT (1T46)",
                "pdb_id": "1T46",
                "description": (
                    "イマチニブは KIT にも結合し消化管間質腫瘍(GIST)に有効だが、"
                    "副作用の原因にもなる。"
                    "/ Imatinib binds KIT; beneficial for GIST but also causes off-target toxicity."
                ),
            },
            {
                "label": "PDGFR-β (3MJG)",
                "pdb_id": "3MJG",
                "description": (
                    "血小板由来増殖因子受容体。ABL阻害薬が交差結合し、浮腫などの副作用を引き起こす。"
                    "/ PDGF receptor — cross-binding by ABL inhibitors causes edema as side effect."
                ),
            },
        ],
    },

    "ALK": {
        "target_label": "ALK",
        "target_pdb": "2XP2",
        "offtargets": [
            {
                "label": "ROS1 (3ZBF)",
                "pdb_id": "3ZBF",
                "description": (
                    "ALK と ROS1 はキナーゼドメインが類似。クリゾチニブは両方を阻害する。"
                    "/ ALK and ROS1 share kinase domain similarity; crizotinib inhibits both."
                ),
            },
            {
                "label": "IGF1R (2OJ9)",
                "pdb_id": "2OJ9",
                "description": (
                    "インスリン様増殖因子受容体。ALK阻害薬との交差反応で血糖調節に影響することがある。"
                    "/ Insulin-like growth factor receptor; cross-reactivity may affect glucose regulation."
                ),
            },
        ],
    },

    # ── GPCRs ────────────────────────────────────────────────────────
    "β2-adrenergic receptor": {
        "target_label": "β2-adrenergic receptor",
        "target_pdb": "2RH1",
        "offtargets": [
            {
                "label": "β1-adrenergic receptor (2VT4)",
                "pdb_id": "2VT4",
                "description": (
                    "β2選択的薬剤（喘息治療）がβ1（心臓）に交差すると頻脈・心悸亢進の副作用が生じる。"
                    "/ β2-selective drugs (asthma) crossing to cardiac β1 cause tachycardia."
                ),
            },
        ],
    },

    "Dopamine D2 receptor": {
        "target_label": "Dopamine D2 receptor",
        "target_pdb": "6CM4",
        "offtargets": [
            {
                "label": "Serotonin 5-HT2A (6A93)",
                "pdb_id": "6A93",
                "description": (
                    "多くの抗精神病薬が D2 と 5-HT2A を同時に阻害する。"
                    "5-HT2A結合が錐体外路症状を緩和するため臨床的に利用される場合もある。"
                    "/ Many antipsychotics co-inhibit D2 and 5-HT2A; sometimes exploited clinically."
                ),
            },
        ],
    },

    "Histamine H1 receptor": {
        "target_label": "Histamine H1 receptor",
        "target_pdb": "3RZE",
        "offtargets": [
            {
                "label": "Muscarinic M1 (5CXV)",
                "pdb_id": "5CXV",
                "description": (
                    "古典的抗ヒスタミン薬はムスカリン受容体にも結合し、口渇・便秘などの副作用を引き起こす。"
                    "/ Classical antihistamines bind muscarinic receptors causing dry mouth and constipation."
                ),
            },
        ],
    },

    # ── Nuclear receptors ─────────────────────────────────────────────
    "Estrogen receptor": {
        "target_label": "Estrogen receptor",
        "target_pdb": "1ERE",
        "offtargets": [
            {
                "label": "Androgen receptor (1T7R)",
                "pdb_id": "1T7R",
                "description": (
                    "ステロイドホルモン受容体間でLBDの類似性があり、選択的モジュレーターの交差活性が問題になる。"
                    "/ Steroid hormone receptor LBD similarity causes cross-activity of selective modulators."
                ),
            },
            {
                "label": "Progesterone receptor (1A28)",
                "pdb_id": "1A28",
                "description": (
                    "エストロゲン受容体作動薬が黄体ホルモン受容体に部分作動し、副作用を生じることがある。"
                    "/ ER agonists may partially activate progesterone receptor causing off-target effects."
                ),
            },
        ],
    },

    "Androgen receptor": {
        "target_label": "Androgen receptor",
        "target_pdb": "1T7R",
        "offtargets": [
            {
                "label": "Estrogen receptor (1ERE)",
                "pdb_id": "1ERE",
                "description": (
                    "AR拮抗薬がERを部分活性化することで女性化乳房などの副作用を引き起こすことがある。"
                    "/ AR antagonists may partially activate ER, causing gynecomastia."
                ),
            },
            {
                "label": "Glucocorticoid receptor (1M2Z)",
                "pdb_id": "1M2Z",
                "description": (
                    "ARとGRはLBDの類似性が高く、AR薬剤のGR交差結合が代謝副作用を生じさせる。"
                    "/ High LBD similarity between AR and GR; cross-binding causes metabolic side effects."
                ),
            },
        ],
    },

    "PPARγ": {
        "target_label": "PPARγ",
        "target_pdb": "2PRG",
        "offtargets": [
            {
                "label": "PPARα (2P54)",
                "pdb_id": "2P54",
                "description": (
                    "PPARファミリーはLBDが類似。PPARγ選択的薬剤がPPARαを活性化すると脂質代謝に影響する。"
                    "/ PPAR family LBD similarity; PPARγ drugs activating PPARα affect lipid metabolism."
                ),
            },
            {
                "label": "PPARδ (2ZNP)",
                "pdb_id": "2ZNP",
                "description": (
                    "PPARδ活性化は骨格筋エネルギー代謝に関与し、選択性の欠如が筋肉副作用を生じる。"
                    "/ PPARδ activation affects skeletal muscle energy metabolism; selectivity loss causes myopathy."
                ),
            },
        ],
    },

    # ── Ion channels ──────────────────────────────────────────────────
    "hERG channel": {
        "target_label": "hERG channel",
        "target_pdb": "5VA1",
        "offtargets": [
            {
                "label": "Nav1.5 (6LQA)",
                "pdb_id": "6LQA",
                "description": (
                    "hERG阻害（QT延長）と Nav1.5阻害（伝導障害）はどちらも心毒性の主要メカニズム。"
                    "/ hERG block (QT prolongation) and Nav1.5 block (conduction failure) are major cardiotoxicity mechanisms."
                ),
            },
        ],
    },

    # ── Proteases ────────────────────────────────────────────────────
    "HIV Protease": {
        "target_label": "HIV Protease",
        "target_pdb": "1HSG",
        "offtargets": [
            {
                "label": "Cathepsin B (1CS8)",
                "pdb_id": "1CS8",
                "description": (
                    "ヒトリソソームシステインプロテアーゼ。HIV プロテアーゼ阻害薬との構造的類似で交差阻害が報告されている。"
                    "/ Human lysosomal cysteine protease; structural similarity with HIV protease causes cross-inhibition."
                ),
            },
            {
                "label": "Renin (2V0Z)",
                "pdb_id": "2V0Z",
                "description": (
                    "ヒトアスパルテートプロテアーゼ。HIV プロテアーゼ阻害薬が血圧調節に影響するリスクがある。"
                    "/ Human aspartate protease; HIV protease inhibitors risk affecting blood pressure regulation."
                ),
            },
        ],
    },

    "Thrombin": {
        "target_label": "Thrombin",
        "target_pdb": "1PPB",
        "offtargets": [
            {
                "label": "Factor Xa (1FJS)",
                "pdb_id": "1FJS",
                "description": (
                    "トロンビンと Factor Xa はどちらも血液凝固カスケードのセリンプロテアーゼ。"
                    "トロンビン阻害薬が Factor Xa に交差すると出血リスクが増大する。"
                    "/ Both are serine proteases in the coagulation cascade; cross-inhibition increases bleeding risk."
                ),
            },
            {
                "label": "Trypsin (1TRN)",
                "pdb_id": "1TRN",
                "description": (
                    "消化酵素のセリンプロテアーゼ。トロンビン阻害薬の trypsin 交差結合は消化器副作用の原因になりうる。"
                    "/ Digestive serine protease; cross-binding of thrombin inhibitors may cause GI side effects."
                ),
            },
        ],
    },

    # ── Epigenetics ──────────────────────────────────────────────────
    "BRD4": {
        "target_label": "BRD4",
        "target_pdb": "4HBV",
        "offtargets": [
            {
                "label": "BRD2 (2DVQ)",
                "pdb_id": "2DVQ",
                "description": (
                    "BET ブロモドメインファミリー。JQ1 などの汎 BET 阻害薬は BRD2/3/4 を同時に阻害する。"
                    "/ BET bromodomain family; pan-BET inhibitors (JQ1 etc.) co-inhibit BRD2/3/4."
                ),
            },
            {
                "label": "BRD3 (4BJX)",
                "pdb_id": "4BJX",
                "description": (
                    "BRD3 阻害は転写調節に関与し、造血器毒性のリスクがある。"
                    "/ BRD3 inhibition affects transcriptional regulation; risk of hematological toxicity."
                ),
            },
        ],
    },

    # ── CYP ──────────────────────────────────────────────────────────
    "CYP3A4": {
        "target_label": "CYP3A4",
        "target_pdb": "1TQN",
        "offtargets": [
            {
                "label": "CYP2D6 (2F9Q)",
                "pdb_id": "2F9Q",
                "description": (
                    "主要薬物代謝酵素。CYP3A4 阻害薬が CYP2D6 にも作用すると、"
                    "他剤の代謝が低下し薬物相互作用が生じる。"
                    "/ Major drug-metabolizing enzyme; dual CYP3A4/2D6 inhibition causes drug-drug interactions."
                ),
            },
        ],
    },
}


def get_offtarget_options(target_key: str) -> list[dict]:
    """ターゲットキーに対応する Off-target リストを返す。存在しない場合は空リスト。"""
    entry = OFFTARGET_DB.get(target_key)
    if entry is None:
        return []
    return entry["offtargets"]


def list_target_keys() -> list[str]:
    """登録済みターゲットキーの一覧を返す。"""
    return list(OFFTARGET_DB.keys())

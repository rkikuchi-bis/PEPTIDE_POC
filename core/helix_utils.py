"""
理想 αヘリックス骨格座標生成ユーティリティ

scripts/mpnn_scorer_receptor.py の _ideal_helix_at_centroid を流用し、
3D Viewer 用の PDB 文字列も生成する。
"""
import math
import numpy as np

AA1_TO_AA3 = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}

# PDB 標準の原子名フォーマット（4文字、左寄せ）
_ATOM_NAME_PDB = {"N": " N  ", "CA": " CA ", "C": " C  ", "O": " O  "}
_ATOM_NAMES = ["N", "CA", "C", "O"]
_ELEMENTS   = ["N", "C",  "C", "O"]


def ideal_helix_coords(length: int, centroid) -> np.ndarray:
    """
    pocket centroid を中心に配置した理想 αヘリックス座標 (length, 4, 3) を返す。
    4 原子: N, CA, C, O。ヘリックス軸は z 方向。
    """
    RISE = 1.5
    RADIUS = 2.3
    DELTA = math.radians(100.0)

    centroid = np.asarray(centroid, dtype=np.float32)
    coords = np.zeros((length, 4, 3), dtype=np.float32)
    z_offset = float(centroid[2]) - (length - 1) * RISE / 2.0

    for i in range(length):
        angle = i * DELTA
        ca = np.array([
            float(centroid[0]) + RADIUS * math.cos(angle),
            float(centroid[1]) + RADIUS * math.sin(angle),
            z_offset + i * RISE,
        ], dtype=np.float32)

        tangent = np.array([
            -math.sin(angle), math.cos(angle),
            RISE / math.sqrt(RADIUS ** 2 + RISE ** 2),
        ], dtype=np.float32)
        tangent /= np.linalg.norm(tangent)
        radial = np.array([math.cos(angle), math.sin(angle), 0.0], dtype=np.float32)

        coords[i, 0] = ca - tangent * 1.46                       # N
        coords[i, 1] = ca                                         # CA
        coords[i, 2] = ca + tangent * 1.52                       # C
        coords[i, 3] = (ca + tangent * 1.52) + radial * 1.24    # O

    return coords


def helix_coords_to_pdb(sequence: str, centroid, chain_id: str = "P") -> str:
    """
    ペプチド配列と pocket centroid から、理想ヘリックス骨格の PDB 文字列を返す。

    chain_id: 受容体と区別するため "P"（Peptide）を使用（デフォルト）
    """
    coords = ideal_helix_coords(len(sequence), centroid)
    lines: list[str] = []
    serial = 1

    for res_i, aa1 in enumerate(sequence.upper()):
        aa3 = AA1_TO_AA3.get(aa1, "ALA")
        res_seq = res_i + 1
        atom_coords = coords[res_i]  # shape (4, 3)

        for atom_j, atom_name in enumerate(_ATOM_NAMES):
            x, y, z = float(atom_coords[atom_j, 0]), float(atom_coords[atom_j, 1]), float(atom_coords[atom_j, 2])
            element = _ELEMENTS[atom_j]
            pdb_atom = _ATOM_NAME_PDB[atom_name]  # 4 chars

            line = (
                "ATOM  "
                f"{serial:5d}"
                " "               # col 12
                f"{pdb_atom}"     # cols 13-16 (atom name)
                " "               # alt loc col 17
                f"{aa3:3s}"       # residue name cols 18-20
                " "               # col 21
                f"{chain_id}"     # chain col 22
                f"{res_seq:4d}"   # resseq cols 23-26
                " "               # insertion code col 27
                "   "             # cols 28-30
                f"{x:8.3f}"       # x cols 31-38
                f"{y:8.3f}"       # y cols 39-46
                f"{z:8.3f}"       # z cols 47-54
                "  1.00"          # occupancy cols 55-60
                "  0.00"          # b-factor cols 61-66
                "          "      # cols 67-76
                f"{element:>2s}"  # element cols 77-78
                "  "
            )
            lines.append(line)
            serial += 1

    # CONECTレコード: 残基内 N-CA, CA-C, C-O と残基間ペプチド結合 C(i)-N(i+1)
    # serial番号: 残基 i の N=4i+1, CA=4i+2, C=4i+3, O=4i+4 (1-indexed)
    n_res = len(sequence)
    for i in range(n_res):
        base = i * 4 + 1  # N のシリアル番号
        n_ser, ca_ser, c_ser, o_ser = base, base + 1, base + 2, base + 3
        lines.append(f"CONECT{n_ser:5d}{ca_ser:5d}")   # N-CA
        lines.append(f"CONECT{ca_ser:5d}{n_ser:5d}{c_ser:5d}")  # CA-N, CA-C
        lines.append(f"CONECT{c_ser:5d}{ca_ser:5d}{o_ser:5d}")  # C-CA, C-O
        lines.append(f"CONECT{o_ser:5d}{c_ser:5d}")   # O-C
        # ペプチド結合: C(i) - N(i+1)
        if i < n_res - 1:
            next_n_ser = (i + 1) * 4 + 1
            lines.append(f"CONECT{c_ser:5d}{next_n_ser:5d}")
            lines.append(f"CONECT{next_n_ser:5d}{c_ser:5d}")

    lines.append("END")
    return "\n".join(lines)

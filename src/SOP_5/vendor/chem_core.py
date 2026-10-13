"""化学核心：SDF 解析、身份核验、描述符、类药性口径、化学型、相互作用指纹。

口径冻结（见 config/frozen_rule.json）：
  描述符  MW = Descriptors.MolWt（不圆整）        cLogP = Crippen.MolLogP
          HBD/HBA = Lipinski.NumHDonors/NumHAcceptors
          TPSA = rdMolDescriptors.CalcTPSA     RotB = Lipinski.NumRotatableBonds
  Ro5 违规  : MW > 500, cLogP > 5, HBD > 5, HBA > 10（严格大于；不先圆整）
  Veber 违规: RotB > 10, TPSA > 140
  PAINS     : RDKit FilterCatalog PAINS（A+B+C）重算
  QED       : rdkit.Chem.QED.qed 默认权重
  化学型     : 2 个 SMARTS 标志（氨基烟腈母核 × 碱性胺）→ 4 类
  AM1-BCC 风险档: 芳香 N ≥ 4 且 N/重原子 ≥ 0.25
  骨架       : Bemis-Murcko（MurckoScaffoldSmiles，无手性/同位素）
  指纹       : ECFP4 = Morgan r=2, 2048 bit（主）；ECFP6 = r=3；MACCS 166（敏感性）
  接触判据   : 重原子间距 < 4.0 Å
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, FilterCatalog, Lipinski, QED
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

# ---------------------------------------------------------------- 常量
CONTACT_D = 4.0
BOX_CENTER = (-0.471, 8.935, 17.159)
BOX_EDGE = 23.042

POCKET16 = {84: "V84", 85: "L85", 168: "F168", 169: "E169", 174: "M174", 177: "M177",
            181: "N181", 246: "W246", 249: "L249", 250: "H250", 253: "N253",
            264: "H264", 267: "L267", 270: "M270", 271: "Y271", 274: "I274"}
ANCHORS = ["N253", "E169", "F168", "W246"]

SMARTS_CORE = "[NX3]c1ncccc1C#N"          # 氨基烟腈母核
SMARTS_AMINE = "[NX3;H2,H1;!$(N-C=O);!$(N-S(=O)=O)]"   # 碱性胺（非酰胺/磺酰胺）

_LFC = rdMolStandardize.LargestFragmentChooser()


def _catalog(kind: str):
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(getattr(FilterCatalog.FilterCatalogParams.FilterCatalogs, kind))
    return FilterCatalog.FilterCatalog(params)


_PAINS = _catalog("PAINS")
_BRENK = _catalog("BRENK")


# ---------------------------------------------------------------- SDF 解析
def molblock_data_lines(block: str) -> list:
    """取 SDF 块的「数据部分」并保证以 `M  END` 结束。

    V2000 molblock 的数据部分 = 标题/程序/注释 3 行 + counts + 原子/键块 + 属性块（M  ...），
    最后必须由 `M  END` 收尾。旧实现按 counts 计算切片，会丢掉 `M  END`，
    写出的 SDF 虽有 $$$$ 分隔但块本身非法（RDKit 报 "Cannot convert '> <' to int"）。
    """
    body = block.split("\n")
    ci = next((i for i, l in enumerate(body) if l.rstrip().endswith("V2000")), None)
    if ci is not None:
        na = int(body[ci][0:3] or 0)
        nb = int(body[ci][3:6] or 0)
        lo = ci + 1 + na + nb
    else:
        lo = 0
    me = next((i for i in range(lo, len(body)) if body[i].strip() == "M  END"), None)
    if me is not None:
        return body[: me + 1]
    return body[:lo] + ["M  END"]


def iter_sdf_blocks(path):
    """按 $$$$ 切分，返回 (block_text, props, title)。不依赖 RDKit，稳健读大文件。"""
    buf: list[str] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("$$$$"):
                if buf:
                    yield _block_payload(buf)
                buf = []
            else:
                buf.append(line.rstrip("\n"))
    if buf and any(x.strip() for x in buf):
        yield _block_payload(buf)


def _block_payload(lines: list[str]):
    props, key = {}, None
    for ln in lines:
        if ln.startswith("> <"):
            key = ln[3:].split(">")[0].strip("< ")
        elif ln.strip() == "":
            key = None
        elif key is not None:
            props[key] = ln.strip()
            key = None
    return "\n".join(lines), props, lines[0].strip()


def block_coords(block: str):
    """返回 (heavy_atom_xyz[N,3], 元素列表, 是否含 H)。仅取 V2000 计数行。"""
    lines = block.split("\n")
    ci = next((i for i, l in enumerate(lines) if l.rstrip().endswith("V2000")), None)
    if ci is None:
        return None
    na = int(lines[ci][0:3])
    xyz, elems = [], []
    for l in lines[ci + 1: ci + 1 + na]:
        p = l.split()
        if len(p) < 4:
            continue
        xyz.append([float(p[0]), float(p[1]), float(p[2])])
        elems.append(p[3])
    xyz = np.asarray(xyz, dtype=float)
    heavy = np.asarray([e != "H" for e in elems])
    return xyz[heavy], elems, bool((~heavy).any())


def mol_from_block(block, sanitize: bool = False, remove_hs: bool = True):
    try:
        m = Chem.MolFromMolBlock(block, sanitize=sanitize, removeHs=False)
        if m is None:
            return None
        return Chem.RemoveHs(m, sanitize=False) if remove_hs else m
    except Exception:
        return None


# ---------------------------------------------------------------- 身份核验
def canon_smiles(mol, iso: bool = True):
    try:
        return Chem.MolToSmiles(mol, isomericSmiles=iso)
    except Exception:
        return None


def inchi_key(mol):
    try:
        return Chem.MolToInchiKey(mol)
    except Exception:
        return None


def template_main(smiles: str):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, 0
    frags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True)
    return max(frags, key=lambda x: x.GetNumHeavyAtoms()), len(frags)


def standardize_main_fragment(mol):
    """结构规范化：去盐 / 取最大片段（保留原电荷态，**不中和**）。

    口径（frozen_rule.step2_druglikeness.standardization）：
      描述符与门所描述的物种 = 对接/模拟所处理的物种。成盐输入（如盐 + 反离子）会使
      MW/cLogP/TPSA/QED 与 Wo5/Veber/BRENK 判定被反离子污染；本函数把分子规范化为
      「最大片段」，并把反离子信息作为独立列保留（可逆、可审计）。

    返回 (mol_main, n_fragments, counter_ion_smiles)：
      - n_fragments == 1 时**原样返回**输入 mol（保证单片段分子的描述符逐位不变）；
      - counter_ion_smiles = 其余片段按重原子数降序的规范 SMILES 拼接（无者为 ""）。
    """
    if mol is None:
        return None, 0, ""
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    n = len(frags)
    if n <= 1:
        return mol, n, ""
    main = _LFC.choose(mol)
    others = [f for f in frags if canon_smiles(f) != canon_smiles(main)]
    others.sort(key=lambda x: (-x.GetNumHeavyAtoms(), canon_smiles(x) or ""))
    return main, n, ";".join(canon_smiles(f) or "" for f in others)


def classify_pose(pose_mol, smiles: str) -> dict:
    """姿态 vs 上游 SMILES 的结构身份分级。

    identical / stereo_only / protonation_diff / graph_mismatch_* / parse_fail
    判定顺序：连接性+键级（去立体）→ 立体（3D 指派后）→ 修正键级后重比。
    """
    out = {"qc_class": "parse_fail", "identity_class": "parse_fail", "stereo_source": "",
           "protonation_note": "", "n_fragments_template": 0, "heavy_atom_delta": None}
    if pose_mol is None:
        return out
    tmpl, nfrag = template_main(smiles)
    out["n_fragments_template"] = nfrag
    if tmpl is None:
        out["qc_class"] = out["identity_class"] = "smiles_parse_fail"
        return out
    out["heavy_atom_delta"] = int(pose_mol.GetNumHeavyAtoms() - tmpl.GetNumHeavyAtoms())
    pose_c = canon_smiles(pose_mol, iso=False)
    tmpl_c = canon_smiles(tmpl, iso=False)
    if pose_c == tmpl_c:
        p2 = Chem.Mol(pose_mol)
        try:
            Chem.AssignStereochemistryFrom3D(p2)
        except Exception:
            pass
        if canon_smiles(p2) == canon_smiles(tmpl):
            out.update(qc_class="identical", identity_class="identical",
                       stereo_source="pose3D")
        else:
            out.update(qc_class="stereo_only", identity_class="stereo_only",
                       stereo_source="pose3D(未在输入中指定)")
        return out
    # 键级/质子化差异 → 尝试模板修正
    fixed = None
    try:
        fixed = AllChem.AssignBondOrdersFromTemplate(tmpl, pose_mol)
    except Exception:
        fixed = None
    if fixed is not None and canon_smiles(fixed, iso=False) == tmpl_c:
        try:
            Chem.AssignStereochemistryFrom3D(fixed)
        except Exception:
            pass
        same_stereo = canon_smiles(fixed) == canon_smiles(tmpl)
        out.update(qc_class="identical_fixed" if same_stereo else "stereo_only_fixed",
                   identity_class="identical" if same_stereo else "stereo_only",
                   stereo_source="pose3D(修正键级后)")
        return out
    # 仍不一致 → 分类：质子化/互变 vs 图不一致
    dH = out["heavy_atom_delta"]
    charge_pose = Chem.GetFormalCharge(pose_mol)
    charge_tmpl = Chem.GetFormalCharge(tmpl)
    taut = _tautomer_flag(tmpl) or _tautomer_flag(pose_mol)
    if dH == 0 and charge_pose != charge_tmpl:
        out.update(qc_class="protonation_diff", identity_class="protonation_diff",
                   protonation_note=(
            f"形式电荷差 {charge_pose - charge_tmpl:+d}；母核含脒/亚胺/胍等可质子化基团={bool(taut)}"),
            stereo_source="pose3D(未在输入中指定)")
    elif _same_skeleton(pose_mol, tmpl):
        out.update(qc_class="protonation_diff", identity_class="protonation_diff",
                   protonation_note=(
            f"InChIKey 骨架层一致（前 14 位），氢/电荷状态不同（ΔH={dH}）"),
            stereo_source="pose3D(未在输入中指定)")
    else:
        out.update(qc_class="graph_mismatch", identity_class="graph_mismatch",
                   protonation_note="")
    return out


def _tautomer_flag(mol) -> bool:
    patt = Chem.MolFromSmarts("[NX3][CX3]=[NX2]")   # 脒/亚胺
    guan = Chem.MolFromSmarts("[NX3][CX3](=[NX2])[NX3]")
    return bool(mol.HasSubstructMatch(patt) or mol.HasSubstructMatch(guan))


def _same_skeleton(a, b) -> bool:
    ka, kb = inchi_key(a), inchi_key(b)
    if not ka or not kb:
        return False
    return ka.split("-")[0] == kb.split("-")[0]


# ---------------------------------------------------------------- 描述符
def base_descriptors(mol) -> dict:
    return {
        "MW": Descriptors.MolWt(mol),
        "cLogP": Crippen.MolLogP(mol),
        "HBD": int(Lipinski.NumHDonors(mol)),
        "HBA": int(Lipinski.NumHAcceptors(mol)),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "RotB": int(Lipinski.NumRotatableBonds(mol)),
        "n_heavy": int(mol.GetNumHeavyAtoms()),
        "n_N": int(sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "N")),
        "n_aromN": int(sum(1 for a in mol.GetAtoms()
                           if a.GetSymbol() == "N" and a.GetIsAromatic())),
        "n_rings": int(rdMolDescriptors.CalcNumRings(mol)),
        "n_aromRings": int(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "fsp3": float(rdMolDescriptors.CalcFractionCSP3(mol)),
        "formal_charge": int(Chem.GetFormalCharge(mol)),
    }


def druglikeness(mol) -> dict:
    d = base_descriptors(mol)
    d["ro5_viol"] = int((d["MW"] > 500) + (d["cLogP"] > 5) + (d["HBD"] > 5) + (d["HBA"] > 10))
    d["ro5_detail"] = ";".join(
        f"{k}={'viol' if v > t else 'ok'}" for k, v, t in
        (("MW", d["MW"], 500), ("cLogP", d["cLogP"], 5),
         ("HBD", d["HBD"], 5), ("HBA", d["HBA"], 10)))
    d["veber_detail"] = ";".join(
        f"{k}={'viol' if v > t else 'ok'}" for k, v, t in
        (("RotB", d["RotB"], 10), ("TPSA", d["TPSA"], 140)))
    d["veber_viol"] = int((d["RotB"] > 10) + (d["TPSA"] > 140))
    d["qed"] = float(QED.qed(mol))
    d["pains"] = int(_PAINS.HasMatch(mol))
    d["pains_names"] = ";".join(sorted({e.GetDescription() for e in _PAINS.GetMatches(mol)}))
    d["brenk"] = int(_BRENK.HasMatch(mol))
    d["brenk_names"] = ";".join(sorted({e.GetDescription() for e in _BRENK.GetMatches(mol)}))
    d["risk_am1bcc"] = int(d["n_aromN"] >= 4 and d["n_N"] / max(d["n_heavy"], 1) >= 0.25)
    return d


def chemotype(mol) -> str:
    core = mol.HasSubstructMatch(Chem.MolFromSmarts(SMARTS_CORE))
    amine = mol.HasSubstructMatch(Chem.MolFromSmarts(SMARTS_AMINE))
    return (f"{'氨基烟腈母核' if core else '无氨基烟腈母核'} + "
            f"{'碱性胺' if amine else '无碱性胺'}")


def scaffold(mol) -> str:
    try:
        s = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        return s
    except Exception:
        return ""


def scaffold_generic(smiles: str) -> str:
    """骨架的“通用化”形式（原子类型/键级抹平），仅用于报告双口径计数。"""
    sc = scaffold(Chem.MolFromSmiles(smiles)) if isinstance(smiles, str) else ""
    if not sc:
        return ""
    ms = Chem.MolFromSmiles(sc)
    if ms is None:
        return ""
    try:
        return Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(ms))
    except Exception:
        return ""


# ---------------------------------------------------------------- 指纹
def fp_ecfp4(mol, n_bits: int = 2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)


def fp_ecfp6(mol, n_bits: int = 2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=n_bits)


def fp_maccs(mol):
    return AllChem.GetMACCSKeysFingerprint(mol)


def fp_to_array(fp):
    a = np.zeros((fp.GetNumBits(),), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fp, a)
    return a


def tanimoto_matrix(fps) -> np.ndarray:
    n = len(fps)
    m = np.zeros((n, n), dtype=np.float32)
    for i in range(n):
        m[i] = DataStructs.BulkTanimotoSimilarity(fps[i], fps)
    return m


def pairwise_tc(fps) -> np.ndarray:
    m = tanimoto_matrix(fps)
    np.fill_diagonal(m, 0.0)
    return m


# ---------------------------------------------------------------- 受体 / IFP
def load_receptor(path, residue_map: dict[int, str] | None = None) -> dict[int, np.ndarray]:
    res = {k: [] for k in (residue_map or POCKET16)}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if ln.startswith(("ATOM", "HETATM")):
                try:
                    rn = int(ln[22:26])
                except ValueError:
                    continue
                if rn in res:
                    res[rn].append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return {k: np.asarray(v, dtype=float) for k, v in res.items() if v}


def block_heavy(block: str):
    """返回 (heavy_xyz[N,3], heavy_elems[N], has_H)：元素与坐标**对齐**的重原子视图。

    与 `block_coords` 的区别：后者返回完整元素列表（含 H），元素与重原子坐标不对齐，
    适合只看坐标的场景；方向性/成键判据需要逐原子元素，故单列本函数。
    """
    lines = block.split("\n")
    ci = next((i for i, l in enumerate(lines) if l.rstrip().endswith("V2000")), None)
    if ci is None:
        return None
    na = int(lines[ci][0:3])
    xyz, elems = [], []
    for l in lines[ci + 1: ci + 1 + na]:
        p = l.split()
        if len(p) < 4:
            continue
        if p[3] == "H":
            continue
        xyz.append([float(p[0]), float(p[1]), float(p[2])])
        elems.append(p[3])
    return np.asarray(xyz, dtype=float), elems, len(elems) < na


def anchor_directionality(coords: np.ndarray, elems: list, res_atoms: dict,
                          polar=("N253", "E169"), bond_d: float = 1.9) -> dict:
    """极性锚点接触的**方向性代理**（无氢几何时的重原子近似）。

    判据：对每个极性锚点 A，取配体最近的 N/O 重原子 D；在配体内部取距 D 最近且 < bond_d 的
    重原子 X 作为「共价键方向」；角度 θ = ∠(X–D–A)，线性排列（理想的 D–H···A）为 180°。
      - θ 接近 180° → 接触方向与 D 的成键方向共线，符合氢键/盐桥几何；
      - θ 明显偏离（< 120°）→ 更像非特异接触（疏水挤压/口袋壁贴合）。
    为什么需要它：`anchor_hits_4A` 只判「4 Å 内有没有」，不判方向；两个都命中 N253 的姿态，
    一个可能是氢键、另一个可能是从侧面蹭过。本函数为姿态质量提供第二维证据（**不参与落门**）。

    返回：`angle_deg_<锚点>`（无极性接触为 nan）、`polar_contact_n`、`angle_deg_min`。
    """
    out = {f"angle_deg_{a}": float("nan") for a in polar}
    vals = []
    if coords is None or coords.size == 0:
        out["polar_contact_n"] = 0
        out["angle_deg_min"] = float("nan")
        return out
    for a in polar:
        ra = res_atoms.get(_name_to_resnum(a))
        if ra is None:
            continue
        el = np.asarray([e != "C" for e in elems], dtype=bool)
        if not el.any():
            continue
        d_mat = np.linalg.norm(coords[:, None, :] - ra[None, :, :], axis=2)
        d_min = d_mat.min(axis=1)
        cand = np.where(el)[0]
        if cand.size == 0:
            continue
        d_idx = int(cand[np.argmin(d_min[cand])])
        if d_min[d_idx] >= CONTACT_D:
            continue
        # 配体内部：D 的最近重原子邻居（< bond_d）
        dd = np.linalg.norm(coords - coords[d_idx], axis=1)
        dd[d_idx] = np.inf
        n_idx = int(np.argmin(dd))
        if not np.isfinite(dd[n_idx]) or dd[n_idx] >= bond_d:
            continue
        a_idx = int(np.argmin(d_mat[d_idx]))
        v1 = coords[n_idx] - coords[d_idx]        # D → X（与 D→A 同向时 θ→180°）
        v2 = ra[a_idx] - coords[d_idx]            # D → A
        c = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12))
        ang = float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))
        out[f"angle_deg_{a}"] = round(ang, 1)
        vals.append(ang)
    out["polar_contact_n"] = len(vals)
    out["angle_deg_min"] = round(min(vals), 1) if vals else float("nan")
    return out


def _name_to_resnum(name: str) -> int:
    for k, v in POCKET16.items():
        if v == name:
            return k
    return -1


def pose_quality_tier(anchor_hits: int, jaccard, thresh: float = 0.60) -> str:
    """姿态质量档（**报告/优先级用，不参与落门**）。

    A：anchor=4 且与共晶接触集的 Jaccard ≥ thresh（既贴 4 个锚点，又与晶体样结合模式相似）
    B：anchor=4 但 Jaccard < thresh
    C：anchor=3
    D：anchor ≤ 2
    动机：`anchor≥3` 一档内的异质性很大（实测 P(active|anchor=3) 低于全集基线），
    必须给出比「命中数」更细的排序量，MD 优先级与人工复核才有依据。
    """
    if anchor_hits >= 4:
        try:
            return "A" if float(jaccard) >= thresh else "B"
        except (TypeError, ValueError):
            return "B"
    return "C" if anchor_hits == 3 else "D"


def read_pdb_heavy_coords(path, resname: str | None = None) -> np.ndarray:
    """读取 PDB 中某配体（HETATM，按残基名过滤）或全部重原子的坐标 [N,3]。

    用于共晶配体（正对照）与受体口袋的几何参照计算：只取非氢原子行。
    """
    xyz = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if not ln.startswith(("ATOM", "HETATM")):
                continue
            if resname is not None and ln[17:20].strip() != resname:
                continue
            el = ln[76:78].strip() or ln[12:16].strip()[:1]
            if el.upper() == "H":
                continue
            try:
                xyz.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            except ValueError:
                continue
    return np.asarray(xyz, dtype=float)


def read_pdb_heavy(path, resname: str | None = None):
    """返回 (heavy_xyz[N,3], heavy_elems[N])：PDB 中某配体的重原子坐标与元素（对齐）。"""
    xyz, elems = [], []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if not ln.startswith(("ATOM", "HETATM")):
                continue
            if resname is not None and ln[17:20].strip() != resname:
                continue
            el = ln[76:78].strip() or ln[12:16].strip()[:1]
            if el.upper() == "H":
                continue
            try:
                xyz.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
                elems.append(el)
            except ValueError:
                continue
    return np.asarray(xyz, dtype=float), elems


def contact_residues(coords: np.ndarray, res_atoms: dict[int, np.ndarray],
                     contact_d: float = CONTACT_D) -> list:
    """返回在 contact_d 内与配体有接触的残基名列表（按 POCKET16 顺序）。"""
    out = []
    for rn, name in POCKET16.items():
        ra = res_atoms.get(rn)
        if ra is None or coords.size == 0:
            continue
        if float(np.linalg.norm(coords[:, None, :] - ra[None, :, :], axis=2).min()) < contact_d:
            out.append(name)
    return out


def inchi_key(mol) -> str | None:
    """标准 InChIKey（用于「去盐后同一物种」的冗余登记；失败返回 None）。"""
    if mol is None:
        return None
    try:
        return Chem.MolToInchiKey(mol)
    except Exception:                                       # noqa: BLE001
        return None


def mol_wt(mol) -> float:
    """相对分子质量（与 druglikeness 口径一致，不圆整）。"""
    return float(Descriptors.MolWt(mol))


def ifp_row(coords: np.ndarray, res_atoms: dict[int, np.ndarray],
            contact_d: float = CONTACT_D) -> dict:
    row, contacts = {}, []
    for rn, name in POCKET16.items():
        ra = res_atoms.get(rn)
        if ra is None or coords.size == 0:
            d = float("nan")
        else:
            d = float(np.linalg.norm(coords[:, None, :] - ra[None, :, :], axis=2).min())
        row[f"d_{name}"] = round(d, 2)
        if d < contact_d:
            contacts.append(name)
    row["contacts_4A"] = ";".join(contacts)
    row["n_contacts_4A"] = len(contacts)
    row["anchor_hits_4A"] = int(sum(1 for a in ANCHORS if row[f"d_{a}"] < contact_d))
    return row


def box_geometry(coords: np.ndarray, center=BOX_CENTER, edge: float = BOX_EDGE) -> dict:
    c = np.asarray(center, dtype=float)
    half = edge / 2.0
    excess = float((np.abs(coords - c).max(axis=1) - half).max()) if coords.size else float("nan")
    return {"box_max_excess": round(excess, 3), "outside_box": int(excess > 0)}


def pose_gate(anchor_hits: int, outside_box: int) -> str:
    """reject: 出盒 或 锚点 ≤1 ；warn: 锚点 = 2 ；pass: 锚点 ≥3"""
    if outside_box or anchor_hits <= 1:
        return "reject"
    return "warn" if anchor_hits == 2 else "pass"


# ---------------------------------------------------------------- 表工具
def add_cns_inputs(df: pd.DataFrame, upstream: pd.DataFrame) -> pd.DataFrame:
    """把 pKa / cLogD7.4 并入候选表（CNS MPO 的两个输入，来源为上游表，非本包重算）。"""
    cols = ["cid", "basic_pka", "clogd74", "pka_status", "pka_basic_center_count",
            "clogd_status", "clogd_confidence", "pka_provider", "pka_provider_version",
            "clogd_provider", "clogd_provider_version"]
    cols = [c for c in cols if c in upstream.columns]
    return df.merge(upstream[cols], on="cid", how="left")


def read_smiles_table(root_input_csv: str) -> pd.DataFrame:
    up = pd.read_csv(root_input_csv, low_memory=False)
    up["ik_upstream"] = [inchi_key(Chem.MolFromSmiles(s)) if isinstance(s, str) else None
                         for s in up["smiles"]]
    return up

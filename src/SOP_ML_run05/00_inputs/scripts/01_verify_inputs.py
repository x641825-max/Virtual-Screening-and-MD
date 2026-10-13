#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步骤 01：输入固化与规则冻结。

做三件事：
  ① 5 个原始输入的 md5/sha256 与 frozen_rule.inputs 对账（不一致 → 退出码 4）；
  ② 登记关键派生输入哈希（受体 PDB、5,098 带标签库、ChEMBL 独立集）；
  ③ 规则冻结：对 frozen_rule.json（剔除 rules_sha256 字段）取规范 JSON 哈希，
     首次运行写入该字段，此后每次校验一致；同时登记全部脚本 sha256（FREEZE.json）。
输出：00_inputs/manifests/INPUT_MANIFEST.md、input_verify.json、FREEZE.json
"""
from __future__ import annotations

import sys
from pathlib import Path

_p = Path(__file__).resolve()
while not (_p / "REPRO_INFO.json").is_file():
    _p = _p.parent
sys.path.insert(0, str(_p))
ROOT = _p

from vendor.io_utils import (StepLog, canonical_hash, md5_file, read_json, sha256_file,
                             write_json, write_text)

RAW = ROOT / "00_inputs/raw"
MAN = ROOT / "00_inputs/manifests"
KEY_DERIVED = {
    "00_inputs/extracted/run01/01_redocking/receptor/4EIY_protein.pdb": "aef775ee83735af8a878d42ca1e5fa85",
}

STEP_DIRS = ["00_inputs", "01_pose_qc", "02_druglikeness", "03_gates", "04_selection",
             "05_robustness", "06_deliverables"]


def _all_code_files(root: Path) -> list[Path]:
    """包内全部可执行代码：步骤脚本 + vendor + tests + run_all.py。"""
    files = []
    for d in STEP_DIRS:
        files += [p for p in (root / d).glob("scripts/*.py")]
    for d in ("vendor", "tests"):
        files += [p for p in (root / d).glob("*.py")]
    if (root / "run_all.py").is_file():
        files.append(root / "run_all.py")
    return files


def main() -> int:
    rule = read_json(ROOT / "config/frozen_rule.json")
    with StepLog(ROOT, "01_verify_inputs") as log:
        rows, ok = [], True
        for rel, spec in rule["inputs"].items():
            path = RAW / ("archives/" + rel if rel.endswith(".zip") else "inputs/" + rel)
            if not path.is_file():
                # 换数据版骨架未含大体积输入（见 RESTORE_LIST.md）：登记缺失并显式失败，
                # 不静默跳过、也不以 traceback 终止
                ok = False
                rows.append({"file": str(path.relative_to(ROOT)), "role": spec["role"],
                             "md5_expected": spec["md5"], "md5_got": "MISSING",
                             "sha256_got": None, "match": False, "size": 0})
                log.info(f"[缺失] {rel}  未找到 {path.relative_to(ROOT)}"
                         f"（恢复方式见 RESTORE_LIST.md）")
                continue
            got_md5, got_sha = md5_file(path), sha256_file(path)
            good = got_md5 == spec["md5"] and got_sha == spec["sha256"]
            ok &= good
            rows.append({"file": str(path.relative_to(ROOT)), "role": spec["role"],
                         "md5_expected": spec["md5"], "md5_got": got_md5,
                         "sha256_got": got_sha, "match": bool(good),
                         "size": path.stat().st_size})
            log.info(f"{'OK  ' if good else 'FAIL'} {rel}  md5={got_md5[:12]}…  "
                     f"{path.stat().st_size / 1e6:.1f} MB")
        for rel, exp in KEY_DERIVED.items():
            p = ROOT / rel
            if not p.is_file():
                log.info(f"[跳过] 派生输入未解包：{rel}（先运行 00_unpack）")
                ok = False
                continue
            got = md5_file(p)
            good = got == exp
            ok &= good
            rows.append({"file": rel, "role": "派生输入（解包产物）", "md5_expected": exp,
                         "md5_got": got, "sha256_got": sha256_file(p), "match": bool(good),
                         "size": p.stat().st_size})
            log.info(f"{'OK  ' if good else 'FAIL'} {rel}  md5={got[:12]}…")

        # 关键成员（标签集 / 独立已知配体集）
        tagged = ROOT / "00_inputs/extracted/run01/02_enrichment_5098/library/quick_5100_library.csv"
        chembl = ROOT / "00_inputs/extracted/screening_kit/07_independent_validation/independent_a2a_chembl.csv"
        n_tag = sum(1 for _ in open(tagged, encoding="utf-8")) - 1 if tagged.is_file() else 0
        n_ch = sum(1 for _ in open(chembl, encoding="utf-8")) - 1 if chembl.is_file() else 0
        log.info(f"带标签富集库 = {n_tag} 条；ChEMBL 独立已知配体 = {n_ch} 条")

        # ---- 外部接口输入（hERG 预测表，登记后校验哈希）
        herg = rule.get("step3_gates", {}).get("herg", {})
        hp = ROOT / herg.get("external_table", "03_gates/herg/herg_predictions.csv")
        herg_entry = {"path": herg.get("external_table"), "present": bool(hp.is_file()),
                      "registered_sha256": herg.get("external_sha256"),
                      "registered_md5": herg.get("external_md5")}
        if hp.is_file():
            herg_entry["sha256"] = sha256_file(hp)
            herg_entry["md5"] = md5_file(hp)
            if herg.get("external_sha256"):
                match = herg_entry["sha256"] == herg["external_sha256"]
                herg_entry["match"] = bool(match)
                if not match:
                    log.info(f"[致命] hERG 预测表与规则登记哈希不一致："
                             f"文件 {herg_entry['sha256'][:16]}… ≠ 登记 "
                             f"{herg['external_sha256'][:16]}…")
                    ok = False
                else:
                    log.info(f"hERG 预测表哈希一致：{herg_entry['sha256'][:16]}…"
                             f"（来源 {herg.get('external_source', '未标注')}）")
            else:
                herg_entry["match"] = None
                log.info(f"hERG 预测表存在但未登记哈希：{herg_entry['sha256'][:16]}…")
        else:
            log.info("hERG 预测表不存在 → 该门保持 pending（不含于输入登记）")

        # ---- 规则冻结
        rules_body = {k: v for k, v in rule.items() if k != "rules_sha256"}
        h = canonical_hash(rules_body)
        if rule.get("rules_sha256", "PENDING") == "PENDING":
            txt = (ROOT / "config/frozen_rule.json").read_text(encoding="utf-8")
            txt = txt.replace('"rules_sha256": "PENDING"', f'"rules_sha256": "{h}"')
            (ROOT / "config/frozen_rule.json").write_text(txt, encoding="utf-8")
            log.info(f"首次冻结：rules_sha256 = {h}")
            frozen_now = True
        else:
            frozen_now = False
            if rule["rules_sha256"] != h:
                log.info(f"[致命] frozen_rule 被修改：文件 {rule['rules_sha256'][:16]}… ≠ "
                         f"重算 {h[:16]}…（改规则必须升 rule_version 并重新冻结）")
                ok = False
        script_hashes = {}
        for s in sorted(_all_code_files(ROOT)):
            script_hashes[str(s.relative_to(ROOT))] = sha256_file(s)
        freeze_path = MAN / "FREEZE.json"
        prev = read_json(freeze_path) if freeze_path.is_file() else None
        write_json(freeze_path, {"rules_sha256": h, "rule_version": rule["rule_version"],
                                 "frozen_now": frozen_now,
                                 "scripts": script_hashes})
        if prev and prev.get("scripts") != script_hashes:
            changed = sorted(k for k in set(prev.get("scripts", {})) | set(script_hashes)
                             if prev.get("scripts", {}).get(k) != script_hashes.get(k))
            log.info(f"[提示] 与上次运行相比脚本哈希变化 {len(changed)} 个：{changed[:6]}")

        write_json(MAN / "input_verify.json",
                   {"all_ok": bool(ok), "rules_sha256": h,
                    "n_inputs": len(rule["inputs"]), "entries": rows,
                    "tagged_library_rows": n_tag, "chembl_known_rows": n_ch,
                    "external_interface_inputs": {"hERG": herg_entry}})

        pkg_name = read_json(ROOT / "REPRO_INFO.json").get("package", "pack")
        md = [f"# 输入固化清单（{pkg_name}）", "",
              f"- 规则集：`{rule['rule_set_id']}` v{rule['rule_version']}",
              f"- rules_sha256：`{h}`", f"- 全部输入哈希一致：**{'是' if ok else '否'}**", "",
              "| 文件 | 角色 | md5 | 校验 | 大小(MB) |", "|---|---|---|---|---|"]
        for r in rows:
            md.append(f"| `{r['file']}` | {r['role']} | `{r['md5_got']}` | "
                      f"{'一致' if r['match'] else '**不一致**'} | {r['size'] / 1e6:.2f} |")
        md += ["", "## 关键派生输入", "",
               f"- 5,098 带标签富集库（479 活性 / 4,619 诱饵）："
               f"`00_inputs/extracted/run01/02_enrichment_5098/library/quick_5100_library.csv`"
               f"（{n_tag} 行）",
               f"- ChEMBL 独立已知 A2A 配体："
               f"`00_inputs/extracted/screening_kit/07_independent_validation/independent_a2a_chembl.csv`"
               f"（{n_ch} 行）", "",
               "> 本包内所有引用均位于包根之下；脚本禁止引用包外路径（`tests/check_isolation.py` 静态扫描）。"]
        write_text(MAN / "INPUT_MANIFEST.md", "\n".join(md))
        log.info(f"输入固化完成：all_ok={ok}；清单 → 00_inputs/manifests/INPUT_MANIFEST.md")
    return 0 if ok else 4


if __name__ == "__main__":
    raise SystemExit(main())

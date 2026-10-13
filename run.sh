#!/usr/bin/env bash
# 一键入口（赛事附件 5 §二(2)）：从随包数据复跑「MD 候选选择 → 联立合并 → 出结果文件」，
# 并与随包冻结的候选清单逐行逐字段核验。
#
#   ./run.sh            默认：复跑 MD 候选路径 + 自检（约 2–4 分钟，纯 CPU）
#   ./run.sh --full     完整链（需要先按 data/DATA_SOURCES.md 恢复上游原始输入；两链实测各 ≈300 s）
#   ./run.sh --check    只做环境自检，不跑链
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
export PYTHONHASHSEED=0                 # 与规则集一致的确定性口径
export SOPML_ROOT="$HERE/src/SOP_ML_run05"   # 另一条链的包根（脚本默认也会按同级解析）

mkdir -p logs
LOG="logs/run.log"
: > "$LOG"
say() { echo "$@" | tee -a "$LOG"; }

say "== [1/5] 环境自检 =="
python3 - <<'PY' 2>&1 | tee -a "$LOG"
import sys, importlib
need = {"numpy": "2.4.6", "pandas": "2.3.3", "scipy": "1.15.3", "sklearn": "1.9.1",
        "matplotlib": "3.11.1", "rdkit": "2025.09.5", "openpyxl": "3.1.5"}
print("python", sys.version.split()[0])
bad = []
for m, want in need.items():
    try:
        v = getattr(importlib.import_module(m), "__version__", "?")
    except Exception as e:
        v = f"MISSING({e})"
    ok = (v == want)
    if not ok:
        bad.append(f"{m}: {v} != {want}")
    print(f"  {m:12s} {v}  {'OK' if ok else 'DIFF(期望 ' + want + ')'}")
if bad:
    print("环境与冻结快照有差异（不阻断运行，仅提示）：", "; ".join(bad))
PY

if [ "${1:-}" = "--check" ]; then say "仅环境自检，结束。"; exit 0; fi

if [ "${1:-}" = "--full" ]; then
  say "== 完整链模式：需要 00_inputs/raw/ 与 00_inputs/extracted/（见 data/DATA_SOURCES.md）=="
  ( cd src/SOP_5        && python3 run_all.py ) 2>&1 | tee -a "$LOG"
  ( cd src/SOP_ML_run05 && python3 run_all.py ) 2>&1 | tee -a "$LOG"
  ( cd src/SOP_5        && python3 10_final_md/scripts/17_final_md_merge.py \
                        && python3 10_final_md/scripts/18_final_md_figures.py ) 2>&1 | tee -a "$LOG"
else
  say "== [2/5] 复跑 SOP_5 目录化合物链 · MD 候选选择层 =="
  ( cd src/SOP_5 && python3 07_md_batch/scripts/15_md_batch_select.py ) 2>&1 | tee -a "$LOG"

  say "== [3/5] 复跑 SOP_ML_run05 生成分子链 · MD 候选选择层 =="
  ( cd src/SOP_ML_run05 && python3 07_md_batch/scripts/15_md_batch_select.py ) 2>&1 | tee -a "$LOG"

  say "== [4/5] 复跑联立合并层（两来源并池重解 20 席）=="
  ( cd src/SOP_5 && python3 10_final_md/scripts/17_final_md_merge.py \
                 && python3 10_final_md/scripts/18_final_md_figures.py ) 2>&1 | tee -a "$LOG"
fi

say "== [5/5] 生成候选清单并核验 =="
python3 predict.py --source repro --out results/reproduced 2>&1 | tee -a "$LOG"
python3 predict.py --verify results/reproduced/results.csv results/results.csv 2>&1 | tee -a "$LOG"
say "完成。复跑清单：results/reproduced/results.csv ；冻结清单：results/results.csv"

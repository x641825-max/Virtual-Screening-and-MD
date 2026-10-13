# data/example · 小规模示例数据

按赛事《代码提交要求》附件 5 §一（演示案例）与 §五（数据受限时提供可复现流程）提供。
三份文件都是**真实数据的切片或副本**，可直接用于冒烟测试与 Notebook 演示。

| 文件 | 行数 | 内容 | 来源 |
|---|---:|---|---|
| `list_main_QED0.50_anchor4_head50.csv` | 50 | 目录链主口径名单（48 候选 + 2 阴性对照），选批层的直接输入 | `src/SOP_5/04_selection/tables/list_main_QED0.50_anchor4.csv` 前 50 行 |
| `candidate_universe_head200.csv` | 200 | 资格池逐分子表（含描述符/骨架/门结果），分位切层与约束校验的输入 | `src/SOP_5/03_gates/tables/candidate_universe.csv` 前 200 行 |
| `final_md_20.csv` | 22 | 联立合并层最终名册（20 候选 + 2 共用对照），含全字段 | `results/merged_final_md/final_md_20.csv` 完整副本 |

## 用途

```bash
# 冒烟测试：确认环境与本包的读取路径可用
python3 - <<'PY'
import pandas as pd
d = pd.read_csv("data/example/list_main_QED0.50_anchor4_head50.csv")
print(d.shape, d.columns[:6].tolist())
print(d.candidate_type.value_counts().to_dict())
PY
```

更完整的演示见 `notebooks/md_candidates_walkthrough.ipynb`：它用 `final_md_20.csv` 复算
「候选 / 对照」构成、来源构成、CNN_VS 分位覆盖与 hERG 双模型分歧，并重跑 `predict.py`
生成统一候选清单。

## 说明

- 示例切片**保留原始列名与数值精度**，未做任何再处理；行数截断只影响规模，不影响列口径。
- 三条链的完整输入与中间产物在 `src/` 下随包发布（体积可控），因此**不需要**示例数据即可复跑
  选择层与合并层；示例数据主要用于环境冒烟与阅读演示。
- 上游大体积原始数据未随包分发，获取方式与哈希见 `data/DATA_SOURCES.md`。

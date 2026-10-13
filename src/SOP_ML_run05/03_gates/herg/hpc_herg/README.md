# hpc_herg/：外部 hERG 预测的原始产出（惰性证据）

本目录存放生成 **落门表** 的那次外部队列运行的原始文件。**链上没有任何脚本读取这里的文件**——
它只作为可审计证据保留；据此重新落门时，把需要的列整理成
`03_gates/herg/herg_predictions.csv`（列 `cid,hERG`）并更新
`config/frozen_rule.json → step3_gates.herg` 的 `external_md5` / `external_sha256` 即可。

| 文件 | 内容 | 说明 |
|---|---|---|
| `molecules.csv` | 3,236 行，列 `cid,smiles` | 四臂并集清单（md5 `3c117b47e6a4c645e229ef1ffe1eb757`） |
| `predictions.csv` | 3,236 行 × 106 列 | ADMET-AI 原始输出，含 `hERG`、`hERG_drugbank_approved_percentile` |
| `output.json` | `{"status":"ok","tool":"admet-ai","n":3236}` | 运行摘要 |

- 产出者：HPC 目录应用 `admet-ai`，job_id `3d9634dd50304d3db84b1f675fb209fc`（2026-09-23）。
- 落门表：`03_gates/herg/herg_predictions.csv`（md5 `41d15798e766be4389c43388877bacc6`，
  sha256 `e52533d8aca3e1f1fd2ab25b4afec4b1b49ae8f6645607293ff78d2401dfa498`）。
- 口径与残余不确定性：见 `03_gates/herg/HERG_INTERFACE.md` 与报告 §4/§9。
- 注意：ADMET-AI 的 `hERG` 与 `hERG_drugbank_approved_percentile` 两列来自**同一个模型**
  （后者只是把预测值换算成 DrugBank 已上市药物的百分位），**不构成独立的第二模型**；
  独立的第二模型交叉验证见 `03_gates/herg/cross_model/`。

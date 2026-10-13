# SOP_5 交付索引（06_deliverables）

- 规则集：`SOP-05-A2A-antagonist-shortlist` v1.7.0（rules_sha256 `a983c4985ebf1fde…`）
- 验收：**硬门 6/9 通过（pending 2）**；次要指标 3 项如实报告（详见 `acceptance_report.json`）

| 文件 | 说明 |
|---|---|
| `INDEX.md` | 本索引 |
| `SOP-05_HERG_THRESHOLD_SOP.md` | hERG 阈值选择 SOP（包内版本，含变更控制） |
| `SOP-05_MD_batch20.csv` | 见对应步骤文档 |
| `SOP-05_MD_batch20.md` | 见对应步骤文档 |
| `SOP-05_MD_batch_excluded.csv` | 见对应步骤文档 |
| `SOP-05_MD_shortlist_main.csv` | 主名单 50（48 候选 + 2 阴性对照），含逐分子门与注释 |
| `SOP-05_MD_shortlist_scaffoldcap3.csv` | 敏感性名单：骨架封顶 ≤3 |
| `SOP-05_MD_shortlist_scaffoldcap5.csv` | 敏感性名单：骨架封顶 ≤5 |
| `SOP-05_MD_shortlist_sens_QED0.40_anchor3.csv` | 敏感性名单：QED≥0.40 ∧ anchor≥3 |
| `SOP-05_MD_shortlist_sens_QED0.40_anchor4.csv` | 敏感性名单：QED≥0.40 ∧ anchor=4 |
| `SOP-05_MD_shortlist_sens_QED0.50_anchor3.csv` | 敏感性名单：QED≥0.50 ∧ anchor≥3 |
| `SOP-05_MD_shortlist_sens_QED0.60_anchor3.csv` | 敏感性名单：QED≥0.60 ∧ anchor≥3 |
| `SOP-05_MD_shortlist_sens_QED0.60_anchor4.csv` | 敏感性名单：QED≥0.60 ∧ anchor=4 |
| `SOP-05_admet_annotation.json` | 注释层证据（覆盖/重标定切点/触发率/分诊层/复核清单） |
| `SOP-05_admet_annotation_main.csv` | 交付名单 ADMET 注释全表（冻结表原精度，50 行） |
| `SOP-05_admet_descriptor_crosscheck.csv` | 描述符独立复算交叉核验（逐分子偏差与判定一致性） |
| `SOP-05_admet_descriptor_crosscheck.json` | 交叉核验汇总（容差与已知差异说明） |
| `SOP-05_admet_layer.md` | 注释层说明（步骤 10b） |
| `SOP-05_admet_md_precheck.csv` | MD 前置复核清单（交付候选中带 ≥1 项警告者） |
| `SOP-05_constraint_trigger_rates.csv` | 选序约束触发率（Tc/骨架/化学型/风险档） |
| `SOP-05_decision_cards_main.csv` | 主名单决策卡（含最近的已知活性、pKa/cLogD、PAINS/Brenk 标注） |
| `SOP-05_determinism_compare.json` | 确定性比对报告（G 臂；最近一次结果） |
| `SOP-05_downstream_MD_plan.md` | 下游 MD 接口与判据（含机时预算单列） |
| `SOP-05_gate_stack.csv` | 逐分子门叠加结果（4,990 条） |
| `SOP-05_herg_cross_model_consistency.csv` | 见对应步骤文档 |
| `SOP-05_herg_cross_model_list_flips.csv` | 见对应步骤文档 |
| `SOP-05_herg_cross_model_rank_consistency.csv` | 见对应步骤文档 |
| `SOP-05_herg_input_union.csv` | hERG 补算清单（四臂并集，外部接口） |
| `SOP-05_herg_retention_curve.csv` | hERG 阈值—保留曲线（0.01 步长） |
| `SOP-05_herg_status.json` | hERG 门状态（缺表时 pending） |
| `SOP-05_herg_threshold_analysis.json` | 阈值锚定证据（含输入 sha256 与复现核对） |
| `SOP-05_herg_threshold_anchors.csv` | hERG 阈值三锚定表（各头部 × min_keep 档取舍线） |
| `SOP-05_md_batch.json` | 见对应步骤文档 |
| `SOP-05_negative_controls.csv` | 阴性对照卡（包内重算 NEG1/NEG2） |
| `SOP-05_poses_main_list.sdf` | 主名单候选起始姿态子集（数据行与上传姿态逐字节一致，块以 M  END 收尾） |
| `SOP-05_report.md` | 最终报告（14_report.py 生成） |
| `SOP-05_robustness.json` | 稳健性 A–G 全量结果 |
| `figures/SOP-05_fig_druglikeness_funnel.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_druglikeness_funnel.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_gate_impact.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_gate_impact.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_gate_stack.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_gate_stack.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_herg_cross_model.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_herg_cross_model.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_list_overlap.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_list_overlap.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_md_batch.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_md_batch.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_pose_gate.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_pose_gate.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_robA_sigma.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_robA_sigma.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_robB_matrix.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_robB_matrix.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_robE_band.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_robE_band.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_robF_clusters.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |
| `figures/SOP-05_fig_robF_clusters.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现） |

脚本与辅助文件：

| 文件 | 说明 |
|---|---|
| `scripts/13_assemble.py` | 本包脚本副本（便于交付目录内追溯生成逻辑） |
| `scripts/14_report.py` | 本包脚本副本（便于交付目录内追溯生成逻辑） |
| `13_assemble.py` | 本包脚本副本（便于交付目录内追溯生成逻辑）（由后续步骤生成） |
| `14_report.py` | 本包脚本副本（便于交付目录内追溯生成逻辑）（由后续步骤生成） |
| `SOP-05_fig_pose_gate.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_pose_gate.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_druglikeness_funnel.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_druglikeness_funnel.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_gate_stack.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_gate_stack.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_list_overlap.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_list_overlap.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_gate_impact.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_gate_impact.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_robA_sigma.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_robA_sigma.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_robB_matrix.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_robB_matrix.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_robE_band.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_robE_band.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_robF_clusters.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_robF_clusters.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_md_batch.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_md_batch.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_herg_cross_model.png` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |
| `SOP-05_fig_herg_cross_model.svg` | 图（PNG+SVG，SVG 已去元数据以保证可复现）（由后续步骤生成） |

# 论文图（6 张）

> 由 `analysis/scripts/make_paper_fig1..6.py` 生成，全部 300 dpi、白底、碰撞扫描 0 处、
> 逐面板经视觉复查。共用样式与配色定义在 `analysis/scripts/_paper_common.py`。
>
> 图表目录按项目规范在 `figures/` 下；论文图放 `figures/paper/` 子目录，
> 与研究过程的 `figures/fig1..fig12` 分开。

| 图 | 文件 | 一句话主张 | 脚本 |
|---|---|---|---|
| Fig 1 | `fig1_conclusion_map.png` | 估计关联集中在选定的非冠脉心脏端点；其他测试端点在预设阈值下未见证据 | `make_paper_fig1.py` |
| Fig 2 | `fig2_replication.png` | 房颤 5 个数据集/定义、心衰 4 个队列方向一致；三个共享端点的芬兰点估计更负，但仅 1/3 差异显著 | `make_paper_fig2.py` |
| Fig 3 | `fig3_diffuse.png` | 600 位点仅 1 个共享因果变异、无通路承载、全库基因集扫描不超随机 | `make_paper_fig3.py` |
| Fig 4 | `fig4_modality.png` | 三篇超声一致保护，DXA 三份既不一致反向也不一致为零 | `make_paper_fig4.py` |
| Fig 5 | `fig5_population.png` | 人群层是预后标志而非共病状态；四个心血管结局无比其他更特异 | `make_paper_fig5.py` |
| Fig 6 | `fig6_robustness.png` | 预设循环性状与体成分模型未见实质衰减；阴性对照未见名义显著 | `make_paper_fig6.py` |

## 与旧图的对应关系

| 论文图 | 旧图 | 处理 |
|---|---|---|
| Fig 1 | `fig11` + `fig9` | 合并为一张四段式森林图：基质 / 负荷冠脉 / 传导 / 心房瓣膜心肌 |
| Fig 2 | `fig12` | 保留三面板；panel c 改为单色 + 印 p 值，避免与"芬兰＝蓝"冲突 |
| Fig 3 | `fig5` + `fig7` | 取共定位、通路分解、基因集扫描三个"弥散"层级 |
| Fig 4 | `fig10` | 两队列并列 + 互测；panel b 加 US/DXA 行标签 |
| Fig 5 | `fig1` | **叙事重写**（见下），面板改为纵向森林以便四个结局横向对照 |
| Fig 6 | `fig8` | 阴性对照保留，中介与体成分改为"衰减百分比"，去掉 LV 与中介的混排 |

降为补充材料：`fig2`、`fig3`（双侧 DEG 不收敛，该路线已否定）、`fig6`
（instrument→cell 富集是 `nearest_genes` 邻近注释上的描述性结果，不是因果定位）。

## Fig 5 的叙事改动（唯一一处旧结论翻转）

旧版把 NHANES 的"四项心血管结局无显著异质性"写成限制——"人群数据不支持心衰特异"。
模块 14 之后这句话反过来读：效应本来就横跨多个非冠脉心脏端点，**无异质性正是支持**。
Fig 5 的 panel c 与正文段落按新读法写，不是加一句注释。

正确表述：NHANES 显示低骨量是**预后**标志而非共病状态；四个终点之间没有差异；
两件事在遗传层定型之后才连成一个论证。

## 配色约定（勿改）

- 蓝 `#1f6fb4`：保护性且显著；在 Fig 2/4 中同时表示"芬兰队列"（同一批保护性估计）
- 灰 `#8a8a8a`：与零相容
- 橙 `#b45f1f`：零但功效不足（空心点）
- 黑 `#111111`：合并值 / 锚点
- 青 `#2b6f6a`：仅 Fig 2 panel c，表示"芬兰减非芬兰"这一差值量

## 复现

```powershell
.venv\Scripts\python.exe analysis\scripts\make_paper_fig1.py
.venv\Scripts\python.exe analysis\scripts\make_paper_fig2.py
.venv\Scripts\python.exe analysis\scripts\make_paper_fig3.py
.venv\Scripts\python.exe analysis\scripts\make_paper_fig4.py
.venv\Scripts\python.exe analysis\scripts\make_paper_fig5.py
.venv\Scripts\python.exe analysis\scripts\make_paper_fig6.py
```

每个脚本只读 `results/tables/` 下的既有表，不重算任何统计量；改图不会改数。

## 已知的图注必须说明的事

- **Fig 1**：18 个端点的例数与功效均在图中标出；"零但功效不足"用空心点区分，不能读作"确认无效应"。
- **Fig 2**：panel c 包含房颤、心衰和非缺血性心肌病；三行中只有一行 p<0.05，且心肌病定义不完全一致。
- **Fig 4**：DXA 三份中只有最大的一份显著为正，且 HERMES 中第三份为负——两面板必须同时看。
- **Fig 6**：panel c 中 `heart failure (positive control)` 是**阳性对照**，不是阴性对照；shared-locus 与 Steiger 检查在结果表/正文中报告，不绘制在当前 panel c。

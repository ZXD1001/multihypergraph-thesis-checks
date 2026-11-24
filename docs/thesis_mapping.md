# 与论文的对应关系

对照版本：2026-11-24 当前论文。以下采用源码文件名及语义标签定位，公式和页码随排版变化时，以标签为准。源码本身不在公开包中。

| 检查 | 对应位置 | 对应程度及用途 |
| --- | --- | --- |
| E01 | `data/chap03.tex`：`prop:ch03-composite`、`eq:ch03-prime-family-counts`、`fig:ch03-state-family-relations` | 核对小规模态族计数、包含/交关系和合数维零多项式；一般命题依赖正文证明 |
| E02 | `data/chap04.tex`：`eq:ch04-running-first-block`、`eq:ch04-running-second-block`、`eq:ch04-running-gap`、`tab:ch04-weight-samples` | 对应 `F=x₀x₁+2x₀x₁²` 的具体测量效应、等权谱和副本数 |
| E03 | `data/chap05.tex`：`ch05hq:lem:record-refinement`、`ch05hq:lem:choi`；`data/appendix_w3.tex`：`ch05hq:app:hiddenmix` | 为实际记录、混合分支及固定态恢复区别提供辅助反例；其中经典密钥例是辅助诊断，不是六方协议的全态模拟 |
| E04 | `data/chap05.tex`：`ch05hq:lem:record-refinement`、`ch05hq:lem:choi`；`data/appendix_w3.tex`：`ch05hq:lem:choi-proof` | 对粗粒化两分支信道作诊断，显示合并记录后不能直接套用细化分支的分类；不声称现稿逐项列出了这组数值 |
| E05 | `data/chap06.tex`：编码等距、门限隐私的定义与定理；`data/appendix_w4.tex`：`app:qss6-privacy-foundations` | 对两组固定相位表检验全部逻辑矩阵单位及禁止子集；任意相位函数的一般结论仍由解析证明给出 |
| E06 | `data/appendix_w4.tex`：`eq:qss6-evolution-five-evaluation-gates`、`eq:qss6-evolution-five-interpolation-gates`、`eq:qss6-evolution-five-residual-gates`、`app:qss6-precision` | 对应三组五方线路矩阵和5次/25次单位根物理相位的区别；模5湮灭与整数湮灭分别检验 |
| L05 | `data/chap05.tex`：`ch05hq:eq:evolution-budget`、`ch05hq:eq:evolution-bell`、`ch05hq:eq:evolution-top-measure`；`data/appendix_w3.tex`：`ch05hq:app:recovery-algorithm` | 核对规定的访问结构、标签分布与模块态表示、实际原始记录和恢复相位；不推导一般通信下界 |
| L346 | `data/chap06.tex`：`tab:qss6-scalarcomparison`；`data/appendix_w4.tex`：`app:qss6-algorithm-b`、`eq:app-qss6-infeasible-witness`、`eq:app-qss6-feasible-witness` | 四个固定编码的全部25位移和625目标；有限计算中的正负证书区分整数与模5关系 |
| CH34 | `data/chap03.tex` 及 `data/appendix_w1.tex` 的四个相位实例；`data/chap04.tex`：`eq:ch04-ordinary-cover-reduction`、`tab:ch04-weight-samples` | 逐层恢复参数、普通超图接受效应的还原、四行谱与副本数；具体标签见 `additional_checks.md` |
| ALG_A | `data/chap06.tex`：`eq:qss6-compilation-weights`、`eq:qss6-local-compiled-operator`、`eq:qss6-designedphase` | 对两组参数核对一般显式公式；相位门使用位移后更新的标签 |

E02 的带幂标签目标与 CH34 的普通超图目标 `F=x₀x₁` 不同。两者分别承担实例核对和共同子类还原，不合并成同一次目标态实验。

脚本中的系统、坐标通常从 0 开始编号；正文参与者编号与求值点的关系须依对应公式判断。L05 的“原始全局记录”包含全局读出组合，不能直接当作某个接收者在协议中实际收到的全部消息。390,625 条历史、15,625 种原始记录及每种记录的25个隐藏原像，是对正文读出映射的辅助枚举结果，并非正文另列的三个数值。

## 未纳入的材料

旧五方联合 POVM 证书对应另一套联盟和相位输入，与现稿六方十寄存器实例不一致，本次不随包发布。论文排版、语言替换、审稿回复生成、签字页面合并及历史 QA 工具也不属于本包的计算范围。

当前准备覆盖上述有限核验，未提供全文全部图表的自动重建程序，未实现一般参数的通用求解软件。已保存的数学函数和断言与原计算材料一致，基线来源的文件哈希列在 `reference_results/provenance.json`。

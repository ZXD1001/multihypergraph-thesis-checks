# 多重超图态与量子秘密共享：有限实例计算核验

本代码包配套张晓东的博士论文《多重超图态及其在量子秘密共享中的应用研究》，提供第三至第六章若干有限实例的可复核计算。版本为 **0.1.0**，发布日期为 **2026-11-24**。代码仓库为 [ZXD1001/multihypergraph-thesis-checks](https://github.com/ZXD1001/multihypergraph-thesis-checks)，版本下载见 [v0.1.0](https://github.com/ZXD1001/multihypergraph-thesis-checks/releases/tag/v0.1.0)。

10 个脚本覆盖参数反求、局域验证效应和样本数、六方三层共享实例，以及固定编码的逻辑相位门库。每项检查的参数和范围均在下表及说明文件中列明。包内不包含论文全文、审稿材料或个人签名文件，也不需要它们才能运行。

## 快速运行

已实测环境：**CPython 3.12.14、NumPy 2.3.5、macOS arm64**。统一入口要求 Python 3.12 或更新版本，依赖版本固定在 `requirements.txt`。计算过程不需要联网；首次安装依赖时可能需要访问软件包源。

在解压后的代码目录执行：

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python run_checks.py
```

Windows 可使用 `py -3.12 -m venv .venv`，随后将 `.venv/bin/python` 换为 `.venv\Scripts\python.exe`。Windows/Linux 尚未实测。

已有符合依赖的 Python 环境时，可直接执行 `python run_checks.py`。每次运行自动建立新的 `results/run_时间戳/` 目录。也可指定一个**尚不存在**的输出目录：

```sh
python run_checks.py --output-dir ../verification-run-01
python run_checks.py --list
python run_checks.py --checks E04 L346 --output-dir ../verification-selected-01
```

全部通过时，标准输出包含 10 行 `PASS`，退出码为 0。详细结果保存在每项检查的子目录，统一摘要为 `run_manifest.json`。任一脚本报错、超时或结果不符，入口均返回非零退出码，并继续运行其他选定项。不要使用 `-O` 或 `PYTHONOPTIMIZE` 跳过断言。

## 检查内容

| ID | 对象与参数 | 主要输出 |
| --- | --- | --- |
| E01 | 五组小维数态族；合数维 `d≤60` 的零多项式 | 态数、包含及交关系的有限检查 |
| E02 | 三个 qutrit 的两层局域测试 | 接受效应、谱重数 `16,10,1`；副本数 `598/600` |
| E03 | 五维 Choi、Clifford 混合及经典密钥实例 | 秩、边缘态、固定态与一般输入恢复的区别 |
| E04 | 作用于五维立方相位态的两分支信道 | Choi 秩 2；输出 Wigner 最小值 `−√5/25` |
| E05 | `(q,k,n,m)=(5,2,3,1)` 与 `(5,3,5,1)` | 共 500 个交叉算符检查及等距性残差 |
| E06 | 五方 `q=5` 编码 | 模证书、方向秩、3,125 个位移恒等式及三组线路矩阵 |
| L05 | 六方、三层、`q=5` | 64 个权限子集、8,000 个完整视图表示检查、存储 10/转交 4 |
| L346 | 四个固定相位、25 个位移、625 个归一化目标 | 62,500 个组合：2,500 个可行、60,000 个排除 |
| CH34 | 四个参数反求实例、两个普通超图测试块、四行样本表 | 66 个回代点、1,458 个矩阵位置、精确整数阈值 |
| ALG_A | `(q,k,m)=(5,2,2)` 和 `(7,3,1)` | 两组插值权、操作次序及 968 行完整相位等式 |

详细范围见 [E01—E06](docs/example_checks.md)、[六方实例及门库](docs/local_checks.md)、[参数恢复及算法 A](docs/additional_checks.md)。与论文的对应关系见 [章节映射](docs/thesis_mapping.md)。

## 如何判断运行结果

统一入口先执行各脚本，再逐字段比较数学结果与 `reference_results/` 中的归档基线。整数、布尔值、列表长度和结构须完全一致；浮点数采用绝对误差 `1e-12`、相对误差 `1e-10`。谱键中的 `-0.0` 与 `0.0` 归一化后比较。参考结果来自既有归档，运行不会自动用新结果覆盖它们。

E04 原脚本仅输出诊断量，统一入口另行验收其秩、完整谱、迹、Wigner 归一化及最小值；E03、E06 也有附加验收条件。检查不仅依赖脚本能否正常退出。

运行时间、Python 版本、来源说明和旧论文文件哈希不是数学输入。明确排除的顶层字段列在 `inputs/checks.json`；归档来源哈希和归一化规则见 `reference_results/provenance.json`。可选源码身份记录在统一入口中关闭，使复跑不依赖作者的本地目录。

## 结果适用范围

这些脚本核对指定参数下的计算，与论文的一般证明共同使用。它们没有重新证明全文定理，也不是任意参数下的通用门库求解器。六方实例的完整视图检查采用论文的模块分解、偏迹和秩论证，没有构建 `5^10` 维系统的全部密度矩阵；L346 保存的是三个指定目标的 75 份负证书及核行，不能将其描述为保存了全部 60,000 份独立负证书。各脚本说明列出了进一步的计算范围。

本包保留原有计算的数学函数和断言，移植涉及路径、输出和来源元数据。第三方依赖通过包管理器安装，不随本包分发。

## 版本、许可与引用

- 代码及随附说明采用 [MIT 许可](LICENSE)；标准许可文本可参阅 [SPDX](https://spdx.org/licenses/MIT.html)。
- 软件引用信息见 [CITATION.cff](CITATION.cff)，已列明仓库地址、版本号和发布日期。当前未分配 DOI。
- `SHA256SUMS.json` 记录冻结文件的 SHA-256；可运行 `python verify_integrity.py` 校验。
- 版本范围及后续维护见 [发布说明](docs/release_notes.md)。

# 测序数据核对工具（FastAPI）

两个入口，同一网页：

1. **双 VCF 变异核对**：同一参考、同一样本下，核对两套变异识别程序输出的单样本 VCF，
   按**归一化后的等位变异**列出两侧一致、仅单侧出现与剂量不符的结果，并保留每条结论的
   VCF 原始行号，供人工复查。
2. **单端 FASTQ 分子共识（UMI 去重）**：把 PCR 复制的读段按 UMI 归回原始分子，
   纠错后为每个合格分子组生成一条共识序列（FASTA 下载），避免重复证据放大。

## 双 VCF 变异核对：功能与判定规则

- **输入限制**
  - 单条染色体参考 FASTA，纯 A/C/G/T 序列，不超过 **1 MiB**；
  - 两份单样本 VCF（恰好 10 列、1 个样本），各不超过 **256 KiB**。
- **接受的变异**
  - A/C/G/T 构成的**单碱基替换**（SNV）；
  - 事件长度不超过 **20 bp** 的**纯插入 / 纯缺失**；
  - GT 必须是**非缺失、未分相的二倍体**（`0/1`、`1/1`，多 ALT 可为 `1/2`）。
  - 其他形态（MNV、符号化 ALT、`*`、分相 `|`、缺失 `.`、非二倍体等）一律判失败。
- **归一化**：对照参考序列对 indel 做左对齐（与 `bcftools norm` / `vt normalize` 等价），
  重复区域中位置写法不同的同一插入 / 缺失归为同一等位变异；SNV 位置不变。
- **剂量计算**：多 ALT 记录按 GT 分别计算每个 ALT 的 0/1/2 剂量，剂量为 0 的变异忽略；
  同一文件内同一归一变异的多条记录剂量累加，超过二倍体上限 2 即判定为**互相矛盾**。
- **整次失败**（返回 422，只显示错误、不展示任何局部结论）：
  - 任一记录 REF 与参考序列不符（含坐标越界、CHROM 不一致）；
  - 同一文件中同一归一变异存在互相矛盾的剂量；
  - 以及任何不满足上述输入规则的情形。
- **结果分组**（均按归一变异排序，并标注两侧原始行号与原始 `POS:REF>ALT` 写法）：
  两侧一致 / 仅 VCF A 或 B 出现 / 剂量不符。

## 单端 FASTQ 分子共识：功能与判定规则

- **输入限制**
  - 单端 FASTQ，不超过 **256 KiB**；
  - 读长 **40–100 bp**，序列仅允许 A/C/G/T/N；
  - 读段名（`@` 后首个空白前的 token）**末尾 6 个字符**为 A/C/G/T 构成的 **UMI**；
  - 质量为 **Phred+33**（可打印字符 `!`–`~`），且与序列等长。
  - 任何结构或质量错误都使**整批失败**（返回 422），不交付部分结果。
- **分组**：同一 UMI 的读段视为同一原始分子的 PCR 复制。
- **UMI 纠错**（方向性、防连锁）：
  - UMI `u` 是子组 `v` 的**合格父组**：`u` 与 `v` 仅差一位（Hamming 距离 1）、
    `count(u) ≥ 3` 且 `count(u) ≥ 3 × count(v)`（原始计数至少三条且为子组至少三倍）；
  - 子组仅有**唯一合格父组**、且该父组本身不是其他 UMI 的子组时，才被接收
    （读段并入父组）；
  - 被**多个父组争用**、或接收后会构成**连锁纠正**的子组一律**隔离**，
    不并入任何组、不生成共识。
- **组内评估**（纠错后的每个最终组）：
  - 组内**读长不一**，或任一位存在两种及以上 Q≥20 的碱基 → **冲突**；
  - 任一位没有至少**两条 Q≥20 读段**支持同一碱基 → **证据不足**；
  - 其余为**合格组**，按位取 Q≥20 的一致碱基作为分子共识。
- **输出**：页面显示各类组数（合格 / 冲突 / 证据不足 / 隔离）与需复查组的明细，
  只为合格组生成一条共识 FASTA 记录（头部标注 `reads=` 与纠错并入的 `corrected=`），
  可在线预览并下载。

## 目录结构

```
app/            FastAPI 应用与核心逻辑（bio.py、umi.py 仅用标准库）
  main.py       路由与上传校验
  bio.py        FASTA/VCF 解析、左对齐归一化、剂量与分组
  umi.py        FASTQ 解析、UMI 分组纠错、组内评估与分子共识
  templates/    上传页与两个结果页（Jinja2）
tests/          95 个单元 / 端到端测试
samples/        可直接试用的演示数据（demo.fa + 两份 VCF + reads.fastq）
scripts/        演示数据生成器（坐标 / 序列由程序自动计算）
Dockerfile, docker-compose.yml, requirements.txt
```

## 用 Docker Compose 启动（推荐）

```bash
docker compose up --build
```

浏览器入口：

- 上传页面：<http://localhost:8000/>
- 交互式 API 文档（Swagger UI）：<http://localhost:8000/docs>

停止：`docker compose down`。容器以只读根文件系统、`no-new-privileges` 运行，
内置 `/` 健康检查。

## 本地运行（无 Docker 时）

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## 试用

打开 <http://localhost:8000/>：

- **双 VCF 核对**：参考 `samples/demo.fa`，VCF A `samples/caller_a.vcf`，
  VCF B `samples/caller_b.vcf`。演示数据覆盖：SNV 一致、重复区内插入 / 缺失的
  两种等价写法归一、多 ALT 按 GT 拆分、`0/0` 剂量 0 被忽略、剂量不符、仅单侧出现。
- **分子共识**：FASTQ 选择 `samples/reads.fastq`。演示数据覆盖：合格组、
  子组并入父组、多父组争用隔离、连锁纠正隔离、碱基冲突、读长不一冲突、证据不足。

也可以走 JSON 接口：

```bash
# 双 VCF 核对
curl -F 'fasta=@samples/demo.fa' \
     -F 'vcf_a=@samples/caller_a.vcf' \
     -F 'vcf_b=@samples/caller_b.vcf' \
     http://localhost:8000/api/compare

# 单端 FASTQ 分子共识
curl -F 'fastq=@samples/reads.fastq' http://localhost:8000/api/consensus
```

## 测试

```bash
# 核心逻辑（仅标准库）
python3 -I -m unittest discover -s tests -t . -p 'test_bio.py' -v
python3 -I -m unittest discover -s tests -t . -p 'test_umi.py' -v

# 全部测试（Web 层需要先安装 requirements.txt）
python3 -I -m unittest discover -s tests -t . -v
```

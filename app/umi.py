"""单端 FASTQ 的 UMI 分组、纠错与分子共识（仅依赖标准库）。

处理流程：
1. 严格解析 4 行结构的单端 FASTQ（Phred+33 质量）：读长须为 40–100 bp，
   读段名（'@' 后首个空白前的 token）末尾 6 个字符须为 A/C/G/T 构成的 UMI。
   任何结构或质量错误都抛出 CheckError，整批失败，不交付部分结果。
2. 按 UMI 把读段分组：同一 UMI 的读段视为同一原始分子的 PCR 复制。
3. UMI 纠错（方向性、防连锁）：
   - UMI u 是子组 v 的合格父组：hamming(u, v) = 1、count(u) >= 3 且
     count(u) >= 3 * count(v)（原始计数至少三条且为子组至少三倍）；
   - 子组仅在有唯一合格父组、且该父组本身不是其他 UMI 的子组时才被接收
     （读段并入父组）；
   - 被多个合格父组争用、或接收它会构成连锁纠正的子组一律隔离，
     不并入任何组，也不生成共识。
4. 组内评估（对纠错后的每个最终组）：
   - 组内读长不一，或任一位存在两种及以上 Q>=20 的碱基 → 冲突；
   - 任一位没有任何碱基获得至少两条 Q>=20 读段支持 → 证据不足；
   - 其余为合格组，按位取 Q>=20 的一致碱基作为分子共识。
5. 只为合格组输出一条共识 FASTA 记录。
"""

from dataclasses import dataclass, field

from .bio import CheckError

# 上传体积上限（字节）
MAX_FASTQ_BYTES = 256 * 1024             # 256 KiB

# 读长范围（碱基）
MIN_READ_LEN = 40
MAX_READ_LEN = 100

# UMI 长度（读段名末尾）
UMI_LEN = 6

# UMI 纠错阈值
MIN_PARENT_COUNT = 3    # 父组原始计数下限
PARENT_RATIO = 3        # 父组计数至少为子组的倍数

# 共识阈值
MIN_BASE_Q = 20         # 高质量碱基的 Phred 阈值
MIN_SUPPORT = 2         # 每位支持同一碱基的高质量读段下限

_SEQ_BASES = frozenset("ACGTN")
_UMI_BASES = frozenset("ACGT")
_QUAL_MIN_ORD = 33      # '!'
_QUAL_MAX_ORD = 126     # '~'

CATEGORY_LABELS = {
    "qualified": "合格",
    "conflict": "冲突",
    "insufficient": "证据不足",
    "quarantined": "隔离",
}


@dataclass(frozen=True)
class Read:
    """一条单端读段。"""

    name: str
    umi: str
    seq: str
    quals: tuple[int, ...]  # 与 seq 等长的 Phred 分数


@dataclass
class GroupResult:
    """纠错后一个最终 UMI 组的评估结果。"""

    umi: str
    reads: int
    category: str                      # qualified / conflict / insufficient / quarantined
    consensus: str | None = None
    corrected: list[str] = field(default_factory=list)  # 被并入的子组 UMI
    note: str = ""

    @property
    def category_label(self) -> str:
        return CATEGORY_LABELS[self.category]


@dataclass
class ConsensusResult:
    """整批共识构建结果。"""

    total_reads: int
    raw_umis: int           # 原始 UMI 组数
    merged_children: int    # 被父组接收的子组数
    groups: list[GroupResult]  # 全部最终组（按 UMI 排序）
    fasta: str              # 合格组的分子共识 FASTA（可能为空字符串）

    @property
    def counts(self) -> dict[str, int]:
        out = {key: 0 for key in CATEGORY_LABELS}
        for g in self.groups:
            out[g.category] += 1
        return out

    @property
    def problem_groups(self) -> list[GroupResult]:
        """需复查的非合格组。"""
        return [g for g in self.groups if g.category != "qualified"]


# --------------------------------------------------------------------------- #
# FASTQ 解析
# --------------------------------------------------------------------------- #
def parse_fastq(data: bytes) -> list[Read]:
    """严格解析 4 行结构的单端 FASTQ；任何结构或质量错误都抛 CheckError。"""
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError as exc:
        raise CheckError(f"FASTQ 不是 ASCII 文本（字节 {exc.start} 处解码失败）。")

    lines = text.splitlines()
    if not lines:
        raise CheckError("FASTQ 为空文件。")
    if len(lines) % 4 != 0:
        raise CheckError(
            f"FASTQ 共 {len(lines)} 行，不是 4 的倍数，记录结构不完整，整批失败。"
        )

    reads: list[Read] = []
    for rec_start in range(0, len(lines), 4):
        rec_no = rec_start // 4 + 1
        name_line, seq_line, plus_line, qual_line = lines[rec_start:rec_start + 4]

        if not name_line.startswith("@"):
            raise CheckError(
                f"第 {rec_no} 条记录（第 {rec_start + 1} 行）：读段名行须以 '@' 开头。"
            )
        tokens = name_line[1:].split()
        if not tokens:
            raise CheckError(f"第 {rec_no} 条记录：'@' 后缺少读段名。")
        name = tokens[0]
        if len(name) < UMI_LEN:
            raise CheckError(
                f"第 {rec_no} 条记录：读段名 {name!r} 不足 {UMI_LEN} 个字符，"
                "末尾须带六碱基 UMI。"
            )
        umi = name[-UMI_LEN:]
        if any(ch not in _UMI_BASES for ch in umi):
            raise CheckError(
                f"第 {rec_no} 条记录：读段名 {name!r} 末尾的 {umi!r} 不是"
                " A/C/G/T 构成的六碱基 UMI。"
            )

        if not (MIN_READ_LEN <= len(seq_line) <= MAX_READ_LEN):
            raise CheckError(
                f"第 {rec_no} 条记录：读长 {len(seq_line)} bp 超出允许范围"
                f" {MIN_READ_LEN}–{MAX_READ_LEN} bp，整批失败。"
            )
        bad = next((ch for ch in seq_line if ch not in _SEQ_BASES), None)
        if bad is not None:
            raise CheckError(
                f"第 {rec_no} 条记录：序列含非法字符 {bad!r}，仅允许 A/C/G/T/N。"
            )

        if not plus_line.startswith("+"):
            raise CheckError(
                f"第 {rec_no} 条记录（第 {rec_start + 3} 行）：分隔行须以 '+' 开头。"
            )

        if len(qual_line) != len(seq_line):
            raise CheckError(
                f"第 {rec_no} 条记录：质量行长度 {len(qual_line)} 与序列长度"
                f" {len(seq_line)} 不一致，整批失败。"
            )
        bad_q = next(
            (ch for ch in qual_line if not (_QUAL_MIN_ORD <= ord(ch) <= _QUAL_MAX_ORD)),
            None,
        )
        if bad_q is not None:
            raise CheckError(
                f"第 {rec_no} 条记录：质量字符 {bad_q!r} 超出 Phred+33 可打印范围"
                "（'!'–'~'），整批失败。"
            )

        reads.append(
            Read(
                name=name,
                umi=umi,
                seq=seq_line,
                quals=tuple(ord(ch) - 33 for ch in qual_line),
            )
        )
    return reads


# --------------------------------------------------------------------------- #
# UMI 分组与纠错
# --------------------------------------------------------------------------- #
def _hamming1(a: str, b: str) -> bool:
    """两个等长 UMI 是否仅差一位。"""
    diff = 0
    for x, y in zip(a, b):
        if x != y:
            diff += 1
            if diff > 1:
                return False
    return diff == 1


def correct_umis(counts: dict[str, int]) -> tuple[dict[str, str], dict[str, str]]:
    """UMI 纠错。

    返回 (merge_into, quarantined)：
    - merge_into: {子组 UMI: 父组 UMI}，子组读段并入父组；
    - quarantined: {被隔离 UMI: 原因}，多父组争用或连锁纠正的子组。
    """
    def qualified_parents(v: str) -> list[str]:
        cv = counts[v]
        return [
            u for u in counts
            if u != v
            and counts[u] >= MIN_PARENT_COUNT
            and counts[u] >= PARENT_RATIO * cv
            and _hamming1(u, v)
        ]

    parents = {v: qualified_parents(v) for v in counts}
    merge_into: dict[str, str] = {}
    quarantined: dict[str, str] = {}
    for v in counts:
        ps = parents[v]
        if not ps:
            continue
        if len(ps) > 1:
            quarantined[v] = (
                f"与 {len(ps)} 个合格父组（{', '.join(sorted(ps))}）均仅差一位，"
                "多父组争用，予以隔离"
            )
            continue
        u = ps[0]
        if parents[u]:
            quarantined[v] = (
                f"唯一合格父组 {u} 本身还需向上纠错，接收将构成连锁纠正，予以隔离"
            )
            continue
        merge_into[v] = u
    return merge_into, quarantined


# --------------------------------------------------------------------------- #
# 组内评估与共识
# --------------------------------------------------------------------------- #
def _evaluate(reads: list[Read]) -> tuple[str, str | None]:
    """评估一个最终组，返回 (类别, 共识序列或 None)。

    冲突判定优先于证据不足：先扫描所有位点的高质量碱基是否一致，
    再检查每位的高质量支持读段数。
    """
    lengths = {len(r.seq) for r in reads}
    if len(lengths) != 1:
        return "conflict", None
    length = lengths.pop()

    hq_per_pos: list[list[str]] = []
    for i in range(length):
        hq = [r.seq[i] for r in reads if r.quals[i] >= MIN_BASE_Q]
        if len(set(hq)) > 1:
            return "conflict", None
        hq_per_pos.append(hq)

    consensus: list[str] = []
    for hq in hq_per_pos:
        if len(hq) < MIN_SUPPORT:
            return "insufficient", None
        consensus.append(hq[0])
    return "qualified", "".join(consensus)


def build_consensus(data: bytes) -> ConsensusResult:
    """整批流程：解析 → 分组 → UMI 纠错 → 组内评估 → 合格组共识 FASTA。"""
    reads = parse_fastq(data)

    raw: dict[str, list[Read]] = {}
    for r in reads:
        raw.setdefault(r.umi, []).append(r)
    counts = {umi: len(rs) for umi, rs in raw.items()}

    merge_into, quarantined = correct_umis(counts)

    final: dict[str, list[Read]] = {}
    corrected_by: dict[str, list[str]] = {}
    for umi, rs in raw.items():
        head = merge_into.get(umi, umi)
        final.setdefault(head, []).extend(rs)
        if head != umi:
            corrected_by.setdefault(head, []).append(umi)

    groups: list[GroupResult] = []
    fasta_lines: list[str] = []
    for umi in sorted(final):
        rs = final[umi]
        corrected = sorted(corrected_by.get(umi, []))
        if umi in quarantined:
            groups.append(
                GroupResult(umi, len(rs), "quarantined", None, corrected,
                            quarantined[umi])
            )
            continue
        category, consensus = _evaluate(rs)
        note = ""
        if category == "conflict":
            note = "组内读长不一，或某位存在 Q≥20 的不同碱基"
        elif category == "insufficient":
            note = "存在没有至少两条 Q≥20 读段支持同一碱基的位点"
        groups.append(GroupResult(umi, len(rs), category, consensus, corrected, note))
        if category == "qualified":
            header = f">UMI:{umi} reads={len(rs)}"
            if corrected:
                header += f" corrected={','.join(corrected)}"
            fasta_lines.append(header)
            fasta_lines.append(consensus)

    fasta = "\n".join(fasta_lines) + "\n" if fasta_lines else ""
    return ConsensusResult(
        total_reads=len(reads),
        raw_umis=len(raw),
        merged_children=len(merge_into),
        groups=groups,
        fasta=fasta,
    )

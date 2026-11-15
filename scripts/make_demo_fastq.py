"""生成演示用单端 FASTQ（samples/reads.fastq）。

运行：python3 scripts/make_demo_fastq.py
场景覆盖：
  合格组 / 子组并入父组（纠错）/ 多父组争用隔离 / 连锁纠正隔离 /
  碱基冲突 / 读长不一冲突 / 证据不足
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "samples" / "reads.fastq"

HI = "I"  # Phred+33 的 Q40
LO = "+"  # Phred+33 的 Q10
READ_LEN = 60


def make_seq(seed: int, length: int = READ_LEN) -> str:
    """确定性的伪随机 60 bp 序列。"""
    bases = "ACGT"
    return "".join(bases[(seed + 3 * i) % 4] for i in range(length))


def flip(base: str) -> str:
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]


def main() -> None:
    recs: list[tuple[str, str, str]] = []

    def add(umi: str, seq: str, qual: str) -> None:
        assert len(seq) == len(qual)
        recs.append((umi, seq, qual))

    # 1) 合格组：AACCGG 6 条相同读段（Q40）
    s1 = make_seq(1)
    for _ in range(6):
        add("AACCGG", s1, HI * READ_LEN)

    # 2) 子组：AACCGT 2 条（与 AACCGG 仅差一位，6 >= 3*2）→ 并入 AACCGG，共 8 条
    for _ in range(2):
        add("AACCGT", s1, HI * READ_LEN)

    # 3) 多父组争用：TTGGCC(5) 与 TTGGCA(5) 都和 TTGGCT(1) 仅差一位 → TTGGCT 隔离
    s2, s3 = make_seq(2), make_seq(3)
    for _ in range(5):
        add("TTGGCC", s2, HI * READ_LEN)
    for _ in range(5):
        add("TTGGCA", s3, HI * READ_LEN)
    add("TTGGCT", s2, HI * READ_LEN)

    # 4) 连锁纠正：GGGAAA(9) ← GGGAAT(3) ← GGGACT(1)
    #    GGGAAT 并入 GGGAAA；GGGACT 的唯一父组本身还要向上纠错 → 隔离
    s4 = make_seq(4)
    for _ in range(9):
        add("GGGAAA", s4, HI * READ_LEN)
    for _ in range(3):
        add("GGGAAT", s4, HI * READ_LEN)
    add("GGGACT", s4, HI * READ_LEN)

    # 5) 冲突（碱基）：CCCAAA 4 条，第 10 位两个不同碱基均为 Q40
    s5 = make_seq(5)
    s5m = s5[:9] + flip(s5[9]) + s5[10:]
    for _ in range(2):
        add("CCCAAA", s5, HI * READ_LEN)
        add("CCCAAA", s5m, HI * READ_LEN)

    # 6) 冲突（读长不一）：AAACCC 两条 60 bp + 一条 50 bp
    s6 = make_seq(6)
    add("AAACCC", s6, HI * READ_LEN)
    add("AAACCC", s6, HI * READ_LEN)
    add("AAACCC", s6[:50], HI * 50)

    # 7) 证据不足：TTTCGG 两条，其中一条全程 Q10，
    #    每位只剩一条 Q>=20 读段支持
    s7 = make_seq(7)
    add("TTTCGG", s7, HI * READ_LEN)
    add("TTTCGG", s7, LO * READ_LEN)

    lines: list[str] = []
    for i, (umi, seq, qual) in enumerate(recs, 1):
        lines += [f"@SIM{i:04d}:{umi}", seq, "+", qual]
    OUT.write_text("\n".join(lines) + "\n", encoding="ascii")
    print(f"wrote {len(recs)} reads -> samples/reads.fastq")


if __name__ == "__main__":
    main()

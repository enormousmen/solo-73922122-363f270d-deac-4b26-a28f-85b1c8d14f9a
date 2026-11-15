"""UMI 分子共识核心逻辑测试，仅依赖标准库：python -I -m unittest -v"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.bio import CheckError  # noqa: E402
from app.umi import (  # noqa: E402
    MAX_READ_LEN,
    MIN_READ_LEN,
    build_consensus,
    correct_umis,
    parse_fastq,
)

HI = "I"  # Phred+33 下的 Q40
LO = "+"  # Phred+33 下的 Q10
SEQ = "ACGT" * 10  # 40 bp，读长下限


def fq(records):
    """records: [(name, seq, qual), ...] -> FASTQ 字节串"""
    return "".join(f"@{n}\n{s}\n+\n{q}\n" for n, s, q in records).encode()


def rec(name, seq=SEQ, qual=None):
    """qual 为 None 时全 Q40；为单字符时重复至读长；否则视为完整质量串。"""
    if qual is None:
        q = HI * len(seq)
    elif len(qual) == 1:
        q = qual * len(seq)
    else:
        q = qual
    return (name, seq, q)


class ParseTests(unittest.TestCase):
    def test_happy_path_extracts_umi_and_quals(self):
        reads = parse_fastq(fq([rec("r1:AAAAAA"), rec("r2:CCCCGG")]))
        self.assertEqual(len(reads), 2)
        self.assertEqual(reads[0].umi, "AAAAAA")
        self.assertEqual(reads[1].umi, "CCCCGG")
        self.assertEqual(reads[0].seq, SEQ)
        self.assertEqual(reads[0].quals, (40,) * len(SEQ))

    def test_name_with_trailing_comment_uses_first_token(self):
        reads = parse_fastq(fq([rec("r1:AAAAAA 1:N:0:1")]))
        self.assertEqual(reads[0].name, "r1:AAAAAA")
        self.assertEqual(reads[0].umi, "AAAAAA")

    def test_crlf_line_endings_ok(self):
        data = ("@r1:AAAAAA\r\n" + SEQ + "\r\n+\r\n" + HI * len(SEQ) + "\r\n").encode()
        reads = parse_fastq(data)
        self.assertEqual(len(reads), 1)

    def test_min_and_max_read_length_accepted(self):
        short = "A" * MIN_READ_LEN
        long = "A" * MAX_READ_LEN
        reads = parse_fastq(fq([rec("r1:AAAAAA", short), rec("r2:AAAAAA", long)]))
        self.assertEqual(len(reads), 2)

    def test_empty_file_rejected(self):
        with self.assertRaises(CheckError):
            parse_fastq(b"")

    def test_non_ascii_rejected(self):
        data = fq([rec("r1:AAAAAA")])[:-1] + b"\xff\n"
        with self.assertRaises(CheckError):
            parse_fastq(data)

    def test_line_count_not_multiple_of_four_rejected(self):
        with self.assertRaises(CheckError) as ctx:
            parse_fastq(b"@r1:AAAAAA\n" + SEQ.encode() + b"\n+\n")
        self.assertIn("4 的倍数", str(ctx.exception))

    def test_name_line_must_start_with_at(self):
        with self.assertRaises(CheckError):
            parse_fastq(b"r1:AAAAAA\n" + SEQ.encode() + b"\n+\n" + (HI * len(SEQ)).encode() + b"\n")

    def test_empty_read_name_rejected(self):
        with self.assertRaises(CheckError):
            parse_fastq(b"@\n" + SEQ.encode() + b"\n+\n" + (HI * len(SEQ)).encode() + b"\n")

    def test_name_shorter_than_umi_rejected(self):
        with self.assertRaises(CheckError) as ctx:
            parse_fastq(fq([rec("ACGT")]))
        self.assertIn("六碱基 UMI", str(ctx.exception))

    def test_umi_must_be_acgt(self):
        for bad_name in ("r1:ACGTAN", "r1:acgtac", "r1:ACGT1C"):
            with self.subTest(name=bad_name):
                with self.assertRaises(CheckError):
                    parse_fastq(fq([rec(bad_name)]))

    def test_read_length_out_of_range_rejected(self):
        with self.assertRaises(CheckError) as ctx:
            parse_fastq(fq([rec("r1:AAAAAA", "A" * (MIN_READ_LEN - 1))]))
        self.assertIn("读长", str(ctx.exception))
        with self.assertRaises(CheckError):
            parse_fastq(fq([rec("r1:AAAAAA", "A" * (MAX_READ_LEN + 1))]))

    def test_bad_base_rejected(self):
        for bad in ("X" + SEQ[1:], SEQ[:5] + "a" + SEQ[6:]):
            with self.subTest(seq=bad):
                with self.assertRaises(CheckError):
                    parse_fastq(fq([rec("r1:AAAAAA", bad)]))

    def test_n_base_accepted(self):
        reads = parse_fastq(fq([rec("r1:AAAAAA", "N" * MIN_READ_LEN)]))
        self.assertEqual(reads[0].seq, "N" * MIN_READ_LEN)

    def test_plus_line_required(self):
        with self.assertRaises(CheckError) as ctx:
            parse_fastq(b"@r1:AAAAAA\n" + SEQ.encode() + b"\n-\n" + (HI * len(SEQ)).encode() + b"\n")
        self.assertIn("'+'", str(ctx.exception))

    def test_qual_length_mismatch_rejected(self):
        with self.assertRaises(CheckError) as ctx:
            parse_fastq(fq([("r1:AAAAAA", SEQ, HI * (len(SEQ) - 1))]))
        self.assertIn("不一致", str(ctx.exception))

    def test_qual_char_out_of_range_rejected(self):
        for bad_q in (" ", "\x7f"):
            with self.subTest(qual=repr(bad_q)):
                with self.assertRaises(CheckError) as ctx:
                    parse_fastq(fq([("r1:AAAAAA", SEQ, bad_q + HI * (len(SEQ) - 1))]))
                self.assertIn("Phred+33", str(ctx.exception))


class CorrectUmiTests(unittest.TestCase):
    def test_unique_parent_accepts_child(self):
        merge, quar = correct_umis({"AAAAAA": 6, "AAAAAT": 2})
        self.assertEqual(merge, {"AAAAAT": "AAAAAA"})
        self.assertEqual(quar, {})

    def test_ratio_boundary_inclusive(self):
        # 父组计数恰好为子组 3 倍即可接收
        merge, _ = correct_umis({"AAAAAA": 6, "AAAAAT": 2})
        self.assertIn("AAAAAT", merge)
        # 差一条不满足 3 倍：不接收
        merge, quar = correct_umis({"AAAAAA": 5, "AAAAAT": 2})
        self.assertEqual(merge, {})
        self.assertEqual(quar, {})

    def test_parent_needs_at_least_three_reads(self):
        # 父组原始计数须 >= 3（3 倍规则隐含，但显式验证边界）
        merge, _ = correct_umis({"AAAAAA": 3, "AAAAAT": 1})
        self.assertEqual(merge, {"AAAAAT": "AAAAAA"})

    def test_hamming_distance_two_not_merged(self):
        merge, quar = correct_umis({"AAAAAA": 10, "AAAATT": 2})
        self.assertEqual(merge, {})
        self.assertEqual(quar, {})

    def test_contested_child_quarantined(self):
        counts = {"AAAAAG": 10, "AAAACG": 1, "AAAACC": 10}
        # AAAAAG 与 AAAACC 相差两位，二者都与 AAAACG 差一位
        merge, quar = correct_umis(counts)
        self.assertEqual(merge, {})
        self.assertIn("AAAACG", quar)
        self.assertIn("多父组争用", quar["AAAACG"])

    def test_chained_correction_quarantined(self):
        counts = {"GGGAAA": 9, "GGGAAT": 3, "GGGACT": 1}
        merge, quar = correct_umis(counts)
        # GGGAAT 正常并入 GGGAAA；GGGACT 的唯一父组 GGGAAT 本身还要向上纠错 → 隔离
        self.assertEqual(merge, {"GGGAAT": "GGGAAA"})
        self.assertIn("GGGACT", quar)
        self.assertIn("连锁纠正", quar["GGGACT"])

    def test_quarantined_parent_does_not_accept_children(self):
        # AAAAAC 被 AAAAAA 与 AAAACC 争用而隔离；AAGAAC 的唯一合格父组是
        # AAAAAC（AAAAAA/AAAACC 与其差两位）→ 连锁纠正，一并隔离
        counts = {"AAAAAA": 10, "AAAACC": 10, "AAAAAC": 3, "AAGAAC": 1}
        merge, quar = correct_umis(counts)
        self.assertEqual(merge, {})
        self.assertIn("AAAAAC", quar)
        self.assertIn("AAGAAC", quar)

    def test_single_umi_no_correction(self):
        merge, quar = correct_umis({"AAAAAA": 5})
        self.assertEqual(merge, {})
        self.assertEqual(quar, {})


class EvaluateTests(unittest.TestCase):
    def test_qualified_group_consensus(self):
        result = build_consensus(fq([rec(f"r{i}:AAAAAA") for i in range(3)]))
        self.assertEqual(result.counts["qualified"], 1)
        (g,) = result.groups
        self.assertEqual(g.category, "qualified")
        self.assertEqual(g.consensus, SEQ)
        self.assertIn(f">UMI:AAAAAA reads=3\n{SEQ}\n", result.fasta)

    def test_conflict_by_mixed_read_lengths(self):
        data = fq([
            rec("r1:AAAAAA", "A" * 60),
            rec("r2:AAAAAA", "A" * 60),
            rec("r3:AAAAAA", "A" * 50),
        ])
        result = build_consensus(data)
        (g,) = result.groups
        self.assertEqual(g.category, "conflict")
        self.assertEqual(result.fasta, "")

    def test_conflict_by_highq_base_disagreement(self):
        seq_a = "A" * 40
        seq_c = "C" + "A" * 39
        data = fq([
            rec("r1:AAAAAA", seq_a), rec("r2:AAAAAA", seq_a),
            rec("r3:AAAAAA", seq_c), rec("r4:AAAAAA", seq_c),
        ])
        result = build_consensus(data)
        (g,) = result.groups
        self.assertEqual(g.category, "conflict")

    def test_lowq_difference_does_not_conflict(self):
        # 第三条读段在首位不同但 Q10：不冲突，共识取高质量一致碱基
        seq_c = "C" + SEQ[1:]
        qual = LO + HI * (len(SEQ) - 1)
        data = fq([
            rec("r1:AAAAAA"), rec("r2:AAAAAA"), rec("r3:AAAAAA", seq_c, qual),
        ])
        result = build_consensus(data)
        (g,) = result.groups
        self.assertEqual(g.category, "qualified")
        self.assertEqual(g.consensus, SEQ)

    def test_insufficient_when_fewer_than_two_highq_supporters(self):
        data = fq([rec("r1:AAAAAA"), rec("r2:AAAAAA", SEQ, LO * len(SEQ))])
        result = build_consensus(data)
        (g,) = result.groups
        self.assertEqual(g.category, "insufficient")
        self.assertEqual(result.fasta, "")

    def test_singleton_group_is_insufficient(self):
        result = build_consensus(fq([rec("r1:AAAAAA")]))
        (g,) = result.groups
        self.assertEqual(g.category, "insufficient")

    def test_conflict_takes_precedence_over_insufficient(self):
        # 首位只有一条高质量读段（证据不足），第 10 位两条高质量读段不一致（冲突）
        seq1 = SEQ
        seq2 = SEQ[:9] + "G" + SEQ[10:]
        qual2 = LO + HI * (len(SEQ) - 1)  # 首位 Q10，其余 Q40
        data = fq([rec("r1:AAAAAA", seq1), rec("r2:AAAAAA", seq2, qual2)])
        result = build_consensus(data)
        (g,) = result.groups
        self.assertEqual(g.category, "conflict")

    def test_two_identical_reads_can_be_qualified(self):
        result = build_consensus(fq([rec("r1:AAAAAA"), rec("r2:AAAAAA")]))
        (g,) = result.groups
        self.assertEqual(g.category, "qualified")


class BuildConsensusTests(unittest.TestCase):
    def test_full_scenario_counts_and_fasta(self):
        s1 = "ACGT" * 10
        s2 = "A" * 40
        s2m = "C" + "A" * 39
        s3 = "G" * 40
        s4 = "T" * 40
        records = []
        records += [rec(f"p{i}:AAAAAA", s1) for i in range(6)]   # 父组
        records += [rec(f"c{i}:AAAAAT", s1) for i in range(2)]   # 子组 → 并入
        records += [rec(f"x{i}:CCCCCC", s2 if i % 2 == 0 else s2m) for i in range(4)]  # 冲突
        records += [rec("y0:GGGGGG", s3)]                        # 证据不足
        records += [rec(f"q{i}:TTTTTT", s3) for i in range(10)]  # 合格
        records += [rec(f"w{i}:TTTTTG", s4) for i in range(10)]  # 合格
        records += [rec("z0:TTTTTA", s3)]                        # 多父组争用 → 隔离
        result = build_consensus(fq(records))

        self.assertEqual(result.total_reads, 6 + 2 + 4 + 1 + 10 + 10 + 1)
        self.assertEqual(result.raw_umis, 7)
        self.assertEqual(result.merged_children, 1)
        self.assertEqual(len(result.groups), 6)
        self.assertEqual(result.counts, {
            "qualified": 3, "conflict": 1, "insufficient": 1, "quarantined": 1,
        })

        by_umi = {g.umi: g for g in result.groups}
        self.assertEqual(by_umi["AAAAAA"].reads, 8)
        self.assertEqual(by_umi["AAAAAA"].corrected, ["AAAAAT"])
        self.assertEqual(by_umi["AAAAAA"].consensus, s1)
        self.assertEqual(by_umi["CCCCCC"].category, "conflict")
        self.assertEqual(by_umi["GGGGGG"].category, "insufficient")
        self.assertEqual(by_umi["TTTTTA"].category, "quarantined")
        self.assertIn("多父组争用", by_umi["TTTTTA"].note)

        # FASTA 只含合格组，按 UMI 排序，头部标注读段数与纠错来源
        lines = result.fasta.splitlines()
        self.assertEqual(lines, [
            ">UMI:AAAAAA reads=8 corrected=AAAAAT",
            s1,
            ">UMI:TTTTTG reads=10",
            s4,
            ">UMI:TTTTTT reads=10",
            s3,
        ])

    def test_no_qualified_groups_gives_empty_fasta(self):
        result = build_consensus(fq([rec("r1:AAAAAA")]))
        self.assertEqual(result.fasta, "")
        self.assertEqual(result.counts["qualified"], 0)

    def test_problem_groups_excludes_qualified(self):
        records = [rec(f"r{i}:AAAAAA") for i in range(3)]
        records += [rec("s1:CCCCCC")]
        result = build_consensus(fq(records))
        self.assertEqual([g.umi for g in result.problem_groups], ["CCCCCC"])


if __name__ == "__main__":
    unittest.main()

"""FastAPI 入口：双 VCF 变异核对 + 单端 FASTQ UMI 分子共识。"""

from pathlib import Path

from fastapi import FastAPI, File, Form, Request, UploadFile
from fastapi.responses import HTMLResponse, JSONResponse, Response
from fastapi.templating import Jinja2Templates

from .bio import (
    MAX_FASTA_BYTES,
    MAX_VCF_BYTES,
    CheckError,
    compare,
    parse_fasta,
    parse_vcf,
)
from .umi import MAX_FASTQ_BYTES, build_consensus

app = FastAPI(
    title="测序数据核对工具",
    description="双 VCF 变异核对，以及单端 FASTQ 的 UMI 分组纠错与分子共识。",
)

templates = Jinja2Templates(directory=str(Path(__file__).parent / "templates"))


def _index_context(**extra):
    ctx = {
        "fasta_limit_kib": MAX_FASTA_BYTES // 1024,
        "vcf_limit_kib": MAX_VCF_BYTES // 1024,
        "fastq_limit_kib": MAX_FASTQ_BYTES // 1024,
    }
    ctx.update(extra)
    return ctx


def _read_limited(file: UploadFile, limit: int, label: str) -> bytes:
    data = file.file.read(limit + 1)
    if len(data) > limit:
        raise CheckError(f"{label} 超过大小上限 {limit} 字节（{limit // 1024} KiB）。")
    if not data:
        raise CheckError(f"{label} 为空文件。")
    return data


def _run_compare(fasta_bytes: bytes, vcf_a_bytes: bytes, vcf_b_bytes: bytes):
    refobj = parse_fasta(fasta_bytes)
    calls_a, sample_a = parse_vcf(vcf_a_bytes, "VCF A", refobj)
    calls_b, sample_b = parse_vcf(vcf_b_bytes, "VCF B", refobj)
    groups = compare(refobj, calls_a, calls_b)
    return refobj, sample_a, sample_b, groups


@app.get("/", response_class=HTMLResponse)
def index(request: Request):
    return templates.TemplateResponse(request, "index.html", _index_context())


# --------------------------------------------------------------------------- #
# 双 VCF 变异核对
# --------------------------------------------------------------------------- #
@app.post("/compare", response_class=HTMLResponse)
def compare_page(
    request: Request,
    fasta: UploadFile = File(...),
    vcf_a: UploadFile = File(...),
    vcf_b: UploadFile = File(...),
):
    try:
        fasta_bytes = _read_limited(fasta, MAX_FASTA_BYTES, "参考 FASTA")
        vcf_a_bytes = _read_limited(vcf_a, MAX_VCF_BYTES, "VCF A")
        vcf_b_bytes = _read_limited(vcf_b, MAX_VCF_BYTES, "VCF B")
        refobj, sample_a, sample_b, groups = _run_compare(
            fasta_bytes, vcf_a_bytes, vcf_b_bytes
        )
    except CheckError as exc:
        # 整次核对失败：只展示错误，不展示任何局部结论
        return templates.TemplateResponse(
            request,
            "index.html",
            _index_context(error=str(exc)),
            status_code=422,
        )

    return templates.TemplateResponse(
        request,
        "result.html",
        {
            "ref_name": refobj.name,
            "ref_length": len(refobj.sequence),
            "sample_a": sample_a,
            "sample_b": sample_b,
            "groups": groups,
        },
    )


@app.post("/api/compare")
def compare_api(
    fasta: UploadFile = File(...),
    vcf_a: UploadFile = File(...),
    vcf_b: UploadFile = File(...),
):
    """JSON 接口，便于脚本化调用；失败返回 422 与错误信息。"""
    try:
        fasta_bytes = _read_limited(fasta, MAX_FASTA_BYTES, "参考 FASTA")
        vcf_a_bytes = _read_limited(vcf_a, MAX_VCF_BYTES, "VCF A")
        vcf_b_bytes = _read_limited(vcf_b, MAX_VCF_BYTES, "VCF B")
        refobj, sample_a, sample_b, groups = _run_compare(
            fasta_bytes, vcf_a_bytes, vcf_b_bytes
        )
    except CheckError as exc:
        return JSONResponse(status_code=422, content={"ok": False, "error": str(exc)})

    def dump_entry(entry):
        def dump_call(call):
            if call is None:
                return None
            return {
                "dosage": call.dosage,
                "line_numbers": call.line_numbers,
                "raw": call.raw_descriptions,
            }

        v = entry["variant"]
        return {
            "variant": {
                "pos": v.pos,
                "ref": v.ref,
                "alt": v.alt,
                "kind": v.kind,
                "label": v.label,
            },
            "left": dump_call(entry["left"]),
            "right": dump_call(entry["right"]),
        }

    return {
        "ok": True,
        "reference": {"name": refobj.name, "length": len(refobj.sequence)},
        "samples": {"left": sample_a, "right": sample_b},
        "result": {key: [dump_entry(e) for e in rows] for key, rows in groups.items()},
    }


# --------------------------------------------------------------------------- #
# 单端 FASTQ UMI 分子共识
# --------------------------------------------------------------------------- #
@app.post("/consensus", response_class=HTMLResponse)
def consensus_page(request: Request, fastq: UploadFile = File(...)):
    try:
        fastq_bytes = _read_limited(fastq, MAX_FASTQ_BYTES, "FASTQ")
        result = build_consensus(fastq_bytes)
    except CheckError as exc:
        # 整批失败：只展示错误，不交付任何部分结果
        return templates.TemplateResponse(
            request,
            "index.html",
            _index_context(consensus_error=str(exc)),
            status_code=422,
        )

    return templates.TemplateResponse(
        request,
        "consensus_result.html",
        {
            "total_reads": result.total_reads,
            "raw_umis": result.raw_umis,
            "merged_children": result.merged_children,
            "final_groups": len(result.groups),
            "counts": result.counts,
            "problem_groups": result.problem_groups,
            "fasta": result.fasta,
        },
    )


@app.post("/consensus/download")
def consensus_download(fasta: str = Form(...)):
    """把结果页内嵌的共识 FASTA 作为附件发回（服务本身不保存任何状态）。"""
    return Response(
        content=fasta,
        media_type="text/x-fasta; charset=utf-8",
        headers={"Content-Disposition": 'attachment; filename="consensus.fa"'},
    )


@app.post("/api/consensus")
def consensus_api(fastq: UploadFile = File(...)):
    """JSON 接口，便于脚本化调用；失败返回 422 与错误信息。"""
    try:
        fastq_bytes = _read_limited(fastq, MAX_FASTQ_BYTES, "FASTQ")
        result = build_consensus(fastq_bytes)
    except CheckError as exc:
        return JSONResponse(status_code=422, content={"ok": False, "error": str(exc)})

    return {
        "ok": True,
        "total_reads": result.total_reads,
        "raw_umis": result.raw_umis,
        "merged_children": result.merged_children,
        "counts": result.counts,
        "groups": [
            {
                "umi": g.umi,
                "reads": g.reads,
                "category": g.category,
                "corrected": g.corrected,
                "note": g.note,
                "consensus": g.consensus,
            }
            for g in result.groups
        ],
        "fasta": result.fasta,
    }

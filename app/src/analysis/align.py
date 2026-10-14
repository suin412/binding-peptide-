"""원본 vs 변이체 서열 정렬 (Biopython PairwiseAligner, BLOSUM62)."""
from __future__ import annotations

import html
from dataclasses import dataclass

from Bio.Align import PairwiseAligner, substitution_matrices


@dataclass
class AlignmentResult:
    aligned_a: str
    aligned_b: str
    midline: str
    identity: float                 # 일치 칸 수 / 정렬 길이 (0~1)
    matches: int
    diff_positions_a: list[int]     # a 기준 1-based, 차이(불일치/gap 대응) 위치
    diff_positions_b: list[int]
    score: float


def _aligner() -> PairwiseAligner:
    al = PairwiseAligner()
    al.mode = "global"
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -10
    al.extend_gap_score = -0.5
    return al


def align_sequences(a: str, b: str) -> AlignmentResult:
    aln = _aligner().align(a, b)[0]
    blocks_a, blocks_b = aln.aligned
    row_a, row_b = [], []
    ia = ib = 0
    for (sa, ea), (sb, eb) in zip(blocks_a, blocks_b):
        # 블록 사이의 gap 구간
        gap_a, gap_b = sa - ia, sb - ib
        if gap_a:
            row_a.append(a[ia:sa]); row_b.append("-" * gap_a)
        if gap_b:
            row_a.append("-" * gap_b); row_b.append(b[ib:sb])
        row_a.append(a[sa:ea]); row_b.append(b[sb:eb])
        ia, ib = ea, eb
    if ia < len(a):
        row_a.append(a[ia:]); row_b.append("-" * (len(a) - ia))
    if ib < len(b):
        row_a.append("-" * (len(b) - ib)); row_b.append(b[ib:])
    aa, bb = "".join(row_a), "".join(row_b)

    mid, matches = [], 0
    diff_a, diff_b = [], []
    pos_a = pos_b = 0
    for x, y in zip(aa, bb):
        if x != "-":
            pos_a += 1
        if y != "-":
            pos_b += 1
        if x == y:
            mid.append("|"); matches += 1
        else:
            mid.append(" ")
            if x != "-":
                diff_a.append(pos_a)
            if y != "-":
                diff_b.append(pos_b)
    length = len(aa) or 1
    return AlignmentResult(aa, bb, "".join(mid), matches / length, matches, diff_a, diff_b, float(aln.score))


def render_html(res: AlignmentResult, name_a: str = "원본", name_b: str = "변이체", width: int = 60) -> str:
    """차이 나는 칸을 강조한 모노스페이스 정렬 HTML."""
    def cell(ch: str, other: str) -> str:
        ch_e = html.escape(ch)
        if ch == other:
            return ch_e
        color = "#ffd54f" if "-" not in (ch, other) else "#ef9a9a"
        return f'<span style="background:{color};color:#000">{ch_e}</span>'

    label_w = max(len(name_a), len(name_b), 4)
    blocks = []
    for start in range(0, len(res.aligned_a), width):
        sa = res.aligned_a[start:start + width]
        sb = res.aligned_b[start:start + width]
        sm = res.midline[start:start + width]
        line_a = "".join(cell(x, y) for x, y in zip(sa, sb))
        line_b = "".join(cell(y, x) for x, y in zip(sa, sb))
        blocks.append(
            f"{html.escape(name_a[:label_w]).ljust(label_w)} {line_a}\n"
            f"{' ' * label_w} {html.escape(sm)}\n"
            f"{html.escape(name_b[:label_w]).ljust(label_w)} {line_b}\n"
        )
    return '<pre style="font-family:monospace;line-height:1.25;white-space:pre">' + "\n".join(blocks) + "</pre>"

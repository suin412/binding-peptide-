"""서열 입력 정리·검증, FASTA 파싱, 3문자 표기 변환."""
from __future__ import annotations

import io
import re

from Bio import SeqIO

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# 비표준 잔기 → 가장 가까운 표준 잔기(단순화용). 구조 예측에서는 원래 화학 구조가 무시된다.
NONSTANDARD_TO_STANDARD = {
    "AIB": "A", "ABU": "A", "DAP": "A", "NVA": "V", "NLE": "L",
    "ORN": "K", "DAB": "K", "PYL": "K", "CIT": "R", "HYP": "P", "SEC": "C",
}

MIN_LEN = 2


class SequenceError(ValueError):
    """서열 검증 실패. 메시지에 문제 위치를 포함한다."""


def clean_sequence(text: str) -> str:
    """공백·줄바꿈·숫자를 제거하고 대문자로 바꾼다 (FASTA 헤더 줄은 무시)."""
    lines = [ln for ln in (text or "").splitlines() if not ln.lstrip().startswith(">")]
    raw = "".join(lines)
    return re.sub(r"[\s\d]+", "", raw).upper()


def validate_sequence(text: str, min_len: int = MIN_LEN) -> str:
    """표준 20 아미노산만 허용. 통과하면 정리된 서열을, 아니면 SequenceError를 던진다."""
    seq = clean_sequence(text)
    if not seq:
        raise SequenceError("서열이 비어 있습니다.")
    bad = [(i + 1, ch) for i, ch in enumerate(seq) if ch not in STANDARD_AA]
    if bad:
        shown = ", ".join(f"{pos}번째 '{ch}'" for pos, ch in bad[:8])
        more = f" 외 {len(bad) - 8}개" if len(bad) > 8 else ""
        raise SequenceError(
            f"표준 20 아미노산(1문자 표기)이 아닌 문자가 있습니다: {shown}{more}. "
            "비표준 잔기는 가장 가까운 표준 잔기로 바꾸고 원문 표기를 메모에 남겨 주세요."
        )
    if len(seq) < min_len:
        raise SequenceError(f"서열이 너무 짧습니다 (최소 {min_len} aa).")
    return seq


def is_valid(text: str) -> bool:
    try:
        validate_sequence(text)
        return True
    except SequenceError:
        return False


def parse_fasta(text: str) -> list[tuple[str, str]]:
    """FASTA 텍스트 → [(이름, 서열)]. '>' 헤더가 없으면 서열 하나로 취급한다. 검증은 하지 않는다."""
    if ">" not in (text or ""):
        seq = clean_sequence(text)
        return [("입력 서열", seq)] if seq else []
    out = []
    for rec in SeqIO.parse(io.StringIO(text), "fasta"):
        seq = re.sub(r"[\s\d]+", "", str(rec.seq)).upper()
        out.append((rec.id or "서열", seq))
    return out


def to_fasta(items: list[tuple[str, str]]) -> str:
    return "".join(f">{name}\n{seq}\n" for name, seq in items)

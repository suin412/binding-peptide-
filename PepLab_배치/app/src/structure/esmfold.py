"""ESMFold 공개 API로 구조 예측. 길이 제한·비표준 문자 차단·캐시."""
from __future__ import annotations

from ..sequence import SequenceError, validate_sequence
from ..sources.http import ApiError, post_text
from .pdbutil import normalize_plddt

URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
MAX_LEN = 400
SHORT_PEPTIDE_LEN = 30


class StructureLimitError(ValueError):
    """ESMFold 호출 전에 걸러지는 입력 문제 (길이 초과 등)."""


def check_sequence(seq: str) -> str:
    try:
        seq = validate_sequence(seq)
    except SequenceError as exc:
        raise StructureLimitError(str(exc)) from exc
    if len(seq) > MAX_LEN:
        raise StructureLimitError(
            f"서열 길이가 {len(seq)} aa입니다. ESMFold 공개 API는 약 {MAX_LEN} aa 이하만 지원하므로 "
            "호출하지 않았습니다. 서열을 잘라서 사용해 주세요."
        )
    return seq


def predict(seq: str, *, post=post_text) -> str:
    """서열 → PDB 텍스트 (B-factor는 0~100 pLDDT로 정규화). post는 테스트용 주입 지점."""
    seq = check_sequence(seq)
    text = post(URL, seq, service="ESMFold")
    if "ATOM" not in text:
        raise ApiError("ESMFold가 구조를 반환하지 않았습니다. 잠시 후 다시 시도해 주세요.")
    return normalize_plddt(text)

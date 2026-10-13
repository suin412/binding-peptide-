"""서열 물성 계산 (Biopython ProtParam) + 소수성 프로파일·조성·소수성 모멘트."""
from __future__ import annotations

import math

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParamData import kd

from ..sequence import STANDARD_AA, validate_sequence


def compute_properties(seq: str) -> dict:
    """길이, 분자량, pI, 순전하(pH 7.4), GRAVY, 불안정성 지수, 방향족성.

    순전하는 Biopython의 charge_at_pH(7.4) 값을 사용한다 (측쇄·말단기 pKa 기반 추정).
    """
    seq = validate_sequence(seq, min_len=1)
    pa = ProteinAnalysis(seq)
    return {
        "length": len(seq),
        "mw": round(pa.molecular_weight(), 2),
        "pI": round(pa.isoelectric_point(), 2),
        "net_charge_pH7.4": round(pa.charge_at_pH(7.4), 2),
        "gravy": round(pa.gravy(), 3),
        "instability_index": round(pa.instability_index(), 1) if len(seq) >= 2 else float("nan"),
        "aromaticity": round(pa.aromaticity(), 3),
    }


def net_charge(seq: str, ph: float = 7.4) -> float:
    return ProteinAnalysis(seq).charge_at_pH(ph)


def hydropathy_profile(seq: str, window: int = 7) -> list[float]:
    """Kyte-Doolittle 슬라이딩 윈도우 평균 (윈도우는 홀수로 보정, 서열보다 길면 축소)."""
    window = max(1, min(window, len(seq)))
    if window % 2 == 0:
        window = max(1, window - 1)
    values = [kd[a] for a in seq]
    half = window // 2
    out = []
    for i in range(len(seq)):
        lo, hi = max(0, i - half), min(len(seq), i + half + 1)
        out.append(sum(values[lo:hi]) / (hi - lo))
    return out


def composition(seq: str) -> dict[str, float]:
    n = len(seq)
    return {aa: seq.count(aa) / n * 100 for aa in STANDARD_AA}


def hydrophobic_moment(seq: str, angle_deg: float = 100.0) -> float:
    """잔기당 평균 소수성 모멘트 (KD 척도를 ±4.5로 정규화, 0~1 근처). 헬릭스 양친매성 지표."""
    if not seq:
        return 0.0
    sin_sum = sum(kd[a] / 4.5 * math.sin(math.radians(angle_deg * i)) for i, a in enumerate(seq))
    cos_sum = sum(kd[a] / 4.5 * math.cos(math.radians(angle_deg * i)) for i, a in enumerate(seq))
    return math.hypot(sin_sum, cos_sum) / len(seq)

"""표면 결합 성향 점수 (학습용 경험적 모델).

**이것은 결합 친화도(Kd)의 예측이 아니다.** 표면 화학별로 문헌에서 알려진 *정성적* 경향
(예: 소수성 고분자에는 방향족·소수성 잔기, 티타늄·실리카에는 양전하 잔기, 하이드록시아파타이트에는 산성 잔기)을
잔기별 가중치로 단순화해, 서열끼리 "이 표면에 붙을 성향이 상대적으로 큰가"를 줄 세우는 용도다.
pH, 이온 세기, 표면 상태, 펩타이드의 입체 구조는 고려하지 않는다. 가중치는 가정이므로 문헌 값과 비교해 보정해야 한다.

점수는 "같은 길이의 무작위 서열(Swiss-Prot 평균 조성)보다 높을 확률"(백분위, 0~100)로 표시한다.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from ..sequence import STANDARD_AA, validate_sequence

# Swiss-Prot 평균 아미노산 조성(%, 근사) — 무작위 서열 기준선
BACKGROUND = {
    "A": 8.25, "R": 5.53, "N": 4.06, "D": 5.45, "C": 1.37, "Q": 3.93, "E": 6.75, "G": 7.07, "H": 2.27, "I": 5.96,
    "L": 9.66, "K": 5.84, "M": 2.42, "F": 3.86, "P": 4.70, "S": 6.56, "T": 5.34, "W": 1.08, "Y": 2.92, "V": 6.87,
}


@dataclass(frozen=True)
class SurfaceClass:
    key: str
    label: str
    weights: dict[str, float]
    rationale: str


CLASSES: dict[str, SurfaceClass] = {c.key: c for c in [
    SurfaceClass(
        "polymer", "소수성 고분자 (PET, PS, PE, PP, PMMA …)",
        {"W": 1.0, "F": 0.85, "Y": 0.7, "H": 0.6, "M": 0.5, "L": 0.4, "I": 0.4, "R": 0.35, "V": 0.25, "C": 0.2,
         "K": 0.1, "A": 0.05, "G": -0.1, "S": -0.1, "T": -0.1, "Q": -0.1, "N": -0.2, "E": -0.4, "D": -0.5},
        "소수성·π 상호작용: 방향족(W·F·Y·H)과 큰 소수성 잔기가 유리하고 산성 잔기는 불리하다고 가정. "
        "플라스틱 결합 펩타이드를 계산으로 설계한 연구에서도 W·H·R이 많이 나타난다.",
    ),
    SurfaceClass(
        "titanium", "티타늄·금속 산화물 (Ti/TiO₂, 스테인리스강)",
        {"R": 1.0, "K": 0.8, "H": 0.6, "S": 0.35, "T": 0.35, "Y": 0.3, "N": 0.2, "D": 0.2, "W": 0.2, "E": 0.1},
        "중성 pH의 산화막 표면은 음전하를 띠어 양전하 잔기(R·K·H)와 수산기(S·T)가 유리하다고 가정. "
        "알려진 Ti 결합 펩타이드(RKLPDA 계열 등)가 R·K를 포함한다.",
    ),
    SurfaceClass(
        "noble", "귀금속 (Au, Ag, Pt)",
        {"H": 0.8, "C": 0.8, "M": 0.6, "W": 0.6, "Y": 0.5, "R": 0.5, "K": 0.3, "F": 0.2, "S": 0.15},
        "황(C·M)·이미다졸(H)·방향족(W·Y) 잔기가 금속 표면에 배위·흡착한다고 가정.",
    ),
    SurfaceClass(
        "silica", "실리카·유리 (SiO₂)",
        {"R": 1.0, "K": 0.9, "S": 0.5, "N": 0.3, "T": 0.3, "Y": 0.3, "H": 0.3, "D": -0.3, "E": -0.3},
        "음전하 실란올 표면에 양전하 잔기(R·K)와 수산기·아마이드 잔기가 유리하고 산성 잔기는 불리하다고 가정.",
    ),
    SurfaceClass(
        "apatite", "하이드록시아파타이트·뼈 (칼슘 인산염)",
        {"D": 1.0, "E": 0.9, "S": 0.4, "T": 0.15, "K": 0.15, "R": 0.1},
        "표면 Ca²⁺와 산성 잔기(D·E, 인산화 Ser)의 상호작용이 유리하다고 가정.",
    ),
    SurfaceClass(
        "carbon", "탄소 소재 (그래핀, 탄소나노튜브)",
        {"W": 1.0, "F": 0.9, "Y": 0.8, "H": 0.6, "R": 0.25, "L": 0.15, "M": 0.15, "D": -0.3, "E": -0.3},
        "π–π 스태킹: 방향족 잔기(W·F·Y·H)가 유리하다고 가정.",
    ),
]}

# 앱의 결합 대상 이름(binding.MATERIALS의 키) → 표면 분류
_MATERIAL_TO_CLASS = {
    "티타늄": "titanium", "스테인리스강": "titanium", "산화아연": "titanium", "산화철": "titanium",
    "금": "noble", "은": "noble", "백금": "noble",
    "실리카": "silica", "실리콘": "silica",
    "하이드록시아파타이트": "apatite", "인산칼슘": "apatite",
    "그래핀": "carbon", "탄소나노튜브": "carbon",
    "폴리스티렌": "polymer", "PET": "polymer", "폴리에틸렌": "polymer", "폴리프로필렌": "polymer", "PMMA": "polymer",
    "PDMS": "polymer", "PTFE": "polymer", "폴리염화비닐": "polymer", "폴리락트산": "polymer", "폴리카보네이트": "polymer",
    "나일론": "polymer", "폴리우레탄": "polymer", "마이크로플라스틱": "polymer",
}


def class_for_material(name: str) -> str | None:
    """결합 대상 이름(예: 'PET (폴리에틸렌 테레프탈레이트)')에서 표면 분류를 추정한다. 모르면 None."""
    name = name or ""
    for key, cls in _MATERIAL_TO_CLASS.items():
        # 한 글자 키(금, 은)는 "금속" 같은 일반 단어에 걸리지 않도록 "금", "금 (Au)"처럼 단독으로 쓰인 경우만 인정한다
        if len(key) == 1:
            if name == key or name.startswith((key + " ", key + "(")):
                return cls
        elif key in name:
            return cls
    return None


@dataclass
class BindingScore:
    class_key: str
    raw: float                  # 잔기당 평균 가중치
    z: float                    # 무작위 서열 대비 표준점수
    percentile: float           # 같은 길이 무작위 서열보다 높을 확률 (0~100) - 강한 후보끼리는 100 근처로 포화한다
    scaled: float               # 잔기당 평균 가중치를 그 표면의 (최소~최대 가중치) 구간에 대응시킨 0~100 - 후보끼리 구분할 때 사용
    contributions: list[float]  # 잔기별 가중치


def _background_moments(weights: dict[str, float]) -> tuple[float, float]:
    total = sum(BACKGROUND.values())
    f = {aa: BACKGROUND[aa] / total for aa in STANDARD_AA}
    mean = sum(f[aa] * weights.get(aa, 0.0) for aa in STANDARD_AA)
    var = sum(f[aa] * weights.get(aa, 0.0) ** 2 for aa in STANDARD_AA) - mean ** 2
    return mean, math.sqrt(max(var, 1e-12))


def score_sequence(seq: str, class_key: str) -> BindingScore:
    seq = validate_sequence(seq, min_len=1)
    cls = CLASSES[class_key]
    contrib = [cls.weights.get(aa, 0.0) for aa in seq]
    raw = sum(contrib) / len(seq)
    mean, sd = _background_moments(cls.weights)
    z = (raw - mean) / (sd / math.sqrt(len(seq)))
    percentile = 50.0 * (1.0 + math.erf(z / math.sqrt(2.0)))
    lo, hi = min(cls.weights.values()), max(cls.weights.values())
    scaled = 100.0 * (raw - lo) / (hi - lo) if hi > lo else 0.0
    return BindingScore(class_key, raw, z, percentile, min(max(scaled, 0.0), 100.0), contrib)

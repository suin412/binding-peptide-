"""문헌에서 보고된 결합 친화도 기록 (Kd, ΔG 등) + 단위 변환 + 본문에서 값 후보 찾기.

예측이 아니라 **논문에 보고된 값**을 서열에 붙여 두고 학습용 예측 점수와 비교하기 위한 모듈이다.
이 모듈은 provenance.py에서 가져다 쓰므로 provenance를 import하지 않는다.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

R_KCAL = 0.0019872041        # kcal/(mol·K)
T_DEFAULT = 298.15           # K

STATUS_CONFIRMED = "확정"
STATUS_SUGGESTED = "제안"

METRIC_KD = "Kd"
METRIC_KA = "Ka"
METRIC_DG = "ΔG"
METRIC_COMPUTED = "계산 결합 에너지"
METRIC_EC50 = "EC50/IC50"
METRIC_ADSORPTION = "흡착량·피복률"
METRIC_OTHER = "기타"
AFFINITY_METRICS = [METRIC_KD, METRIC_KA, METRIC_DG, METRIC_COMPUTED, METRIC_EC50, METRIC_ADSORPTION, METRIC_OTHER]

UNITS = {
    METRIC_KD: ["pM", "nM", "µM", "mM", "M"],
    METRIC_KA: ["M⁻¹", "µM⁻¹", "nM⁻¹"],
    METRIC_DG: ["kcal/mol", "kJ/mol"],
    METRIC_COMPUTED: ["kcal/mol", "kJ/mol"],
    METRIC_EC50: ["pM", "nM", "µM", "mM", "M"],
    METRIC_ADSORPTION: ["ng/cm²", "µg/cm²", "%", "기타"],
    METRIC_OTHER: ["기타"],
}

CALCULATED_METHOD = "계산(분자동역학·PepBD 등)"
AFFINITY_METHODS = [
    "QCM-D", "QCM", "SPR (표면 플라즈몬 공명)", "BLI", "AFM 힘 분광", "엘립소메트리", "형광/발광 흡착 정량",
    "ELISA형 결합 어세이", "ITC", "FP/FRET", "MST", "유세포분석", "접촉각·XPS 등 표면 분석 (정성)",
    CALCULATED_METHOD, "기타·불명",
]

# 본문에서 측정 방법 이름을 찾을 때 쓰는 키워드 (제안용 — 사용자가 확정)
_METHOD_PATTERNS = {
    "QCM-D": r"QCM-?D\b|quartz crystal microbalance with dissipation",
    "QCM": r"\bQCM\b|quartz crystal microbalance",
    "SPR (표면 플라즈몬 공명)": r"\bSPR\b|surface plasmon resonance",
    "BLI": r"\bBLI\b|bio-?layer interferometry",
    "AFM 힘 분광": r"\bAFM\b|atomic force",
    "엘립소메트리": r"ellipsometr",
    "형광/발광 흡착 정량": r"fluorescen\w+ (?:assay|intensity|quantif)|FITC",
    "ELISA형 결합 어세이": r"\bELISA\b",
    "ITC": r"\bITC\b|isothermal titration",
    "FP/FRET": r"fluorescence polari[sz]ation|\bFRET\b",
    "MST": r"microscale thermophoresis|\bMST\b",
    "유세포분석": r"flow cytometr",
    CALCULATED_METHOD: r"molecular dynamics|\bMD simulation|PepBD|free energy|Monte Carlo|docking",
}


def to_molar(value: float, unit: str) -> float | None:
    factors = {"pM": 1e-12, "nM": 1e-9, "µM": 1e-6, "mM": 1e-3, "M": 1.0}
    return value * factors[unit] if unit in factors else None


def _kcal(value: float, unit: str) -> float | None:
    if unit == "kcal/mol":
        return value
    if unit == "kJ/mol":
        return value / 4.184
    return None


def delta_g_kcal(metric: str, value: float, unit: str, temperature: float = T_DEFAULT) -> float | None:
    """결합 자유에너지(kcal/mol)로 환산한 비교용 값. 더 음수일수록 강한 결합. 환산 불가면 None.

    Kd: ΔG = RT ln(Kd[M]) · Ka: Kd = 1/Ka · ΔG/계산 결합 에너지: 단위만 변환.
    EC50/IC50은 겉보기 값이라 같은 식을 쓰지 않는다(None).
    """
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    if metric == METRIC_KD:
        molar = to_molar(value, unit)
        return R_KCAL * temperature * math.log(molar) if molar and molar > 0 else None
    if metric == METRIC_KA:
        factors = {"M⁻¹": 1.0, "µM⁻¹": 1e6, "nM⁻¹": 1e9}
        ka = value * factors.get(unit, math.nan)
        return -R_KCAL * temperature * math.log(ka) if ka and ka > 0 and not math.isnan(ka) else None
    if metric in (METRIC_DG, METRIC_COMPUTED):
        return _kcal(value, unit)
    return None


@dataclass
class Affinity:
    metric: str = METRIC_KD
    value: float | None = None
    unit: str = "nM"
    method: str = ""
    target: str = ""
    conditions: str = ""          # pH, 완충액, 이온 세기, 온도 …
    source: str = ""              # PMID / DOI 등
    status: str = STATUS_CONFIRMED
    note: str = ""

    @property
    def computed(self) -> bool:
        return self.metric == METRIC_COMPUTED or self.method == CALCULATED_METHOD

    @property
    def delta_g(self) -> float | None:
        return delta_g_kcal(self.metric, self.value, self.unit) if self.value is not None else None

    def label(self) -> str:
        if self.value is None:
            return "미보고"
        text = f"{self.metric} {self.value:g} {self.unit}"
        if self.method and not self.computed:        # 계산값은 아래에서 "계산값"으로 표시하므로 방법 이름을 반복하지 않는다
            text += f" ({self.method})"
        if self.computed:
            text += " · 계산값"
        return text


# ---- 본문에서 친화도 후보 찾기 --------------------------------------------------------

@dataclass
class AffinityCandidate:
    metric: str
    value: float
    unit: str
    context: str
    methods: list[str] = field(default_factory=list)      # 근처에서 발견된 측정 방법 이름 (제안)

    def label(self) -> str:
        return f"{self.metric} {self.value:g} {self.unit}"


_NUM = r"(\d+(?:\.\d+)?)"
_SCI = rf"{_NUM}(?:\s*[×x]\s*10\s*\^?\s*([−\-–]?\s*\d+))?"
_PREFIX = r"(?:=|:|of|was|were|is|are|≈|~|about|approximately|around|\s)*"
_KD = re.compile(rf"\bK\s*_?\s*[Dd]\b{_PREFIX}{_SCI}\s*(?:±\s*\d+(?:\.\d+)?\s*)?(pM|nM|[µμu]M|mM|M)\b")
_KA = re.compile(rf"\bK\s*_?\s*[Aa]\b{_PREFIX}{_SCI}\s*(?:±\s*\d+(?:\.\d+)?\s*)?(M\s*[⁻−\-–^]\s*1|[µμu]M\s*[⁻−\-–^]\s*1)")
_DG = re.compile(
    rf"(?:ΔG(?:\s*_?\s*bind\w*)?|dG|binding (?:free )?energy)\s*{_PREFIX}([−\-–]?\s*\d+(?:\.\d+)?)\s*(kcal\s*/\s*mol|kJ\s*/\s*mol)",
    re.IGNORECASE,
)


def _num(mant: str, exp: str | None) -> float:
    value = float(mant)
    if exp:
        value *= 10 ** int(exp.replace("−", "-").replace("–", "-").replace(" ", ""))
    return value


def _unit(u: str) -> str:
    return "µM" if u in ("μM", "uM", "µM") else u


def _methods_near(text: str, start: int, end: int, window: int = 400) -> list[str]:
    snippet = text[max(0, start - window): end + window]
    return [name for name, pat in _METHOD_PATTERNS.items() if re.search(pat, snippet, re.IGNORECASE)]


def extract_affinity_candidates(text: str) -> list[AffinityCandidate]:
    """본문에서 Kd·Ka·ΔG 형태의 값을 후보로 찾는다. 어느 펩타이드의 값인지는 문맥을 보고 사용자가 판단한다."""
    text = text or ""
    found: list[tuple[int, AffinityCandidate]] = []
    for m in _KD.finditer(text):
        found.append((m.start(), AffinityCandidate(METRIC_KD, _num(m.group(1), m.group(2)), _unit(m.group(3)),
                                                   _ctx(text, m.start(), m.end()), _methods_near(text, m.start(), m.end()))))
    for m in _KA.finditer(text):
        unit = "µM⁻¹" if m.group(3).lstrip().startswith(("µ", "μ", "u")) else "M⁻¹"
        found.append((m.start(), AffinityCandidate(METRIC_KA, _num(m.group(1), m.group(2)), unit,
                                                   _ctx(text, m.start(), m.end()), _methods_near(text, m.start(), m.end()))))
    for m in _DG.finditer(text):
        value = float(m.group(1).replace("−", "-").replace("–", "-").replace(" ", ""))
        unit = "kJ/mol" if "kj" in m.group(2).lower() else "kcal/mol"
        found.append((m.start(), AffinityCandidate(METRIC_DG, value, unit, _ctx(text, m.start(), m.end()),
                                                   _methods_near(text, m.start(), m.end()))))
    found.sort(key=lambda x: x[0])
    seen, out = set(), []
    for _, cand in found:
        key = (cand.metric, cand.value, cand.unit, cand.context)
        if key not in seen:
            seen.add(key)
            out.append(cand)
    return out


def _ctx(text: str, start: int, end: int, pad: int = 90) -> str:
    return re.sub(r"\s+", " ", text[max(0, start - pad): end + pad]).strip()


def rank_correlation(xs: list[float], ys: list[float]) -> float | None:
    """스피어만 순위상관 (동률은 평균 순위). 점이 3개 미만이거나 분산이 없으면 None."""
    n = len(xs)
    if n < 3 or n != len(ys):
        return None

    def ranks(values: list[float]) -> list[float]:
        order = sorted(range(n), key=lambda i: values[i])
        out = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and values[order[j + 1]] == values[order[i]]:
                j += 1
            for k in range(i, j + 1):
                out[order[k]] = (i + j) / 2 + 1
            i = j + 1
        return out

    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx == 0 or syy == 0:
        return None
    return sum((a - mx) * (b - my) for a, b in zip(rx, ry)) / math.sqrt(sxx * syy)

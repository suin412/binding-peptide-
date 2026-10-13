"""가상 스크리닝: 필터 깔때기 → 학습용 점수 → 순위화.

점수식은 **학습용 휴리스틱**이며 실제 활성·결합을 예측하지 않는다.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from .align import align_sequences
from .binding_score import score_sequence
from .properties import compute_properties, hydrophobic_moment


@dataclass
class FilterConfig:
    min_len: int = 5
    max_len: int = 100
    min_charge: float = -50.0
    max_charge: float = 50.0
    min_gravy: float = -5.0
    max_gravy: float = 5.0
    max_cys: int = 99
    max_instability: float | None = None      # None이면 단계 생략
    forbidden_motifs: list[str] = field(default_factory=list)


@dataclass
class ScoreConfig:
    target_charge: float = 4.0     # 양이온성 양친매성 펩타이드 가정(학습용)
    w_identity: float = 0.4
    w_moment: float = 0.3
    w_charge: float = 0.3
    w_binding: float = 0.0                 # 표면 결합 성향 가중치 (0이면 사용 안 함, 기존 점수와 동일)
    surface_class: str | None = None       # binding_score.CLASSES의 키


def build_table(library: list[tuple[str, str]]) -> pd.DataFrame:
    """[(변이 표기, 서열)] → 물성 표."""
    rows = []
    for label, seq in library:
        p = compute_properties(seq)
        rows.append({
            "variant": label, "sequence": seq, "length": p["length"], "charge": p["net_charge_pH7.4"],
            "gravy": p["gravy"], "instability": p["instability_index"], "cys": seq.count("C"),
            "pI": p["pI"], "mw": p["mw"],
        })
    return pd.DataFrame(rows)


def run_funnel(df: pd.DataFrame, cfg: FilterConfig) -> tuple[pd.DataFrame, list[tuple[str, int]]]:
    """단계별 필터를 순차 적용. (통과 표, [(단계명, 통과 수)])를 반환."""
    stages = [("라이브러리 전체", len(df))]
    cur = df
    cur = cur[(cur["length"] >= cfg.min_len) & (cur["length"] <= cfg.max_len)]
    stages.append((f"길이 {cfg.min_len}-{cfg.max_len}", len(cur)))
    cur = cur[(cur["charge"] >= cfg.min_charge) & (cur["charge"] <= cfg.max_charge)]
    stages.append((f"순전하 {cfg.min_charge:g}~{cfg.max_charge:g}", len(cur)))
    cur = cur[(cur["gravy"] >= cfg.min_gravy) & (cur["gravy"] <= cfg.max_gravy)]
    stages.append((f"GRAVY {cfg.min_gravy:g}~{cfg.max_gravy:g}", len(cur)))
    cur = cur[cur["cys"] <= cfg.max_cys]
    stages.append((f"Cys ≤ {cfg.max_cys}", len(cur)))
    motifs = [m.strip().upper() for m in cfg.forbidden_motifs if m.strip()]
    if motifs:
        mask = cur["sequence"].apply(lambda s: not any(m in s for m in motifs))
        cur = cur[mask]
        stages.append((f"금지 모티프 없음 ({', '.join(motifs)})", len(cur)))
    if cfg.max_instability is not None:
        cur = cur[cur["instability"] <= cfg.max_instability]
        stages.append((f"불안정성 지수 ≤ {cfg.max_instability:g}", len(cur)))
    return cur.reset_index(drop=True), stages


def _identity_to_parent(parent: str, seq: str) -> float:
    if len(parent) == len(seq):
        return sum(a == b for a, b in zip(parent, seq)) / len(seq)
    return align_sequences(parent, seq).identity


def score_table(df: pd.DataFrame, parent: str, cfg: ScoreConfig) -> pd.DataFrame:
    """부모 동일성 + 소수성 모멘트 + 전하 근접도를 합산한 0~100 학습용 점수로 정렬."""
    out = df.copy()
    if out.empty:
        out["identity"] = out["moment"] = out["score"] = []
        return out
    out["identity"] = [round(_identity_to_parent(parent, s), 3) for s in out["sequence"]]
    out["moment"] = [round(hydrophobic_moment(s), 3) for s in out["sequence"]]
    target = cfg.target_charge if cfg.target_charge else 1.0
    charge_term = (1 - (out["charge"] - cfg.target_charge).abs() / abs(target)).clip(lower=0, upper=1)
    moment_term = out["moment"].clip(upper=1.0)
    use_binding = cfg.w_binding > 0 and cfg.surface_class
    if use_binding:
        out["binding"] = [round(score_sequence(s, cfg.surface_class).scaled, 1) for s in out["sequence"]]
    total_w = (cfg.w_identity + cfg.w_moment + cfg.w_charge + (cfg.w_binding if use_binding else 0.0)) or 1.0
    combined = cfg.w_identity * out["identity"] + cfg.w_moment * moment_term + cfg.w_charge * charge_term
    if use_binding:
        combined = combined + cfg.w_binding * out["binding"] / 100.0
    out["score"] = (100 * combined / total_w).round(1)
    return out.sort_values("score", ascending=False).reset_index(drop=True)

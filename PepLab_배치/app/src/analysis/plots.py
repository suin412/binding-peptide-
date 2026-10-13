"""matplotlib 시각화: 소수성 프로파일, 아미노산 조성, helical wheel, 스크리닝 깔때기."""
from __future__ import annotations

import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from Bio.SeqUtils.ProtParamData import kd  # noqa: E402

from matplotlib import font_manager  # noqa: E402

from .properties import composition, hydropathy_profile  # noqa: E402

_KOREAN_FONTS = ["Malgun Gothic", "AppleGothic", "NanumGothic", "Noto Sans CJK KR", "Noto Sans KR"]


def _setup_korean_font() -> None:
    """설치된 한글 폰트가 있으면 사용한다 (없으면 기본 폰트 — 한글이 네모로 보일 수 있음)."""
    installed = {f.name for f in font_manager.fontManager.ttflist}
    for name in _KOREAN_FONTS:
        if name in installed:
            matplotlib.rcParams["font.family"] = name
            break
    matplotlib.rcParams["axes.unicode_minus"] = False


_setup_korean_font()

_HYDROPHOBIC = set("AVLIMFWC")
_POSITIVE = set("KRH")
_NEGATIVE = set("DE")


def hydropathy_plot(series: dict[str, str], window: int = 7):
    """{이름: 서열} → Kyte-Doolittle 프로파일 (여러 서열 겹쳐 그리기)."""
    fig, ax = plt.subplots(figsize=(7, 3))
    for name, seq in series.items():
        ax.plot(range(1, len(seq) + 1), hydropathy_profile(seq, window), label=name, linewidth=1.8)
    ax.axhline(0, color="gray", linewidth=0.8)
    ax.set_xlabel("잔기 위치")
    ax.set_ylabel(f"소수성 (KD, 윈도우 {window})")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def composition_plot(series: dict[str, str]):
    """아미노산 조성(%) 막대그래프. 여러 서열은 나란히 그린다."""
    import numpy as np

    aas = list("ACDEFGHIKLMNPQRSTVWY")
    fig, ax = plt.subplots(figsize=(7, 3))
    n = len(series)
    width = 0.8 / max(n, 1)
    x = np.arange(len(aas))
    for i, (name, seq) in enumerate(series.items()):
        comp = composition(seq)
        ax.bar(x + i * width, [comp[a] for a in aas], width=width, label=name)
    ax.set_xticks(x + width * (n - 1) / 2)
    ax.set_xticklabels(aas)
    ax.set_ylabel("조성 (%)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def helical_wheel(seq: str, title: str = ""):
    """Helical wheel (100°/잔기). 소수성=주황, 양전하=파랑, 음전하=빨강, 그 외=회색."""
    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(plt.Circle((0, 0), 1.0, fill=False, color="lightgray"))
    pts = []
    for i, aa in enumerate(seq):
        ang = math.radians(100 * i - 90)
        r = 1.0
        x, y = r * math.cos(ang), -r * math.sin(ang)
        pts.append((x, y))
        if aa in _POSITIVE:
            color = "#4c78a8"
        elif aa in _NEGATIVE:
            color = "#e45756"
        elif aa in _HYDROPHOBIC:
            color = "#f58518"
        else:
            color = "#bab0ac"
        ax.add_patch(plt.Circle((x, y), 0.11, color=color, zorder=3))
        ax.text(x, y, aa, ha="center", va="center", fontsize=8, zorder=4)
        ax.text(x * 1.22, y * 1.22, str(i + 1), ha="center", va="center", fontsize=6, color="gray")
    for (x1, y1), (x2, y2) in zip(pts, pts[1:]):
        ax.plot([x1, x2], [y1, y2], color="#dddddd", linewidth=0.8, zorder=1)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_title(title or "Helical wheel", fontsize=10)
    fig.tight_layout()
    return fig


def funnel_plot(stages: list[tuple[str, int]]):
    """스크리닝 깔때기: 단계별 통과 수 가로 막대."""
    fig, ax = plt.subplots(figsize=(7, 0.55 * max(len(stages), 2) + 0.8))
    names = [s[0] for s in stages][::-1]
    counts = [s[1] for s in stages][::-1]
    bars = ax.barh(names, counts, color="#4c78a8")
    top = max(counts) if counts and max(counts) else 1
    for bar, c in zip(bars, counts):
        ax.text(bar.get_width() + top * 0.01, bar.get_y() + bar.get_height() / 2, str(c), va="center", fontsize=9)
    ax.set_xlim(0, top * 1.15)
    ax.set_xlabel("통과한 후보 수")
    fig.tight_layout()
    return fig


def residue_contribution_plot(seq: str, contributions: list[float], title: str = ""):
    """잔기별 결합 성향 가중치 막대 (양수=유리, 음수=불리)."""
    fig, ax = plt.subplots(figsize=(max(5.0, 0.32 * len(seq) + 1.5), 2.8))
    colors = ["#4c78a8" if c > 0 else "#e45756" if c < 0 else "#bab0ac" for c in contributions]
    ax.bar(range(len(seq)), contributions, color=colors)
    ax.set_xticks(range(len(seq)))
    ax.set_xticklabels(list(seq), fontsize=8)
    ax.axhline(0, color="gray", linewidth=0.8)
    ax.set_ylabel("잔기 가중치")
    ax.set_title(title or "잔기별 기여", fontsize=10)
    fig.tight_layout()
    return fig


def affinity_scatter_plot(points: list[dict], rho: float | None = None):
    """예측 점수(x) vs 보고된 결합 세기(y, -ΔG). points: {name, x, dg, computed}. 위쪽일수록 강한 결합."""
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for computed, color, label in ((False, "#4c78a8", "측정값"), (True, "#f58518", "계산값")):
        pts = [p for p in points if p["computed"] == computed]
        if pts:
            ax.scatter([p["x"] for p in pts], [p["dg"] for p in pts], color=color, label=label, s=45, zorder=3)
    for p in points:
        ax.annotate(p["name"][:16], (p["x"], p["dg"]), fontsize=7, xytext=(4, 4), textcoords="offset points")
    ax.set_xlim(-2, 102)
    ax.invert_yaxis()                       # ΔG가 더 음수(강한 결합)일수록 위쪽
    ax.set_xlabel("결합 성향 점수 (가중치 척도 0–100)")
    ax.set_ylabel("ΔG (kcal/mol) — 위쪽일수록 강한 결합")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    if rho is not None:
        ax.set_title(f"순위상관 ρ = {rho:.2f} (n={len(points)})", fontsize=10)
    fig.tight_layout()
    return fig
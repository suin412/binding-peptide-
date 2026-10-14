"""PDB 텍스트 유틸: pLDDT(B-factor) 추출·정규화, Cα 개수, 신뢰도 구간."""
from __future__ import annotations

PLDDT_BANDS = [
    (90, "매우 높음", "#0053d6"),
    (70, "높음", "#65cbf3"),
    (50, "낮음", "#ffdb13"),
    (0, "매우 낮음", "#ff7d45"),
]


def ca_records(pdb_text: str) -> list[tuple[int, float]]:
    """[(잔기 번호, B-factor)] — Cα 원자 기준."""
    out = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                out.append((int(line[22:26]), float(line[60:66])))
            except ValueError:
                continue
    return out


def normalize_plddt(pdb_text: str) -> str:
    """B-factor가 0~1 범위(ESMFold)면 0~100으로 환산해 B-factor 열을 다시 쓴다. 이미 0~100이면 그대로."""
    recs = ca_records(pdb_text)
    if not recs or max(b for _, b in recs) > 1.0:
        return pdb_text
    lines = []
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
            try:
                value = float(line[60:66]) * 100
                line = f"{line[:60]}{value:6.2f}{line[66:]}"
            except ValueError:
                pass
        lines.append(line)
    return "\n".join(lines) + "\n"


def plddt_per_residue(pdb_text: str) -> list[tuple[int, float]]:
    return ca_records(normalize_plddt(pdb_text))


def mean_plddt(pdb_text: str) -> float:
    recs = plddt_per_residue(pdb_text)
    return sum(b for _, b in recs) / len(recs) if recs else 0.0


def band(value: float) -> tuple[str, str]:
    """pLDDT 값 → (구간 이름, 색)."""
    for threshold, label, color in PLDDT_BANDS:
        if value >= threshold:
            return label, color
    return PLDDT_BANDS[-1][1], PLDDT_BANDS[-1][2]


def residue_count(pdb_text: str) -> int:
    return len(ca_records(pdb_text))

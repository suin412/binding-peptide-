"""변이 파서·적용, 라이브러리 생성, 융합, 보존 치환 제안.

모든 위치는 1-based이며 '원본 서열 기준'이다 (여러 변이를 한 번에 적용해도 위치가 밀리지 않는다).
"""
from __future__ import annotations

import random
import re

from ..sequence import STANDARD_AA

MAX_LIBRARY = 2000

CONSERVATIVE_GROUPS = [
    "AVLIM",   # 지방족 소수성
    "FWY",     # 방향족
    "KRH",     # 염기성
    "DE",      # 산성
    "STNQ",    # 극성
    "GP",      # 특수
]


class MutationError(ValueError):
    """변이 입력 오류. 메시지는 사용자에게 그대로 보여준다."""


# ---- 파서 -------------------------------------------------------------

_POINT = re.compile(r"^([A-Za-z])(\d+)([A-Za-z])$")
_DEL = re.compile(r"^del\s*(\d+)(?:\s*-\s*(\d+))?$", re.IGNORECASE)
_INS = re.compile(r"^ins\s*(\d+)\s+([A-Za-z]+)$", re.IGNORECASE)
_TRUNC = re.compile(r"^trunc\s*(\d+)\s*-\s*(\d+)$", re.IGNORECASE)


def parse_mutations(text: str) -> list[tuple]:
    """'K12A, del 5-9, ins 10 GGS, trunc 1-20' → 연산 목록. 구분자는 쉼표/세미콜론/줄바꿈."""
    ops = []
    for raw in re.split(r"[,;\n]+", text or ""):
        item = raw.strip()
        if not item:
            continue
        if m := _POINT.match(item):
            ops.append(("point", m.group(1).upper(), int(m.group(2)), m.group(3).upper()))
        elif m := _DEL.match(item):
            start = int(m.group(1))
            end = int(m.group(2)) if m.group(2) else start
            ops.append(("del", start, end))
        elif m := _INS.match(item):
            ops.append(("ins", int(m.group(1)), m.group(2).upper()))
        elif m := _TRUNC.match(item):
            ops.append(("trunc", int(m.group(1)), int(m.group(2))))
        else:
            raise MutationError(
                f"'{item}'을(를) 해석할 수 없습니다. 예: K12A / del 5-9 / ins 10 GGS / trunc 1-20"
            )
    if not ops:
        raise MutationError("변이를 입력해 주세요. 예: K12A, L15P")
    return ops


def _check_range(seq: str, start: int, end: int) -> None:
    if start < 1 or end > len(seq) or start > end:
        raise MutationError(f"위치 {start}-{end}이(가) 서열 범위(1-{len(seq)})를 벗어났습니다.")


def apply_mutations(seq: str, text: str) -> tuple[str, str]:
    """변이 문자열을 적용해 (새 서열, 정규화된 변이 표기)를 반환한다."""
    ops = parse_mutations(text)
    items: list[list] = [[i + 1, aa] for i, aa in enumerate(seq)]  # [원본위치|None, 잔기]
    labels = []
    for op in ops:
        kind = op[0]
        if kind == "point":
            _, orig, pos, new = op
            if pos < 1 or pos > len(seq):
                raise MutationError(f"위치 {pos}이(가) 서열 범위(1-{len(seq)})를 벗어났습니다.")
            if seq[pos - 1] != orig:
                raise MutationError(f"{pos}번째 잔기는 {seq[pos - 1]}인데 {orig}로 입력했습니다.")
            if new not in STANDARD_AA:
                raise MutationError(f"'{new}'는 표준 아미노산이 아닙니다.")
            target = next((it for it in items if it[0] == pos), None)
            if target is None:
                raise MutationError(f"{pos}번째 잔기는 이미 삭제되어 변이를 적용할 수 없습니다.")
            target[1] = new
            labels.append(f"{orig}{pos}{new}")
        elif kind == "del":
            _, start, end = op
            _check_range(seq, start, end)
            items = [it for it in items if it[0] is None or not (start <= it[0] <= end)]
            labels.append(f"del{start}-{end}" if end != start else f"del{start}")
        elif kind == "ins":
            _, pos, ins = op
            if pos < 0 or pos > len(seq):
                raise MutationError(f"삽입 위치 {pos}이(가) 범위(0-{len(seq)})를 벗어났습니다 (0=맨 앞).")
            bad = [c for c in ins if c not in STANDARD_AA]
            if bad:
                raise MutationError(f"삽입 서열에 표준 아미노산이 아닌 문자가 있습니다: {''.join(bad)}")
            idx = 0 if pos == 0 else next((i + 1 for i, it in enumerate(items) if it[0] == pos), None)
            if idx is None:
                raise MutationError(f"{pos}번째 잔기는 이미 삭제되어 그 뒤에 삽입할 수 없습니다.")
            items[idx:idx] = [[None, c] for c in ins]
            labels.append(f"ins{pos}{ins}")
        elif kind == "trunc":
            _, start, end = op
            _check_range(seq, start, end)
            items = [it for it in items if it[0] is None or start <= it[0] <= end]
            labels.append(f"trunc{start}-{end}")
    result = "".join(it[1] for it in items)
    if not result:
        raise MutationError("변이를 적용하면 서열이 비어 버립니다.")
    return result, "+".join(labels)


# ---- 라이브러리 생성 -----------------------------------------------------

def _cap(items: list, limit: int) -> list:
    return items[:limit]


def alanine_scan(seq: str, limit: int = MAX_LIBRARY) -> list[tuple[str, str]]:
    """각 위치를 Ala로 치환 (이미 Ala인 위치는 건너뜀)."""
    out = [(f"{aa}{i}A", seq[:i - 1] + "A" + seq[i:]) for i, aa in enumerate(seq, start=1) if aa != "A"]
    return _cap(out, limit)


def site_saturation(seq: str, positions: list[int] | None = None, limit: int = MAX_LIBRARY) -> list[tuple[str, str]]:
    """지정 위치(기본: 전체)를 나머지 19개 아미노산으로 치환."""
    positions = positions or list(range(1, len(seq) + 1))
    out = []
    for pos in positions:
        if pos < 1 or pos > len(seq):
            raise MutationError(f"위치 {pos}이(가) 서열 범위(1-{len(seq)})를 벗어났습니다.")
        orig = seq[pos - 1]
        for aa in STANDARD_AA:
            if aa != orig:
                out.append((f"{orig}{pos}{aa}", seq[:pos - 1] + aa + seq[pos:]))
    return _cap(out, limit)


def random_variants(seq: str, n: int, max_mutations: int = 3, seed: int | None = None,
                    limit: int = MAX_LIBRARY) -> list[tuple[str, str]]:
    """무작위로 1~max_mutations개 위치를 치환한 변이체 n개 (중복 제외)."""
    rng = random.Random(seed)
    n = min(n, limit)
    seen, out = {seq}, []
    attempts = 0
    while len(out) < n and attempts < n * 50:
        attempts += 1
        k = rng.randint(1, min(max_mutations, len(seq)))
        positions = sorted(rng.sample(range(1, len(seq) + 1), k))
        chars = list(seq)
        labels = []
        for pos in positions:
            new = rng.choice([a for a in STANDARD_AA if a != seq[pos - 1]])
            labels.append(f"{seq[pos - 1]}{pos}{new}")
            chars[pos - 1] = new
        variant = "".join(chars)
        if variant not in seen:
            seen.add(variant)
            out.append(("+".join(labels), variant))
    return out


def truncations(seq: str, min_len: int = 5, limit: int = MAX_LIBRARY) -> list[tuple[str, str]]:
    """N-말단/C-말단 절단체 세트."""
    out = []
    for i in range(2, len(seq) - min_len + 2):
        out.append((f"ΔN{i - 1}", seq[i - 1:]))
    for j in range(len(seq) - 1, min_len - 1, -1):
        out.append((f"ΔC{len(seq) - j}", seq[:j]))
    return _cap(out, limit)


def reverse_sequence(seq: str) -> str:
    return seq[::-1]


def shuffle_sequence(seq: str, seed: int | None = None) -> str:
    chars = list(seq)
    random.Random(seed).shuffle(chars)
    return "".join(chars)


def fuse(seq_a: str, seq_b: str, linker: str = "") -> str:
    bad = [c for c in linker.upper() if c not in STANDARD_AA]
    if bad:
        raise MutationError(f"링커에 표준 아미노산이 아닌 문자가 있습니다: {''.join(bad)}")
    return seq_a + linker.upper() + seq_b


def conservative_options(aa: str) -> list[str]:
    """같은 성질 그룹의 다른 아미노산 목록 (보존 치환 후보)."""
    for group in CONSERVATIVE_GROUPS:
        if aa in group:
            return [c for c in group if c != aa]
    return []


def conservative_suggestions(seq: str, pos: int) -> list[tuple[str, str]]:
    if pos < 1 or pos > len(seq):
        raise MutationError(f"위치 {pos}이(가) 서열 범위(1-{len(seq)})를 벗어났습니다.")
    orig = seq[pos - 1]
    return [(f"{orig}{pos}{new}", seq[:pos - 1] + new + seq[pos:]) for new in conservative_options(orig)]

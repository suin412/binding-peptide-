"""서열 → 논문: 이 서열이 언급된 논문을 찾고, 어디(초록/본문)에 어떤 문맥으로 나오는지 확인한다.

Europe PMC는 본문 색인에서 서열 문자열을 그대로 검색할 수 있다(예: 'WEWWHVFHHRLR' → 해당 논문 1편).
검색은 단어(토큰) 단위 일치라서 정확히 같은 서열만 찾고, 더 긴 서열에 들어 있는 경우는 와일드카드(*서열*)로 따로 찾는다.
검색된 논문마다 초록과 (PMC에 공개된) 본문을 읽어 실제 위치와 문맥을 확인한다. 못 찾으면 "위치 미확인"으로 표시한다.
"""
from __future__ import annotations

import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import Callable

from ..sequence import THREE_TO_ONE, SequenceError, validate_sequence
from . import europepmc
from .europepmc import Paper
from .http import ApiError
from .seqextract import extract_candidates

MIN_LEN = 5
SHORT_WARN_LEN = 7
FULLTEXT_WORKERS = 3

KIND_EXACT = "정확히 일치"
KIND_CONTAINED = "더 긴 서열에 포함"
KIND_THREE = "3문자 표기로 일치"
WHERE_ABSTRACT = "초록"
WHERE_FULLTEXT = "전문"
WHERE_UNKNOWN = "위치 미확인"


@dataclass
class Occurrence:
    kind: str
    where: str            # 초록 | 전문
    context: str


@dataclass
class SeqMatch:
    paper: Paper
    occurrences: list[Occurrence] = field(default_factory=list)

    @property
    def located(self) -> bool:
        return bool(self.occurrences)

    @property
    def best_kind(self) -> str:
        order = [KIND_EXACT, KIND_THREE, KIND_CONTAINED]
        kinds = {o.kind for o in self.occurrences}
        return next((k for k in order if k in kinds), WHERE_UNKNOWN)

    @property
    def where(self) -> str:
        places = sorted({o.where for o in self.occurrences}, key=[WHERE_ABSTRACT, WHERE_FULLTEXT].index)
        return "/".join(places) if places else WHERE_UNKNOWN


@dataclass
class SeqSearchResult:
    sequence: str
    queries: list[str]
    hit_count: int
    matches: list[SeqMatch]
    fulltext_read: int = 0
    fulltext_failed: int = 0
    warnings: list[str] = field(default_factory=list)


def three_letter_phrase(seq: str) -> str:
    one_to_three = {v: k.title() for k, v in THREE_TO_ONE.items()}
    return " ".join(one_to_three[a] for a in seq)


def build_query(seq: str, *, contained: bool = False, three_letter: bool = False) -> str:
    """Europe PMC 검색식. 기본은 정확한 서열 토큰, 옵션으로 *서열*(더 긴 서열 포함)과 3문자 표기 구문을 OR로 더한다."""
    parts = [f'"{seq}"']
    if contained:
        parts.append(f"*{seq}*")
    if three_letter:
        parts.append(f'"{three_letter_phrase(seq)}"')
    return "(" + " OR ".join(parts) + ")" if len(parts) > 1 else parts[0]


def _context(text: str, start: int, end: int, pad: int = 110) -> str:
    return re.sub(r"\s+", " ", text[max(0, start - pad): end + pad]).strip()


def locate(seq: str, text: str, where: str, *, max_hits: int = 5) -> list[Occurrence]:
    """텍스트에서 서열이 나오는 위치를 찾는다 (1문자 정확 일치, 더 긴 서열에 포함, 3문자 표기)."""
    out: list[Occurrence] = []
    seen: set[str] = set()

    def add(kind: str, start: int, end: int) -> None:
        ctx = _context(text, start, end)
        if ctx not in seen and len(out) < max_hits:
            seen.add(ctx)
            out.append(Occurrence(kind, where, ctx))

    for m in re.finditer(rf"(?<![A-Za-z0-9]){re.escape(seq)}(?![A-Za-z0-9])", text):
        add(KIND_EXACT, m.start(), m.end())
    for m in re.finditer(rf"(?<![A-Za-z0-9])[A-Z]*{re.escape(seq)}[A-Z]*(?![A-Za-z0-9])", text):
        if m.group(0) != seq:
            add(KIND_CONTAINED, m.start(), m.end())
    if len(out) < max_hits:
        for cand in extract_candidates(text, min_len=len(seq)):
            if cand.kind != "1문자" and seq in cand.sequence:
                kind = KIND_THREE
                ctx = re.sub(r"\s+", " ", cand.context)
                if ctx not in seen and len(out) < max_hits:
                    seen.add(ctx)
                    out.append(Occurrence(kind, where, ctx))
    return out


def find_papers_for_sequence(
    sequence: str, *, max_papers: int = 50, contained: bool = False, three_letter: bool = False,
    verify_fulltext: bool = True, max_fulltext: int = 40, progress: Callable[[float, str], None] | None = None,
    search_fn=None, fulltext_fn=None,
) -> SeqSearchResult:
    """서열이 언급된 논문을 찾는다. search_fn/fulltext_fn은 테스트용 주입 지점."""
    seq = validate_sequence(sequence, min_len=MIN_LEN)
    search_fn = search_fn or europepmc.search_many
    fulltext_fn = fulltext_fn or europepmc.fetch_fulltext

    def report(frac: float, msg: str) -> None:
        if progress:
            progress(min(max(frac, 0.0), 1.0), msg)

    query = build_query(seq, contained=contained, three_letter=three_letter)
    report(0.0, "Europe PMC 본문 색인에서 서열을 검색하는 중…")
    found = search_fn(query, max_papers=max_papers, sort=europepmc.SORT_RELEVANCE)
    papers = found.papers
    warnings: list[str] = []
    if len(seq) < SHORT_WARN_LEN:
        warnings.append(f"{len(seq)} aa의 짧은 서열은 다른 논문의 우연한 일치가 많을 수 있어 문맥을 꼭 확인하세요.")
    if found.hit_count > len(papers):
        warnings.append(f"전체 일치 {found.hit_count:,}편 중 상위 {len(papers)}편만 확인했습니다.")

    matches = [SeqMatch(p, locate(seq, f"{p.title}. {p.abstract}", WHERE_ABSTRACT)) for p in papers]
    by_key = {m.paper.key: m for m in matches}

    read = failed = 0
    if verify_fulltext:
        todo = [m.paper for m in matches if m.paper.pmcid and not m.located][:max_fulltext]
        done = 0
        if todo:
            with ThreadPoolExecutor(max_workers=FULLTEXT_WORKERS) as pool:
                futures = {pool.submit(fulltext_fn, p.pmcid): p for p in todo}
                results: dict[str, str | None] = {}
                for fut in as_completed(futures):
                    p = futures[fut]
                    done += 1
                    try:
                        results[p.key] = fut.result()
                    except ApiError:
                        results[p.key] = None
                        failed += 1
                    report(0.2 + 0.75 * done / len(todo), f"본문 {done}/{len(todo)} 확인 중…")
            for p in todo:                              # 결과 순서를 고정
                text = results.get(p.key)
                if text:
                    read += 1
                    by_key[p.key].occurrences.extend(locate(seq, text, WHERE_FULLTEXT))

    # 위치가 확인되고 정확히 일치하는 논문을 위로, 위치 미확인은 아래로 (같으면 검색 순서 유지)
    rank = {KIND_EXACT: 0, KIND_THREE: 1, KIND_CONTAINED: 2, WHERE_UNKNOWN: 3}
    matches.sort(key=lambda m: rank[m.best_kind])
    report(1.0, "완료")
    return SeqSearchResult(seq, [query], found.hit_count, matches, read, failed, warnings)


def validate(sequence: str) -> str:
    """UI용: 정리된 서열을 돌려주거나 SequenceError를 던진다 (최소 길이 포함)."""
    return validate_sequence(sequence, min_len=MIN_LEN)


__all__ = ["find_papers_for_sequence", "validate", "SequenceError"]

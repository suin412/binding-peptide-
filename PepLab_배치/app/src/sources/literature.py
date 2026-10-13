"""키워드 하나로 여러 논문의 초록(+선택: 오픈액세스 전문)에서 서열 후보를 바로 모은다.

논문 하나하나를 고르지 않아도 되도록, 모든 논문의 후보를 서열 단위로 묶어 보여주기 위한 모듈.
저장하는 것은 서열·출처 식별자·짧은 문맥뿐이며 논문 본문은 보관하지 않는다.
"""
from __future__ import annotations

from collections import Counter
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import Callable

from . import europepmc
from .binding import TargetSuggestion, material_matches_terms, suggest_targets
from .europepmc import Paper
from .http import ApiError
from .query_assist import expand_query
from .seqextract import Candidate, extract_candidates

SRC_ABSTRACT = "초록"
SRC_FULLTEXT = "전문"
DEFAULT_MIN_LEN = 6


@dataclass
class SeqHit:
    candidate: Candidate
    paper: Paper
    source: str                              # 초록 | 전문
    targets: list[TargetSuggestion] = field(default_factory=list)   # 그 논문 텍스트에서 제안된 결합 대상 후보
    score: float = 0.0                                              # 찾는 키워드와의 관련도 (논문 점수 + 서열 바로 옆에 키워드가 나오는지)

    @property
    def sequence(self) -> str:
        return self.candidate.sequence


@dataclass
class SequenceGroup:
    sequence: str
    hits: list[SeqHit]

    @property
    def paper_count(self) -> int:
        return len({h.paper.key for h in self.hits})

    @property
    def score(self) -> float:
        return max(h.score for h in self.hits)

    @property
    def years(self) -> list[int]:
        return sorted({int(h.paper.year) for h in self.hits if h.paper.year.isdigit()})

    @property
    def has_warning(self) -> bool:
        return any(h.candidate.warnings for h in self.hits)

    @property
    def uncertain(self) -> bool:
        """모든 출현이 '서열인지 불확실'일 때만 불확실로 본다 (한 번이라도 서열 문맥이면 신뢰)."""
        return all(h.candidate.uncertain for h in self.hits)

    @property
    def best(self) -> SeqHit:
        """대표 출현: 서열 문맥이 확실한 것을 우선, 같으면 먼저 발견된 것."""
        return next((h for h in self.hits if not h.candidate.uncertain), self.hits[0])


@dataclass
class LiteratureResult:
    papers: list[Paper]
    hit_count: int
    hits: list[SeqHit]
    groups: list[SequenceGroup]
    fulltext_checked: int = 0
    fulltext_failed: int = 0
    year_counts: dict[int, int] = field(default_factory=dict)       # 가져온 논문의 연도 분포
    fulltext_read: int = 0                                          # 본문을 실제로 읽은 논문 수
    unread: list[Paper] = field(default_factory=list)               # 본문을 읽지 못한(초록만 확인한) 논문
    effective_query: str = ""                                       # Europe PMC에 실제로 보낸 검색식
    expansion_notes: list[str] = field(default_factory=list)        # 약어 확장 등 검색어 보정 안내


def group_hits(hits: list[SeqHit]) -> list[SequenceGroup]:
    """같은 서열을 한 줄로 묶는다. 문맥이 확실한 서열 → 키워드와 관련도가 높은 서열 → 여러 논문에 나온 서열 → 먼저 발견된 순."""
    order: dict[str, int] = {}
    by_seq: dict[str, list[SeqHit]] = {}
    for i, hit in enumerate(hits):
        order.setdefault(hit.sequence, i)
        by_seq.setdefault(hit.sequence, []).append(hit)
    groups = [SequenceGroup(seq, items) for seq, items in by_seq.items()]
    groups.sort(key=lambda g: (g.uncertain, -g.score, -g.paper_count, order[g.sequence]))
    return groups


_PEPTIDE_WORD = re.compile(r"peptid", re.IGNORECASE)
_PLASTIC_WORDS = re.compile(r"plastic|polymer|polystyrene|polyethylene|polypropylene|terephthalate", re.IGNORECASE)
FULLTEXT_WORKERS = 3      # NCBI는 초당 3회 요청 제한이 있어 많이 늘리지 않는다


def _targets(text: str, local: str = "", seq: str = "") -> list[TargetSuggestion]:
    """논문 전체 텍스트에서 제안한 결합 대상 후보. local(서열 바로 주변 문맥)에서 나온 후보를 앞에 둔다."""
    plastic = bool(_PLASTIC_WORDS.search(text))
    paper_level = suggest_targets(text, plastic_context=plastic)
    if not local:
        return paper_level
    near = suggest_targets(local, plastic_context=plastic, anchor=local.find(seq) if seq in local else None)
    names = {t.name for t in near}
    return near + [t for t in paper_level if t.name not in names]


_BIND_WORD = re.compile(r"bind|adsor|affinity|attach|anchor|immobili", re.IGNORECASE)
_SCREEN_WORD = re.compile(r"phage display|biopanning|screen|in silico|computational design|selected", re.IGNORECASE)


def relevance_score(paper: Paper, terms: list[str]) -> float:
    """제목·초록만 보고 "이 키워드에 결합하는 펩타이드를 다룬 논문일 가능성"을 점수화한다 (본문 읽기 우선순위용).

    peptide·binding 언급 + 사용자가 찾는 재료/대상 용어(PET이면 polyethylene terephthalate 등)가 제목이나 초록에 나오는지 본다.
    """
    title, abstract = paper.title.lower(), paper.abstract.lower()
    score = 0.0
    score += 4 if _PEPTIDE_WORD.search(title) else 2 if _PEPTIDE_WORD.search(abstract) else 0
    score += 3 if _BIND_WORD.search(title) else 1.5 if _BIND_WORD.search(abstract) else 0
    for term in terms:
        pattern = re.compile(rf"(?<![a-z0-9]){re.escape(term)}(?![a-z0-9])")
        if pattern.search(title):
            score += 4
        elif pattern.search(abstract):
            score += 2 + min(len(pattern.findall(abstract)) - 1, 2) * 0.5
    score += 1 if _SCREEN_WORD.search(abstract) else 0
    return score


def _term_pattern(term: str) -> re.Pattern:
    return re.compile(rf"(?<![a-z0-9]){re.escape(term)}(?![a-z0-9])")


def hit_score(paper: Paper, terms: list[str], context: str, seq: str = "") -> float:
    """논문 관련도 + 서열 바로 옆 문맥 가산점.

    서열 근처에 소재가 언급돼 있으면 **가장 가까운 소재가 찾는 소재일 때만** 가산한다 (표에 여러 플라스틱의 서열이 나란히 있을 때
    다른 플라스틱의 서열이 같은 점수를 받지 않도록). 소재 언급이 없으면 찾는 용어가 문맥에 있는지만 본다.
    """
    anchor = context.find(seq) if seq and seq in context else None
    nearby = suggest_targets(context, plastic_context=True, anchor=anchor)
    if nearby:
        bonus = 5.0 if material_matches_terms(nearby[0].name, terms) else 0.0
    else:
        near = context.lower()
        bonus = 5.0 if any(_term_pattern(t).search(near) for t in terms) else 0.0
    return relevance_score(paper, terms) + bonus


def _fulltext_order(papers: list[Paper], terms: list[str] | None = None) -> list[Paper]:
    """본문을 읽을 논문: PMCID가 있는 것 중 관련도 점수가 높은 순(같으면 검색 순위 순)."""
    with_pmc = [p for p in papers if p.pmcid]
    return sorted(with_pmc, key=lambda p: -relevance_score(p, terms or []))


def find_sequences(
    query: str, *, max_papers: int = 100, year_from: int | None = None, year_to: int | None = None,
    sort: str = europepmc.SORT_BALANCED, use_fulltext: bool = False, max_fulltext: int = 20,
    min_len: int = DEFAULT_MIN_LEN, progress: Callable[[float, str], None] | None = None,
    expand_aliases: bool = True, search_fn=None, fulltext_fn=None,
) -> LiteratureResult:
    """키워드로 논문을 모아 서열 후보를 추출한다. search_fn/fulltext_fn은 테스트용 주입 지점.

    1단계: 모든 논문의 제목·초록에서 추출. 2단계(선택): PMCID가 있는 논문의 전문에서 추가 추출
    (오픈액세스 표시가 없는 저자 원고도 PMC에 본문이 있을 수 있어 시도한다. 실패는 건너뛴다).
    """
    search_fn = search_fn or europepmc.search_many
    fulltext_fn = fulltext_fn or europepmc.fetch_fulltext

    def report(frac: float, msg: str) -> None:
        if progress:
            progress(min(max(frac, 0.0), 1.0), msg)

    parsed = expand_query(query)
    expanded = parsed if expand_aliases else None
    effective = expanded.text if expanded else query
    report(0.0, "Europe PMC에서 논문을 검색하는 중…")
    found = search_fn(effective, max_papers=max_papers, year_from=year_from, year_to=year_to, sort=sort)
    papers = found.papers
    hits: list[SeqHit] = []
    abstract_seqs: dict[str, set[str]] = {}
    total = max(len(papers), 1)
    span = 0.4 if use_fulltext else 0.9           # 초록 단계가 진행률에서 차지하는 비율

    for i, paper in enumerate(papers):
        report(0.05 + span * i / total, f"초록 {i + 1}/{len(papers)} 분석 중…")
        text = f"{paper.title}. {paper.abstract}".strip()
        seen = abstract_seqs.setdefault(paper.key, set())
        for cand in extract_candidates(text, min_len=min_len):
            seen.add(cand.sequence)
            hits.append(SeqHit(cand, paper, SRC_ABSTRACT, _targets(text, cand.context, cand.sequence),
                               hit_score(paper, parsed.terms, cand.context, cand.sequence)))

    ft_checked = ft_failed = 0
    results: dict[str, str | None] = {}
    if use_fulltext:
        targets_for_ft = _fulltext_order(papers, parsed.terms)[:max_fulltext]
        ft_checked = len(targets_for_ft)
        done = 0
        if targets_for_ft:
            with ThreadPoolExecutor(max_workers=FULLTEXT_WORKERS) as pool:
                futures = {pool.submit(fulltext_fn, p.pmcid): p for p in targets_for_ft}
                for fut in as_completed(futures):
                    paper = futures[fut]
                    done += 1
                    try:
                        results[paper.key] = fut.result()
                    except ApiError:
                        results[paper.key] = None
                        ft_failed += 1
                    report(0.45 + 0.5 * done / len(targets_for_ft), f"전문 {done}/{len(targets_for_ft)} 확인 중…")
        for paper in targets_for_ft:            # 순서를 고정해 결과가 매번 같게 한다
            full = results.get(paper.key)
            if not full:
                continue
            base = f"{paper.title}. {paper.abstract}"
            for cand in extract_candidates(full, min_len=min_len):
                if cand.sequence not in abstract_seqs.get(paper.key, set()):
                    hits.append(SeqHit(cand, paper, SRC_FULLTEXT, _targets(f"{base} {full}", cand.context, cand.sequence),
                                   hit_score(paper, parsed.terms, cand.context, cand.sequence)))

    year_counts = Counter(int(p.year) for p in papers if p.year.isdigit())
    read_keys = {k for k, v in results.items() if v} if use_fulltext else set()
    unread = [p for p in papers if p.key not in read_keys] if use_fulltext else []
    report(1.0, "완료")
    return LiteratureResult(papers, found.hit_count, hits, group_hits(hits), ft_checked, ft_failed,
                            dict(sorted(year_counts.items())), len(read_keys), unread,
                            effective, expanded.notes if expanded else [])
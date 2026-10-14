"""Europe PMC 논문 검색 (메타데이터 + 초록) + 오픈액세스 전문 텍스트.

Europe PMC는 기본 정렬이 최신 논문 쪽으로 치우쳐 있어(상위 100편이 모두 최근 1~2년인 경우가 많다),
전 기간을 보려면 연도 구간별 검색(search_balanced)이나 인용순 정렬을 써야 한다.
"""
from __future__ import annotations

import html
import math
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from datetime import date

from .http import ApiError, get_json, get_text

BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest"
SEARCH_URL = f"{BASE}/search"
EARLIEST_YEAR = 1975        # 이보다 오래된 논문은 초록이 거의 없다
MAX_PAGE = 1000

SORT_RELEVANCE = "관련도 (Europe PMC 기본, 최근 논문 쪽으로 치우침)"
SORT_BALANCED = "전 기간 고르게 (5년 단위 구간별 상위)"
SORT_CITED = "인용 많은 순 (전 기간)"
SORT_NEWEST = "최신순"
SORT_MODES = [SORT_BALANCED, SORT_CITED, SORT_NEWEST, SORT_RELEVANCE]
_SORT_PARAM = {SORT_CITED: "CITED desc", SORT_NEWEST: "P_PDATE_D desc"}


@dataclass
class Paper:
    pmid: str
    pmcid: str
    doi: str
    title: str
    authors: str
    journal: str
    year: str
    abstract: str
    open_access: bool
    source: str = ""        # Europe PMC 출처 코드: MED(PubMed), PMC, PPR(프리프린트), PAT(특허) 등
    ext_id: str = ""        # 그 출처 안의 ID (특허번호 등)

    @property
    def is_patent(self) -> bool:
        return self.source == "PAT"

    @property
    def identifier(self) -> str:
        if self.pmid:
            return f"PMID {self.pmid}"
        if self.pmcid:
            return self.pmcid
        if self.doi:
            return f"DOI {self.doi}"
        if self.source and self.ext_id:
            label = "특허" if self.is_patent else self.source
            return f"{label} {self.ext_id}"
        return ""

    @property
    def key(self) -> str:
        return self.pmid or self.pmcid or self.doi or (f"{self.source}:{self.ext_id}" if self.ext_id else "") or self.title

    @property
    def url(self) -> str:
        if self.pmid:
            return f"https://pubmed.ncbi.nlm.nih.gov/{self.pmid}/"
        if self.pmcid:
            return f"https://europepmc.org/article/PMC/{self.pmcid}"
        if self.doi:
            return f"https://doi.org/{self.doi}"
        if self.source and self.ext_id:
            return f"https://europepmc.org/article/{self.source}/{self.ext_id}"
        return ""


@dataclass
class SearchResult:
    papers: list[Paper]
    hit_count: int          # Europe PMC에 있는 전체 일치 논문 수 (가져온 수와 다를 수 있음)


def _strip_tags(text: str) -> str:
    return html.unescape(re.sub(r"<[^>]+>", "", text or "")).strip()


def parse_results(payload: dict) -> list[Paper]:
    papers = []
    for r in payload.get("resultList", {}).get("result", []):
        journal = r.get("journalInfo", {}).get("journal", {}).get("title", "") or r.get("journalTitle", "")
        papers.append(Paper(
            pmid=str(r.get("pmid", "") or ""),
            pmcid=str(r.get("pmcid", "") or ""),
            doi=str(r.get("doi", "") or ""),
            title=_strip_tags(r.get("title", "")),
            authors=r.get("authorString", ""),
            journal=journal,
            year=str(r.get("pubYear", "")),
            abstract=_strip_tags(r.get("abstractText", "")),
            open_access=r.get("isOpenAccess") == "Y",
            source=str(r.get("source", "") or ""),
            ext_id=str(r.get("id", "") or ""),
        ))
    return papers


def search(query: str, *, page_size: int = 10) -> list[Paper]:
    """논문 제목·키워드 검색 (소량). 서열 탐색에는 search_many/search_balanced를 쓴다."""
    if not query.strip():
        return []
    payload = get_json(
        SEARCH_URL,
        {"query": query, "format": "json", "pageSize": page_size, "resultType": "core"},
        service="Europe PMC",
    )
    return parse_results(payload)


def build_query(query: str, year_from: int | None = None, year_to: int | None = None) -> str:
    q = f"({query.strip()})"
    if year_from is not None or year_to is not None:
        lo = year_from if year_from is not None else EARLIEST_YEAR
        hi = year_to if year_to is not None else date.today().year
        q += f" AND (PUB_YEAR:[{lo} TO {hi}])"
    return q


def _fetch(query: str, count: int, sort: str | None) -> SearchResult:
    """cursorMark 페이지네이션으로 count편을 가져온다."""
    papers: list[Paper] = []
    seen: set[str] = set()
    cursor = "*"
    hit_count = 0
    while len(papers) < count:
        size = min(MAX_PAGE, count - len(papers), 100 if count > 100 else count)
        params = {"query": query, "format": "json", "pageSize": max(size, 1), "resultType": "core", "cursorMark": cursor}
        if sort:
            params["sort"] = sort
        payload = get_json(SEARCH_URL, params, service="Europe PMC", timeout=40)
        hit_count = int(payload.get("hitCount", hit_count) or 0)
        page = parse_results(payload)
        if not page:
            break
        for p in page:
            if p.key not in seen:
                seen.add(p.key)
                papers.append(p)
        nxt = payload.get("nextCursorMark")
        if not nxt or nxt == cursor:
            break
        cursor = nxt
    return SearchResult(papers[:count], hit_count)


def year_buckets(year_from: int, year_to: int, width: int = 5) -> list[tuple[int, int]]:
    """[(시작, 끝)] — 최신 구간부터. 마지막 구간이 가장 오래된 해를 포함한다."""
    out, hi = [], year_to
    while hi >= year_from:
        lo = max(year_from, hi - width + 1)
        out.append((lo, hi))
        hi = lo - 1
    return out


def search_many(query: str, *, max_papers: int = 100, year_from: int | None = None, year_to: int | None = None,
                sort: str = SORT_BALANCED) -> SearchResult:
    """키워드로 최대 max_papers편을 가져온다. sort 모드는 SORT_MODES 참고."""
    if not query.strip():
        return SearchResult([], 0)
    lo = year_from if year_from is not None else EARLIEST_YEAR
    hi = year_to if year_to is not None else date.today().year
    if sort == SORT_BALANCED:
        return search_balanced(query, max_papers=max_papers, year_from=lo, year_to=hi)
    return _fetch(build_query(query, year_from, year_to), max_papers, _SORT_PARAM.get(sort))


def search_balanced(query: str, *, max_papers: int, year_from: int, year_to: int, bucket_years: int = 5) -> SearchResult:
    """연도 구간마다 관련도 상위 논문을 같은 수만큼 가져와 전 기간이 고르게 포함되게 한다."""
    # 전체 일치 논문이 상한 이하이면 구간으로 나누지 않고 전부 가져온다 (순위 밖으로 밀려 놓치는 논문이 없도록)
    probe = _fetch(build_query(query, year_from, year_to), 1, None)
    if 0 < probe.hit_count <= max_papers:
        return _fetch(build_query(query, year_from, year_to), probe.hit_count, None)
    buckets = year_buckets(year_from, year_to, bucket_years)
    per = max(1, math.ceil(max_papers / len(buckets)))
    papers: list[Paper] = []
    seen: set[str] = set()
    total = 0
    for lo, hi in buckets:
        res = _fetch(build_query(query, lo, hi), per, None)
        total += res.hit_count
        for p in res.papers:
            if p.key not in seen:
                seen.add(p.key)
                papers.append(p)
    return SearchResult(papers[:max_papers], total)


# ---- 오픈액세스 전문 ------------------------------------------------------

def xml_to_text(xml: str) -> str:
    """JATS XML → 평문. 참고문헌 목록(ref-list)은 오탐을 줄이려고 제외한다."""
    xml = re.sub(r"<!DOCTYPE[^>]*>", "", xml, count=1)
    try:
        root = ET.fromstring(xml)
        for parent in list(root.iter()):
            for child in list(parent):
                if child.tag == "ref-list":
                    parent.remove(child)
        return re.sub(r"\s+", " ", " ".join(root.itertext())).strip()
    except ET.ParseError:
        stripped = re.sub(r"<ref-list>.*?</ref-list>", " ", xml, flags=re.DOTALL)
        return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", stripped))).strip()


NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
MIN_FULLTEXT_CHARS = 3000       # 이보다 짧으면 본문 없이 서지정보·초록만 있는 것으로 본다


def fetch_fulltext(pmcid: str, *, text_getter=get_text) -> str:
    """PMC 논문 전문 텍스트 (참고문헌 제외). 저작권을 고려해 디스크에 저장(캐시)하지 않는다.

    Europe PMC를 먼저 시도하고, 실패하면 NCBI PMC(efetch)로 다시 시도한다.
    저자 원고(NIH manuscript)처럼 오픈액세스 표시가 없는 논문은 Europe PMC가 500을 돌려주는 경우가 있지만
    NCBI에서는 본문이 받아지는 경우가 있다. 본문이 없는(초록만 있는) 응답은 실패로 취급한다.
    """
    if not pmcid:
        raise ApiError("PMCID가 없어 전문을 가져올 수 없습니다.")
    errors: list[str] = []
    for label, url, params in (
        ("Europe PMC", f"{BASE}/{pmcid}/fullTextXML", None),
        ("NCBI PMC", NCBI_EFETCH, {"db": "pmc", "id": pmcid, "rettype": "xml"}),
    ):
        try:
            xml = text_getter(url, params, ttl=0, service=f"{label} 전문", timeout=60)
        except ApiError as exc:
            errors.append(str(exc))
            continue
        text = xml_to_text(xml)
        if len(text) >= MIN_FULLTEXT_CHARS:
            return text
        errors.append(f"{label}: 본문이 제공되지 않는 논문입니다.")
    raise ApiError("전문을 가져오지 못했습니다. " + " / ".join(errors))

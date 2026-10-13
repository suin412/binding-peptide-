"""UniProt REST 클라이언트: 검색 + 항목 상세(서열, 기능, 참고문헌, 성숙 펩타이드 구간)."""
from __future__ import annotations

import re
from dataclasses import dataclass, field

from .http import get_json

BASE = "https://rest.uniprot.org/uniprotkb"
ACCESSION_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)

PE_LABELS = {
    "1": "단백질 수준 증거", "2": "전사체 수준 증거", "3": "상동성으로 추정",
    "4": "예측됨", "5": "불확실",
}


@dataclass
class SearchHit:
    accession: str
    entry_id: str
    name: str
    organism: str
    length: int


@dataclass
class Region:
    label: str
    start: int
    end: int
    kind: str  # Chain / Peptide


@dataclass
class UniProtEntry:
    accession: str
    entry_id: str
    name: str
    organism: str
    sequence: str
    function: str = ""
    protein_existence: str = ""
    keywords: list[str] = field(default_factory=list)
    pdb_ids: list[str] = field(default_factory=list)
    references: list[dict] = field(default_factory=list)  # {pmid, title, journal, year}
    regions: list[Region] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def url(self) -> str:
        return f"https://www.uniprot.org/uniprotkb/{self.accession}"

    @property
    def direct_sequencing(self) -> bool:
        return "Direct protein sequencing" in self.keywords


def _protein_name(entry: dict) -> str:
    desc = entry.get("proteinDescription", {})
    for key in ("recommendedName",):
        name = desc.get(key, {}).get("fullName", {}).get("value")
        if name:
            return name
    for key in ("submissionNames", "alternativeNames"):
        items = desc.get(key) or []
        if items:
            name = items[0].get("fullName", {}).get("value")
            if name:
                return name
    return entry.get("uniProtkbId", "")


def _query(text: str, reviewed_only: bool) -> str:
    text = text.strip()
    if ACCESSION_RE.match(text.upper()):
        return f"accession:{text.upper()}"
    q = f"({text})"
    return q + " AND reviewed:true" if reviewed_only else q


def search(query: str, *, size: int = 10, reviewed_only: bool = True) -> list[SearchHit]:
    if not query.strip():
        return []
    payload = get_json(
        f"{BASE}/search",
        {"query": _query(query, reviewed_only), "size": size, "format": "json",
         "fields": "accession,id,protein_name,organism_name,length"},
        service="UniProt",
    )
    return parse_search(payload)


def parse_search(payload: dict) -> list[SearchHit]:
    hits = []
    for item in payload.get("results", []):
        hits.append(SearchHit(
            accession=item.get("primaryAccession", ""),
            entry_id=item.get("uniProtkbId", ""),
            name=_protein_name(item),
            organism=item.get("organism", {}).get("scientificName", ""),
            length=item.get("sequence", {}).get("length", 0),
        ))
    return hits


def fetch(accession: str) -> UniProtEntry:
    payload = get_json(f"{BASE}/{accession.strip().upper()}.json", service="UniProt")
    return parse_entry(payload)


def parse_entry(data: dict) -> UniProtEntry:
    seq = data.get("sequence", {}).get("value", "")
    functions = []
    for comment in data.get("comments", []):
        if comment.get("commentType") == "FUNCTION":
            functions += [t.get("value", "") for t in comment.get("texts", [])]

    references = []
    for ref in data.get("references", []):
        cit = ref.get("citation", {})
        pmid = next((x.get("id") for x in cit.get("citationCrossReferences", [])
                     if x.get("database") == "PubMed"), "")
        references.append({
            "pmid": pmid or "",
            "title": cit.get("title", ""),
            "journal": cit.get("journal", ""),
            "year": (cit.get("publicationDate") or "")[:4],
        })

    regions = []
    for feat in data.get("features", []):
        if feat.get("type") in ("Chain", "Peptide"):
            start = feat["location"]["start"].get("value")
            end = feat["location"]["end"].get("value")
            if start and end:
                regions.append(Region(
                    label=feat.get("description") or feat["type"],
                    start=int(start), end=int(end), kind=feat["type"],
                ))

    return UniProtEntry(
        accession=data.get("primaryAccession", ""),
        entry_id=data.get("uniProtkbId", ""),
        name=_protein_name(data),
        organism=data.get("organism", {}).get("scientificName", ""),
        sequence=seq,
        function=" ".join(functions),
        protein_existence=str(data.get("proteinExistence", "")),
        keywords=[k.get("name", "") for k in data.get("keywords", [])],
        pdb_ids=[x["id"] for x in data.get("uniProtKBCrossReferences", []) if x.get("database") == "PDB"],
        references=references,
        regions=regions,
    )

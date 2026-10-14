"""검색어 보정: 재료 약어(PET 등)가 다른 분야 약어와 겹쳐 엉뚱한 논문이 나오는 문제를 줄인다.

예) "PET binding peptide"의 PET은 폴리에틸렌 테레프탈레이트(플라스틱)인데, Europe PMC에서는
양전자방출단층촬영(PET 영상) 논문이 압도적으로 많이 걸린다. 그래서 약어를 정식 명칭으로 확장하고
겹치는 다른 분야의 용어를 제외한다. 사용자가 따옴표·괄호·AND/OR/NOT를 직접 쓰면 건드리지 않는다.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Alias:
    abbr: str
    meaning: str                 # 사용자에게 보여줄 한글 설명
    expansions: tuple[str, ...]  # Europe PMC 검색 구문 (따옴표 포함)
    excludes: tuple[str, ...]    # 겹치는 다른 분야 용어
    exclude_meaning: str


ALIASES: tuple[Alias, ...] = (
    Alias("PET", "폴리에틸렌 테레프탈레이트(플라스틱)",
          ('"polyethylene terephthalate"', '"poly(ethylene terephthalate)"', '"polyethyleneterephthalate"', '"PET plastic"', '"PET film"', '"PET surface"', '"PET microplastic*"',
           '"PET bottle*"', '"PET fiber*"', '"PET substrate*"'),
          ('"positron emission"', "tomograph*", "radiotracer*", "radioligand*", "radiolabel*", '"PET imaging"',
           '"PET/CT"', '"PET/MR*"', "radiopharmaceutical*", '"18F"', '"68Ga"'),
          "양전자방출단층촬영(PET 영상)"),
    Alias("PS", "폴리스티렌",
          ('polystyrene', '"PS surface"', '"PS plate*"', '"PS microplastic*"'),
          ("phosphatidylserine", "photosystem", '"PS-II"'), "포스파티딜세린·광계"),
    Alias("PE", "폴리에틸렌",
          ('polyethylene', '"PE film"', '"PE surface"', '"PE microplastic*"'),
          ('"pulmonary embol*"', "phycoerythrin", "preeclampsia", '"PE-conjugated"'), "폐색전증·피코에리트린"),
    Alias("PP", "폴리프로필렌",
          ('polypropylene', '"PP surface"', '"PP microplastic*"'),
          ('"protein phosphatase"', "phosphopeptide*", '"PP2A"'), "단백질 탈인산화효소"),
    Alias("PC", "폴리카보네이트",
          ('polycarbonate', '"PC surface"'),
          ("phosphatidylcholine", '"prostate cancer"', '"pancreatic cancer"', '"PC12"'), "포스파티딜콜린·암"),
    Alias("PVC", "폴리염화비닐", ('"polyvinyl chloride"', '"poly(vinyl chloride)"', '"PVC surface"'), (), ""),
    Alias("PLA", "폴리락트산",
          ('"polylactic acid"', '"poly(lactic acid)"', '"polylactide"'),
          ("phospholipase", '"PLA2"'), "포스포리파아제"),
    Alias("PMMA", "폴리메틸메타크릴레이트", ('"poly(methyl methacrylate)"', "polymethylmethacrylate", '"PMMA surface"'), (), ""),
    Alias("PDMS", "폴리디메틸실록세인", ("polydimethylsiloxane", '"PDMS surface"'), (), ""),
    Alias("HA", "하이드록시아파타이트",
          ("hydroxyapatite", '"HA surface"', '"HA crystal*"'),
          ("hyaluronic", "hyaluronan", "hemagglutinin", '"HA tag"'), "히알루론산·헤마글루티닌"),
)

_ADVANCED = re.compile(r'["()]|\b(?:AND|OR|NOT)\b')


@dataclass
class ExpandedQuery:
    original: str
    text: str                                   # Europe PMC에 실제로 보낼 검색식
    notes: list[str] = field(default_factory=list)
    terms: list[str] = field(default_factory=list)  # 논문 제목·초록에서 찾을 핵심 용어 (본문 읽기 우선순위 점수용)

    @property
    def changed(self) -> bool:
        return self.text != self.original


def _token_pattern(abbr: str) -> re.Pattern:
    # 앞뒤가 영숫자가 아닌 대문자 약어. 'PET-binding'처럼 뒤에 붙은 하이픈도 함께 소비한다.
    return re.compile(rf"(?<![A-Za-z0-9]){re.escape(abbr)}(?![A-Za-z0-9])-?")


# 핵심 용어에서 제외할 일반 단어 (peptide·binding은 별도 점수로 다룬다)
_GENERIC = {"binding", "bind", "binds", "bound", "peptide", "peptides", "protein", "proteins", "and", "or", "not", "the",
            "for", "of", "to", "on", "in", "with", "affinity", "specific", "selective", "sequence", "sequences",
            "phage", "display", "design", "surface", "surfaces", "material", "materials"}


def _plain_terms(query: str, skip: set[str] = frozenset()) -> list[str]:
    words = re.findall(r"[A-Za-z][A-Za-z0-9\-]{2,}", query)
    return [w.lower() for w in words if w.lower() not in _GENERIC and w not in skip]


def expand_query(query: str) -> ExpandedQuery:
    """재료 약어를 정식 명칭으로 확장하고 겹치는 분야를 제외한 검색식을 만든다."""
    q = (query or "").strip()
    if not q or _ADVANCED.search(q):
        return ExpandedQuery(q, q, terms=_plain_terms(q))
    notes: list[str] = []
    excludes: list[str] = []
    used: list[Alias] = []
    for alias in ALIASES:
        pattern = _token_pattern(alias.abbr)
        if not pattern.search(q):
            continue
        # 이미 정식 명칭을 같이 적었다면 사용자의 의도를 존중한다.
        if re.search(re.escape(alias.expansions[0].strip('"')), q, re.IGNORECASE):
            continue
        group = "(" + " OR ".join(alias.expansions) + ") "
        q = pattern.sub(group, q, count=1).strip()
        used.append(alias)
        note = f"'{alias.abbr}'을(를) {alias.meaning}로 확장했습니다"
        if alias.excludes:
            excludes += list(alias.excludes)
            note += f" ({alias.exclude_meaning} 논문은 제외)"
        notes.append(note)
    if not notes:
        return ExpandedQuery(query.strip(), query.strip(), terms=_plain_terms(query))
    text = f"({q})" + (f" NOT ({' OR '.join(dict.fromkeys(excludes))})" if excludes else "")
    terms = [e.strip('"').replace("*", "").lower() for a in used for e in a.expansions]
    terms += [a.abbr.lower() for a in used]
    terms += _plain_terms(query, skip={a.abbr for a in used})
    return ExpandedQuery(query.strip(), text, notes, list(dict.fromkeys(terms)))

"""소재 키워드 사전으로 결합 대상 '후보'를 제안한다 (확정은 항상 사용자가 한다).

결합 친화도(Kd 등)와 측정 방법은 v1 범위 밖이다.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

# 정식 명칭 → (분류, [(정규식, 대소문자 무시 여부)])
MATERIALS: dict[str, tuple[str, list[tuple[str, bool]]]] = {
    "티타늄 (Ti / TiO₂)": ("금속", [(r"titanium", True), (r"titania", True), (r"TiO2|TiO₂", False), (r"\bTi\b", False)]),
    "금 (Au)": ("금속", [(r"\bgold\b", True), (r"\bAu\b", False), (r"Au\(\d{3}\)", False)]),
    "은 (Ag)": ("금속", [(r"\bsilver\b", True), (r"\bAg\b", False)]),
    "백금 (Pt)": ("금속", [(r"platinum", True), (r"\bPt\b", False)]),
    "스테인리스강": ("금속", [(r"stainless[\s-]steel", True)]),
    "실리카 (SiO₂)": ("무기 소재", [(r"silica", True), (r"SiO2|SiO₂", False), (r"silicon dioxide", True), (r"\bglass\b", True)]),
    "실리콘 (Si)": ("무기 소재", [(r"\bsilicon\b", True), (r"\bSi\b", False)]),
    "하이드록시아파타이트": ("생체 소재", [(r"hydroxyapatite", True), (r"\bHAp?\b", False)]),
    "인산칼슘": ("무기 소재", [(r"calcium phosphate", True)]),
    "산화아연 (ZnO)": ("무기 소재", [(r"zinc oxide", True), (r"\bZnO\b", False)]),
    "산화철 (Fe₃O₄)": ("무기 소재", [(r"iron oxide", True), (r"magnetite", True), (r"Fe3O4", False)]),
    "폴리스티렌": ("고분자", [(r"polystyrene", True)]),
    "PMMA": ("고분자", [(r"\bPMMA\b", False), (r"poly\(methyl methacrylate\)", True)]),
    "PDMS": ("고분자", [(r"\bPDMS\b", False), (r"polydimethylsiloxane", True)]),
    # PET은 의학 영상(양전자방출단층촬영)과 겹치므로 단독 약어는 쓰지 않고 소재 문맥이 있는 표현만 인식한다
    "PET (폴리에틸렌 테레프탈레이트)": ("고분자", [
        (r"polyethylene terephthalate|poly\(ethylene terephthalate\)|polyethyleneterephthalate|polyethylene-terephthalate", True),
        (r"\bPET[- ](?:plastic|film|surface|bottle|fib(?:er|re)|substrate|microplastic|polymer|binding|granule|flake)", False),
    ]),
    "폴리에틸렌": ("고분자", [(r"polyethylene(?! glycol| terephthalate| ?-?terephthalate)", True)]),
    "폴리프로필렌": ("고분자", [(r"polypropylene", True)]),
    "폴리염화비닐 (PVC)": ("고분자", [(r"polyvinyl chloride|poly\(vinyl chloride\)", True), (r"\bPVC\b", False)]),
    "폴리락트산 (PLA)": ("고분자", [(r"polylactic acid|poly\(lactic acid\)|polylactide", True)]),
    "폴리카보네이트": ("고분자", [(r"polycarbonate", True)]),
    "나일론": ("고분자", [(r"\bnylon\b", True), (r"polyamide", True)]),
    "폴리우레탄": ("고분자", [(r"polyurethane", True)]),
    "마이크로플라스틱": ("고분자", [(r"micro-?plastics?", True), (r"nano-?plastics?", True)]),
    "PTFE (테플론)": ("고분자", [(r"\bPTFE\b", False), (r"teflon", True)]),
    "그래핀": ("탄소 소재", [(r"graphene", True), (r"graphite", True)]),
    "탄소나노튜브": ("탄소 소재", [(r"carbon nanotube", True), (r"\bCNTs?\b", False)]),
    "콜라겐": ("생체 소재", [(r"collagen", True)]),
    "키틴/키토산": ("생체 소재", [(r"chitosan", True), (r"chitin", True)]),
    "셀룰로오스": ("생체 소재", [(r"cellulose", True)]),
    "세포막/세포": ("세포·조직", [(r"cell membrane", True), (r"cancer cells?", True), (r"cell surface", True)]),
}

_BINDING_WORDS = re.compile(r"bind|bound|adsor|affinity|attach|immobili[sz]|specific|selected|screen", re.IGNORECASE)
_SENTENCE_SPLIT = re.compile(r"(?<=[.!?])\s+")


@dataclass
class TargetSuggestion:
    name: str
    category: str
    evidence: str       # 근거 문장 1개 (짧은 발췌)
    mentions: int
    binding_context: bool


def _trim(sentence: str, limit: int = 220) -> str:
    sentence = re.sub(r"\s+", " ", sentence).strip()
    return sentence if len(sentence) <= limit else sentence[: limit - 1] + "…"


_BARE_PET = re.compile(r"\bPET\b")


def material_matches_terms(name: str, terms: list[str]) -> bool:
    """MATERIALS 항목(name)이 사용자가 찾는 용어(예: polyethylene terephthalate) 중 하나와 같은 소재인지."""
    entry = MATERIALS.get(name)
    if not entry:
        return False
    for pattern, ignore_case in entry[1]:
        rx = re.compile(pattern, re.IGNORECASE if ignore_case else 0)
        if any(rx.search(t) for t in terms):
            return True
    return False


def suggest_targets(text: str, *, plastic_context: bool = False, anchor: int | None = None) -> list[TargetSuggestion]:
    """텍스트에서 소재 키워드를 찾아 결합 대상 후보를 반환한다. 언급이 없으면 빈 리스트.

    plastic_context=True(논문이 플라스틱·고분자를 다루는 경우)면 단독 약어 `PET`도 폴리에틸렌 테레프탈레이트로 본다.
    기본값에서는 PET 영상(양전자방출단층촬영)과 구분할 수 없으므로 소재 문맥이 있는 표현만 인식한다.
    anchor(텍스트 안의 위치)를 주면 그 위치에 가장 가까이 언급된 소재부터 정렬한다 (표의 행처럼 서열 바로 옆 소재가 대상일 때).
    """
    if not text or not text.strip():
        return []
    sentences = [s for s in _SENTENCE_SPLIT.split(text) if s.strip()]
    results: list[TargetSuggestion] = []
    distances: dict[str, int] = {}
    for name, (category, patterns) in MATERIALS.items():
        compiled = [re.compile(p, re.IGNORECASE if ci else 0) for p, ci in patterns]
        if plastic_context and name.startswith("PET"):
            compiled.append(_BARE_PET)
        hit_sentences = [s for s in sentences if any(c.search(s) for c in compiled)]
        if not hit_sentences:
            continue
        # 결합 관련 단어가 함께 나온 문장을 근거로 우선 사용
        binding_hits = [s for s in hit_sentences if _BINDING_WORDS.search(s)]
        evidence = (binding_hits or hit_sentences)[0]
        mentions = sum(len(c.findall(s)) for s in hit_sentences for c in compiled)
        results.append(TargetSuggestion(name, category, _trim(evidence), mentions, bool(binding_hits)))
        if anchor is not None:
            distances[name] = min((abs(m.start() - anchor) for c in compiled for m in c.finditer(text)), default=10**9)
    if anchor is not None:
        results.sort(key=lambda r: distances.get(r.name, 10**9))
    else:
        results.sort(key=lambda r: (r.binding_context, r.mentions), reverse=True)
    return results

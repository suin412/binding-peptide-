"""논문/특허 텍스트에서 펩타이드 서열 후보를 추출한다.

자동 추출은 오탐이 있을 수 있으므로 결과는 '후보'일 뿐이고, 채택은 항상 사용자가 직접 한다.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

from ..sequence import NONSTANDARD_TO_STANDARD, STANDARD_AA, THREE_TO_ONE

_STD3 = "|".join(k.title() for k in THREE_TO_ONE)
_NON3 = "|".join(k.title() for k in NONSTANDARD_TO_STANDARD)
_TOKEN = rf"(?:[DL]-)?(?:{_STD3}|{_NON3})"
_THREE_RUN = re.compile(rf"\b{_TOKEN}(?:[-–\s]{{1,2}}{_TOKEN}){{3,}}\b", re.IGNORECASE)
_TOKEN_RE = re.compile(rf"([DL]-)?({_STD3}|{_NON3})", re.IGNORECASE)

_ONE_RUN = re.compile(r"(?<![A-Za-z0-9])[ACDEFGHIKLMNPQRSTVWYX]{5,100}(?![A-Za-z0-9])")
_SEQ_ID = re.compile(r"SEQ\.?\s*ID\.?\s*(?:NO|No)\.?s?\.?\s*:?\s*(\d+)", re.IGNORECASE)
_LISTING = re.compile(r"<400>\s*(\d+)\s*(.*?)(?=<\s*2\d0\s*>|<\s*400\s*>|\Z)", re.DOTALL)

# 일반 영단어 중 표준 아미노산 문자로만 이뤄진 흔한 대문자 단어 (오탐 방지용)
_STOPWORDS = {
    "PEPTIDE", "PEPTIDES", "DISPLAY", "PHAGE", "CLASS", "STEEL", "SILVER", "GLASS", "PLASTIC",
    "SPECIES", "SCHEME", "THESIS", "FIGURE", "TABLE", "CHAIN", "LIGHT", "HEAVY", "SHEET",
}


@dataclass
class Candidate:
    sequence: str
    kind: str                       # "1문자" | "3문자" | "서열목록"
    seq_id_no: str = ""
    context: str = ""
    position: int = 0
    warnings: list[str] = field(default_factory=list)
    nonstandard_notes: str = ""
    uncertain: bool = False         # 1문자 나열이 서열 표기인지 일반 약어인지 문맥으로 확신하기 어려운 경우

    @property
    def valid(self) -> bool:
        return bool(self.sequence) and all(c in STANDARD_AA for c in self.sequence)


# 서열 바로 곁에 나오는 단어(강): 가까이 있으면 대문자 나열을 서열로 신뢰한다.
_STRONG_WORDS = re.compile(
    r"peptid|sequence|motif|residue|amino acid|aptamer|\baa\b|\d-?mer\b|SEQ\.?\s*ID|"
    r"polypeptide|epitope|designated|termed|abbreviated|named",
    re.IGNORECASE,
)
# 서열 근처에 흔하지만 약어 주변에도 흔한 단어(약): 더 가까이 있어야 신뢰한다.
_WEAK_WORDS = re.compile(r"clone|phage|display|ligand|bind|bound|fragment|domain|tag\b", re.IGNORECASE)
# 소프트웨어·시약·제조사 이름 주변이면 서열이 아니다 (RIDASCREEN®, CHARMM, "... GmbH, Germany")
_VENDOR_NEAR = re.compile(
    r"GmbH|\bInc\b|\bLtd\b|\bCorp\b|\bCo\.,|Germany|Switzerland|\bUSA\b|\bUK\b|Japan|Taiwan|Sigma|Merck|Thermo|"
    r"software|version|database|\bkit\b|server|package|program|plugin|toolkit",
    re.IGNORECASE,
)
_TRADEMARK_AFTER = re.compile(r"\s?[®™©]")
_TRUST_MIN_LEN = 10

# 대문자 제목·약어·소프트웨어 이름 중 표준 아미노산 문자로만 이뤄진 흔한 것
_NON_SEQUENCE_WORDS = {
    "REFERENCES", "INTEREST", "INTERFERENCE", "MATERIALS", "SAMPLE", "SAMPLES", "IMPACT", "SMILES", "CHARMM",
    "MEDLINE", "PRISMA", "TEGDMA", "RAMPAGE", "NIPALS", "ESTHER", "AVESTIN", "TESCAN", "YASARA", "METTLER",
    "NVIDIA", "CEECIND", "ACTEGA", "ACPYPE", "RIDASCREEN", "INMARE", "NLMPMC", "PARTICLE", "PARTICLES",
    "STATEMENT", "DECLARATION", "AVAILABILITY", "SUPPLEMENTARY", "ACKNOWLEDGEMENTS", "FINANCIAL", "DISCLOSURE",
    # 논문 초록·본문에서 실제로 서열로 잘못 잡힌 약어·이름 (고분자, 유전자, 시약, 소프트웨어, 일반 단어)
    "CRISPR", "PNIPAM", "DMAEMA", "PDMAEMA", "PEGMA", "ADAMTS", "TWEAKR", "SIGNIFICANCE", "FILAMENT", "PERSIA",
    "SPARTAN", "VISTEC", "RTQPCR", "HEPES", "FAPESP", "QIAGEN", "ASSAY", "ISSN", "MENDEP",
}


def _context_trusted(text: str, start: int, end: int, seq: str) -> bool:
    """대문자 나열이 서열일 가능성이 높은지 판단한다 (MEDLINE, CHARMM 같은 약어·이름 걸러내기용).

    신뢰하지 않는 경우: 알려진 비서열 단어 / ®·™가 붙음 / 제조사·소프트웨어 문맥.
    신뢰하는 경우(위에 걸리지 않을 때): 길이가 충분히 김 / 괄호·따옴표로 감싸짐 /
    바로 곁에 peptide·sequence 같은 강한 단어 / 더 가까이에 binding·display 같은 약한 단어.
    """
    if seq in _NON_SEQUENCE_WORDS or _TRADEMARK_AFTER.match(text[end: end + 2]):
        return False
    if _VENDOR_NEAR.search(text[max(0, start - 30): start] + " " + text[end: end + 25]):
        return False
    if len(seq) >= _TRUST_MIN_LEN:
        return True
    before, after = text[max(0, start - 3): start], text[end: end + 2]
    if re.search(r"[(\[\"'“‘]\s*$", before) and re.match(r"\s*[)\]\"'”’,;.]", after):
        return True
    if _STRONG_WORDS.search(text[max(0, start - 60): start] + " " + text[end: end + 30]):
        return True
    return bool(_WEAK_WORDS.search(text[max(0, start - 25): start] + " " + text[end: end + 15]))


def _context(text: str, start: int, end: int, pad: int = 60) -> str:
    snippet = text[max(0, start - pad): end + pad]
    return re.sub(r"\s+", " ", snippet).strip()


def _convert_tokens(tokens: list[tuple[str, str]]) -> tuple[str, list[str]]:
    """[(d접두, 3문자코드)] → (1문자 서열, 비표준 메모 목록)."""
    out, notes = [], []
    for i, (dflag, code) in enumerate(tokens, start=1):
        up = code.upper()
        if up in THREE_TO_ONE:
            out.append(THREE_TO_ONE[up])
            if dflag.upper() == "D-":
                notes.append(f"{i}번째 D-{code.title()} (D형 → 구조 예측에서는 L형으로 취급)")
        else:
            std = NONSTANDARD_TO_STANDARD[up]
            out.append(std)
            notes.append(f"{i}번째 {code.title()} → {std}로 단순화")
    return "".join(out), notes


def _terminal_notes(text: str, start: int, end: int) -> list[str]:
    notes = []
    before = text[max(0, start - 4): start]
    after = text[end: end + 5]
    if re.search(r"Ac-$", before, re.IGNORECASE):
        notes.append("N-말단 아세틸화(Ac) — 구조 예측에서 무시됨")
    if re.match(r"-NH2", after, re.IGNORECASE):
        notes.append("C-말단 아마이드화(-NH2) — 구조 예측에서 무시됨")
    return notes


def _is_heading_like(text: str, start: int, end: int) -> bool:
    """주변이 거의 전부 대문자(제목·표 머리글)면 서열이 아닐 가능성이 높다."""
    # 바로 앞에 대문자 영단어가 둘 이상 이어지면 제목(예: "... ARRAYING INORGANIC PARTICLE")이다.
    # 서열 목록은 보통 쉼표·세미콜론으로 구분되므로 구두점이 붙은 토큰은 해당하지 않는다.
    prev_tokens = text[max(0, start - 60): start].split()[-2:]
    if len(prev_tokens) == 2 and all(re.fullmatch(r"[A-Z][A-Z\-]+", t) for t in prev_tokens) \
            and max(len(t) for t in prev_tokens) >= 4:
        return True
    window = text[max(0, start - 30): start] + text[end: end + 30]
    letters = [c for c in window if c.isalpha()]
    if len(letters) < 8:
        return False
    lower = sum(1 for c in letters if c.islower())
    return lower / len(letters) < 0.15


_KNOWN3 = {*THREE_TO_ONE, *NONSTANDARD_TO_STANDARD}


def _parse_listing_block(block: str) -> tuple[str, list[str]] | None:
    """<400> 블록 본문(3문자 나열, 말단기 포함 가능 또는 1문자 나열) → (1문자 서열, 메모)."""
    cleaned = re.sub(r"(?<![A-Za-z])\d+(?![A-Za-z])", " ", block)  # 줄 번호만 제거 (NH2의 2는 보존)
    tokens = [t for t in re.split(r"[\s\-–]+", cleaned) if t]
    if not tokens:
        return None
    notes: list[str] = []
    body: list[tuple[str, str]] = []
    has_xaa = False
    i = 0
    while i < len(tokens):
        t, up = tokens[i], tokens[i].upper()
        if up == "AC" and not body:
            notes.append("N-말단 아세틸화(Ac) — 구조 예측에서 무시됨")
        elif up in ("NH2", "OH", "H") :
            if up == "NH2":
                notes.append("C-말단 아마이드화(-NH2) — 구조 예측에서 무시됨")
        elif up in ("D", "L") and i + 1 < len(tokens) and tokens[i + 1].upper() in _KNOWN3:
            body.append((f"{up}-", tokens[i + 1]))
            i += 1
        elif up in _KNOWN3:
            body.append(("", t))
        elif up == "XAA":
            has_xaa = True
            body.append(("", "Xaa"))
        else:
            body = []
            break
        i += 1
    if body:
        if has_xaa:
            seq = "".join("X" if c.upper() == "XAA" else THREE_TO_ONE.get(c.upper(), NONSTANDARD_TO_STANDARD.get(c.upper(), "?"))
                          for _, c in body)
            return seq, notes + ["Xaa(미지정 잔기) 포함 — 수정 필요"]
        seq, conv_notes = _convert_tokens(body)
        return seq, conv_notes + notes
    joined = "".join(tokens).upper()
    if joined.isalpha():
        return joined, []
    return None


def extract_candidates(text: str, *, min_len: int = 5) -> list[Candidate]:
    """텍스트에서 서열 후보를 찾아 위치 순으로 반환한다."""
    text = text or ""
    found: list[Candidate] = []
    taken: list[tuple[int, int]] = []

    def overlaps(s: int, e: int) -> bool:
        return any(s < b and a < e for a, b in taken)

    # 1) 특허 서열목록 <400> 블록
    for m in _LISTING.finditer(text):
        parsed = _parse_listing_block(m.group(2))
        if not parsed:
            continue
        seq, notes = parsed
        if len(seq) < min_len:
            continue
        found.append(Candidate(seq, "서열목록", seq_id_no=m.group(1), context=_context(text, m.start(), m.end(), 20),
                               position=m.start(), nonstandard_notes="; ".join(notes)))
        taken.append((m.start(), m.end()))

    # 2) 3문자 표기
    for m in _THREE_RUN.finditer(text):
        if overlaps(m.start(), m.end()):
            continue
        tokens = [(d or "", c) for d, c in _TOKEN_RE.findall(m.group(0))]
        seq, notes = _convert_tokens(tokens)
        notes += _terminal_notes(text, m.start(), m.end())
        found.append(Candidate(seq, "3문자", context=_context(text, m.start(), m.end()),
                               position=m.start(), nonstandard_notes="; ".join(notes)))
        taken.append((m.start(), m.end()))

    # 3) 1문자 표기
    for m in _ONE_RUN.finditer(text):
        s, e = m.start(), m.end()
        seq = m.group(0)
        if overlaps(s, e) or seq in _STOPWORDS or len(seq) < min_len:
            continue
        if len(set(seq)) < 2 or set(seq) <= set("ACGT"):
            continue
        if _is_heading_like(text, s, e):
            continue
        notes = _terminal_notes(text, s, e)
        cand = Candidate(seq, "1문자", context=_context(text, s, e), position=s,
                         nonstandard_notes="; ".join(notes), uncertain=not _context_trusted(text, s, e, seq))
        found.append(cand)
        taken.append((s, e))

    # SEQ ID NO 연결 (서열목록 블록은 이미 번호가 있다)
    mentions = [(m.start(), m.end(), m.group(1)) for m in _SEQ_ID.finditer(text)]
    for cand in found:
        if cand.seq_id_no or not mentions:
            continue
        best, best_dist = None, None
        for ms, me, num in mentions:
            if cand.position >= me:
                dist = cand.position - me
                limit = 120
            else:
                dist = ms - cand.position
                limit = 80
            if dist <= limit and (best_dist is None or dist < best_dist):
                best, best_dist = num, dist
        if best:
            cand.seq_id_no = best

    # 경고 추가 + 중복 제거 (같은 서열이면 SEQ ID NO가 있는 쪽을 우선)
    uniq: dict[str, Candidate] = {}
    for cand in sorted(found, key=lambda c: c.position):
        if "X" in cand.sequence or "?" in cand.sequence:
            cand.warnings.append("표준 20 아미노산이 아닌 문자(X 등)가 있어 수정 후 채택해야 합니다.")
        if cand.uncertain:
            cand.warnings.append("주변에 서열임을 알려주는 단어가 없어 일반 약어(예: MEDLINE)일 수 있습니다. 문맥을 확인하세요.")
        elif len(cand.sequence) < 8:
            cand.warnings.append("짧은 서열이라 일반 단어·약어와 혼동됐을 수 있습니다.")
        if cand.nonstandard_notes:
            cand.warnings.append("비표준 요소가 있어 표준 서열로 단순화했습니다 (원문 표기는 메모 참고).")
        prev = uniq.get(cand.sequence)
        if prev is None or (not prev.seq_id_no and cand.seq_id_no) or (prev.uncertain and not cand.uncertain):
            uniq[cand.sequence] = cand
    return sorted(uniq.values(), key=lambda c: c.position)

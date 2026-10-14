"""서열 + 출처 + 발굴 방법 + 결합 대상 메타데이터 (SeqRecord)."""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import date

from .affinity import Affinity

NOT_REPORTED = "미보고"

# 발굴(스크리닝) 방법 선택 목록 — 사용자가 직접 고른다. 자동 추출은 하지 않는다.
SCREENING_METHODS = [
    "파지 디스플레이 (biopanning)",
    "효모·박테리아 표면 디스플레이",
    "mRNA·리보솜 디스플레이",
    "DEL·OBOC 비드 라이브러리",
    "합성 펩타이드 라이브러리 (어레이, 알라닌 스캐닝 등)",
    "기능 기반 HTS (MIC, 세포·효소 어세이)",
    "천연물 분획 (bioactivity-guided)",
    "인실리코 스크리닝",
    "기타·불명",
]

BINDING_CATEGORIES = [
    "금속", "무기 소재", "고분자", "탄소 소재", "생체 소재",
    "세포·조직", "단백질·핵산", "기타",
]

# 결합 대상 상태
STATUS_NONE = "미보고"
STATUS_SUGGESTED = "제안"
STATUS_CONFIRMED = "확정"

SOURCE_UNIPROT = "UniProt"
SOURCE_PAPER = "논문"
SOURCE_PATENT = "특허"
SOURCE_MANUAL = "수기 입력"
SOURCE_PRESET = "프리셋"
SOURCE_DESIGNED = "설계"


@dataclass
class Provenance:
    source_type: str = SOURCE_MANUAL
    identifier: str = ""          # UniProt AC, PMID/DOI, 특허번호 + SEQ ID NO 등
    url: str = ""
    retrieved: str = field(default_factory=lambda: date.today().isoformat())
    title: str = ""

    def short(self) -> str:
        if self.identifier:
            return f"{self.source_type} · {self.identifier}"
        return self.source_type


@dataclass
class BindingTarget:
    name: str = ""
    category: str = ""
    note: str = ""
    status: str = STATUS_NONE

    @property
    def reported(self) -> bool:
        return bool(self.name.strip())

    def label(self) -> str:
        if not self.reported:
            return NOT_REPORTED
        tag = "" if self.status == STATUS_CONFIRMED else f" ({self.status})"
        return f"{self.name}{tag}"


@dataclass
class SeqRecord:
    name: str
    sequence: str
    provenance: Provenance = field(default_factory=Provenance)
    screening_methods: list[str] = field(default_factory=list)
    screening_note: str = ""
    binding_target: BindingTarget = field(default_factory=BindingTarget)
    nonstandard_notes: str = ""
    affinities: list[Affinity] = field(default_factory=list)   # 문헌에서 보고된 결합 친화도 (예측값이 아님)
    parent: str | None = None          # 부모 SeqRecord 이름
    mutation: str = ""                 # 변이 표기 (예: K12A)
    is_designed: bool = False
    reference_info: dict = field(default_factory=dict)  # UniProt PE, 키워드 등 "참고 정보"

    @property
    def length(self) -> int:
        return len(self.sequence)

    def affinity_label(self) -> str:
        """대표 친화도 한 줄 (측정값 우선, 없으면 계산값, 아무것도 없으면 미보고)."""
        if not self.affinities:
            return NOT_REPORTED
        ordered = sorted(self.affinities, key=lambda a: (a.computed, a.status != "확정"))
        extra = f" 외 {len(ordered) - 1}건" if len(ordered) > 1 else ""
        return ordered[0].label() + extra

    def methods_label(self) -> str:
        return ", ".join(self.screening_methods) if self.screening_methods else NOT_REPORTED

    def summary(self, parent: "SeqRecord | None" = None) -> str:
        """서열 헤더에 항상 표시하는 요약 문장 (마크다운)."""
        if self.is_designed:
            base = "**설계된 서열 (비실험)** — 결합·활성은 측정되지 않았습니다."
            if parent is not None:
                base += (
                    f" 부모 `{parent.name}` 참고 정보: 결합 대상 {parent.binding_target.label()}, "
                    f"발굴 방법 {parent.methods_label()} (변이가 이를 유지한다는 뜻이 아님)."
                )
            return base
        # 굵게 표시 뒤에 한글 조사를 바로 붙이면 마크다운이 깨지므로(**...)에 처럼) 항상 구분 기호로 끊는다
        affinity = self.affinity_label()
        return (f"결합 대상: **{self.binding_target.label()}** · 발굴 방법: **{self.methods_label()}** · "
                f"친화도: **{affinity}**")

    def make_variant(self, name: str, sequence: str, mutation: str = "") -> "SeqRecord":
        """변이체 생성: 부모 참조를 설정하고 출처는 '설계'로 표시하되 부모 식별자를 남긴다."""
        prov = Provenance(
            source_type=SOURCE_DESIGNED,
            identifier=f"부모 {self.name}" + (f" · {mutation}" if mutation else ""),
            url=self.provenance.url,
            title=self.provenance.title,
        )
        return SeqRecord(
            name=name,
            sequence=sequence,
            provenance=prov,
            nonstandard_notes=self.nonstandard_notes,
            parent=self.name,
            mutation=mutation,
            is_designed=True,
        )

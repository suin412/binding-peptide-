"""여러 탭이 공유하는 UI 조각: 서열 헤더, 발굴 방법·결합 대상 입력."""
from __future__ import annotations

import streamlit as st

from ..provenance import (
    BINDING_CATEGORIES, NOT_REPORTED, SCREENING_METHODS, STATUS_CONFIRMED, STATUS_NONE, STATUS_SUGGESTED,
    BindingTarget, SeqRecord,
)
from ..sources.binding import TargetSuggestion
from .state import records

DIRECT = "직접 입력"
NONE_OPTION = "(미보고)"


def excel_download_button(key: str, *, target=st, label: str = "📥 전체 결과 엑셀(.xlsx) 내보내기") -> None:
    """지금까지의 모든 결과(내 서열·친화도·점수·모든 추천 실행 기록·논문 검색 결과)를 한 엑셀 파일로 내려받는다."""
    from datetime import datetime

    from ..export import build_export

    data = build_export(records(), st.session_state.get("rec_result"), st.session_state.get("imp_result"),
                       runs=st.session_state.get("run_log") or None, lit_state=st.session_state.get("lit_result"))
    target.download_button(label, data, file_name=f"peplab_결과_{datetime.now():%Y%m%d_%H%M}.xlsx", key=key,
                           mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                           help="한 파일에 모든 결과가 시트로 들어갑니다: 안내 · 내 서열 · 결합 성향 점수 · 결합 친화도 · 예측 vs 보고값 · 개선형 추천(실행 기록 전체) · 전체 서열 기반 추천(실행 기록 전체) · 참고한 알려진 서열 · 논문 검색 결과")


def quick_add_form(key: str, *, expanded: bool = False) -> None:
    """분석 탭 위에 두는 '새 서열 바로 입력' 폼: 저장하면 바로 분석 대상이 된다."""
    from ..provenance import SOURCE_MANUAL, Provenance
    from ..sequence import SequenceError, validate_sequence
    from .state import add_record, flash, unique_name

    with st.expander("➕ 새 서열 직접 입력 (저장하면 바로 분석·개선형 추천 가능)", expanded=expanded):
        seq_text = st.text_area("서열 (1문자 표기, 표준 20 아미노산)", height=70, key=f"qa_seq_{key}",
                                placeholder="예: WEWWHVFHHRLR")
        name = st.text_input("이름 (비워 두면 서열로 정함)", key=f"qa_name_{key}")
        seq_ok = None
        if seq_text.strip():
            try:
                seq_ok = validate_sequence(seq_text)
                st.caption(f"{len(seq_ok)} aa — 검증 통과")
            except SequenceError as exc:
                st.error(str(exc))
        if st.button("내 서열에 추가", type="primary", key=f"qa_add_{key}", disabled=seq_ok is None):
            rec = SeqRecord(name=unique_name(name.strip() or seq_ok[:14]), sequence=seq_ok,
                            provenance=Provenance(SOURCE_MANUAL, "직접 입력"))
            add_record(rec)
            flash(f"'{rec.name}'을(를) 추가했습니다. 개선형 추천·서열 분석 탭에서 바로 사용할 수 있습니다.")
            st.rerun()


def render_header(rec: SeqRecord) -> None:
    """서열 이름·길이·출처와 요약 문장을 항상 보여준다."""
    parent = records().get(rec.parent) if rec.parent else None
    badge = "🧪 설계된 서열 (비실험)" if rec.is_designed else f"📄 {rec.provenance.short()}"
    st.markdown(f"### {rec.name}")
    st.caption(f"{rec.length} aa · {badge}" + (f" · 변이 {rec.mutation}" if rec.mutation else ""))
    st.markdown(rec.summary(parent))
    if rec.nonstandard_notes:
        st.warning(f"비표준 요소 메모: {rec.nonstandard_notes}  \n표준 서열로 단순화되어 있으며, 분석·추천에서는 화학 변형이 무시됩니다.")
    with st.expander("서열 · 출처 · 참고 정보"):
        st.code(rec.sequence, language=None, wrap_lines=True)
        prov = rec.provenance
        if prov.url:
            st.markdown(f"출처: [{prov.short()}]({prov.url}) · 확인일 {prov.retrieved}")
        else:
            st.markdown(f"출처: {prov.short()} · 확인일 {prov.retrieved}")
        if prov.title:
            st.caption(prov.title)
        ref = rec.reference_info
        if ref:
            lines = []
            if ref.get("protein_existence"):
                lines.append(f"UniProt 단백질 존재 수준: {ref['protein_existence']}")
            if ref.get("direct_sequencing"):
                lines.append("키워드 'Direct protein sequencing' 있음")
            if ref.get("pdb_ids"):
                lines.append("PDB: " + ", ".join(ref["pdb_ids"][:8]))
            if ref.get("improved_from"):
                lines.append(f"개선형 추천: '{ref['improved_from']}'에서 변이")
            if ref.get("recommended_for"):
                lines.append(f"자동 추천 대상: {ref['recommended_for']} (예상 결합 성향 {ref.get('recommend_score')}, 학습용 점수)")
            if lines:
                st.caption("참고 정보 (발굴 방법을 대신하지 않음) — " + " · ".join(lines))
        if rec.screening_note:
            st.caption(f"발굴 방법 메모: {rec.screening_note}")
        if rec.binding_target.note:
            st.caption(f"결합 대상 메모: {rec.binding_target.note}")


def context_editor(key: str, *, suggestions: list[TargetSuggestion] | None = None,
                   initial: SeqRecord | None = None) -> tuple[list[str], str, BindingTarget]:
    """발굴 방법(선택 목록) + 결합 대상(제안 → 사용자 확정) 입력. 반환값으로 레코드를 채운다."""
    suggestions = suggestions or []
    init_methods = initial.screening_methods if initial else []
    methods = st.multiselect(
        "발굴(스크리닝) 방법 — 논문을 읽고 직접 선택", SCREENING_METHODS, default=init_methods,
        key=f"{key}_methods", placeholder="선택하지 않으면 '미보고'로 표시됩니다",
    )
    note = st.text_input("발굴 방법 메모 (라운드 수, 선택 압력 등)", value=initial.screening_note if initial else "",
                         key=f"{key}_method_note")

    st.markdown("**결합 대상**")
    init_target = initial.binding_target if initial else BindingTarget()
    options = [NONE_OPTION] + [s.name for s in suggestions] + [DIRECT]
    default_index = 0
    if init_target.reported:
        default_index = options.index(init_target.name) if init_target.name in options else len(options) - 1
    choice = st.selectbox(
        "논문 텍스트에서 제안된 후보" if suggestions else "결합 대상", options, index=default_index, key=f"{key}_target_choice",
    )
    name, category, evidence = "", "", ""
    if choice == DIRECT:
        name = st.text_input("결합 대상 이름 (예: 티타늄, 금 나노입자, 특정 수용체)",
                             value=init_target.name if init_target.reported and init_target.name not in options else "",
                             key=f"{key}_target_name")
        category = st.selectbox("분류", BINDING_CATEGORIES, key=f"{key}_target_cat",
                                index=BINDING_CATEGORIES.index(init_target.category) if init_target.category in BINDING_CATEGORIES else 0)
    elif choice != NONE_OPTION:
        sug = next(s for s in suggestions if s.name == choice) if suggestions else None
        name = choice if sug is None else sug.name
        category = sug.category if sug else init_target.category
        if sug:
            st.caption(f"근거 문장: “{sug.evidence}”")
            st.caption(f"분류: {sug.category} · 언급 {sug.mentions}회 — 소재 키워드 매칭에 의한 **제안**이며 정확하지 않을 수 있습니다.")
    tnote = ""
    status = STATUS_NONE
    if name.strip():
        tnote = st.text_input("형태·조건 메모 (박막/나노입자, pH, 완충액 등 — 논문에 있는 것만)",
                              value=init_target.note, key=f"{key}_target_note")
        confirmed = st.checkbox("논문에서 직접 확인했으므로 '확정'으로 기록", value=(choice == DIRECT or init_target.status == STATUS_CONFIRMED),
                                key=f"{key}_target_confirm")
        status = STATUS_CONFIRMED if confirmed else STATUS_SUGGESTED
    else:
        st.caption(f"결합 대상은 '{NOT_REPORTED}'로 기록됩니다 (추측으로 채우지 않습니다).")
    return methods, note, BindingTarget(name.strip(), category, tnote, status)


def apply_context(rec: SeqRecord, methods: list[str], note: str, target: BindingTarget) -> None:
    rec.screening_methods = methods
    rec.screening_note = note
    rec.binding_target = target

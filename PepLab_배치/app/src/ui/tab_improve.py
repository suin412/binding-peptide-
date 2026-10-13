"""개선형 추천 탭: 채택한 서열 하나를 고르면 그 서열의 몇 개 잔기만 바꾼 개선형을 자동으로 추천한다."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from ..analysis import plots
from ..analysis.binding_score import CLASSES, class_for_material
from ..design.improve import improve, parse_positions
from ..provenance import SOURCE_MANUAL, Provenance, SeqRecord
from ..sequence import SequenceError, validate_sequence
from .common import excel_download_button
from .state import add_record, current_record, flash, log_run, records, set_current, unique_name
from .structure_view import properties_table
from .tab_recommend import reference_options

DISCLAIMER = ("**점수가 오르는 방향의 변이를 제안할 뿐, 실제 친화도(Kd) 개선을 예측하지 않습니다.** 표면 가중치·알려진 서열과의 유사도·"
              "용해도 제약으로 점수를 매기고, 바꾼 잔기 하나마다 2점의 비용을 물려 적게 바꾸는 쪽을 우선합니다. "
              "서열의 순서·입체 구조·pH는 거의 반영하지 못하니 후보를 실험으로 확인하세요.")


MODE_RECORDS = "🔎 검색해서 채택한 서열에서 선택"
MODE_DIRECT = "✍️ 서열 직접 입력"


def render() -> None:
    recs = records()
    st.markdown("**검색해서 채택한 서열 하나**를 고르거나 **서열을 직접 입력**하면, 그 서열에서 **몇 개 잔기만 바꿔** "
                "결합 성향 점수를 높일 수 있는 **개선형**을 자동으로 추천합니다.")
    st.warning(DISCLAIMER)
    modes = [MODE_RECORDS, MODE_DIRECT]
    mode = st.radio("개선할 서열 입력 방식", modes, index=0 if recs else 1, horizontal=True,
                    key=f"imp_mode_{'r' if recs else 'e'}")
    if mode == MODE_RECORDS:
        parent = _pick_from_records(recs)
    else:
        parent = _type_sequence()
    if parent is None:
        return
    parent_name = parent.name
    try:
        validate_sequence(parent.sequence, min_len=5)
        if len(parent.sequence) > 40:
            raise SequenceError("개선형 추천은 40 aa 이하 서열만 지원합니다.")
        valid = True
    except SequenceError as exc:
        st.error(str(exc))
        valid = False
    if parent.is_designed:
        st.caption("설계된 서열(비실험)입니다. 개선형도 같은 방식으로 추천하지만 둘 다 측정된 서열이 아닙니다.")

    guess = class_for_material(parent.binding_target.name) if parent.binding_target.reported else None
    keys = list(CLASSES)
    class_key = st.selectbox("결합 대상 (표면 종류)", keys, index=keys.index(guess) if guess else 0,
                             format_func=lambda k: CLASSES[k].label, key=f"imp_class_{parent_name}_{guess}")
    if guess:
        st.caption(f"이 서열의 결합 대상 '{parent.binding_target.name}' → {CLASSES[guess].label.split(' (')[0]}")
    else:
        st.caption("이 서열에 결합 대상이 기록돼 있지 않거나 표면 분류를 추정하지 못했습니다. 개선할 표면을 직접 골라 주세요.")
    cls = CLASSES[class_key]
    st.caption(f"가정: {cls.rationale}")

    c1, c2 = st.columns(2)
    max_mut = c1.slider("최대 변이 수 (바꿀 잔기 개수)", 1, 5, 3, key="imp_maxmut")
    n = c2.slider("추천 개수", 3, 15, 8, key="imp_n")
    fixed_text = st.text_input("바꾸지 않을 위치 (선택, 예: 1-3, 7)", key=f"imp_fixed_{parent_name}",
                               help="알려진 활성 모티프 등 지키고 싶은 위치를 적으면 그 위치는 변이시키지 않습니다. 1부터 센 위치입니다.")
    fixed: list[int] = []
    fixed_ok = True
    if fixed_text.strip():
        try:
            fixed = parse_positions(fixed_text, len(parent.sequence))
            st.caption(f"고정 위치: {', '.join(f'{parent.sequence[p - 1]}{p}' for p in fixed)}")
        except ValueError as exc:
            st.error(str(exc))
            fixed_ok = False

    options = {k: s for k, s in reference_options(class_key, False).items() if s.sequence != parent.sequence}
    with st.expander("참고할 알려진 서열 · 고급 설정"):
        refs = st.multiselect(f"알려진 서열 {len(options)}개 (개선 방향이 이 서열들과 닮도록 가산)" if options else "참고할 알려진 서열이 없습니다 (표면 가정 점수만 사용)",
                              list(options), default=list(options), key=f"imp_refs_{class_key}_{len(options)}_{parent_name}", disabled=not options)
        weight = st.slider("알려진 서열과 닮은 정도의 반영 비중", 0.0, 0.8, 0.2, step=0.05, key="imp_sim")

    run = st.button("개선형 추천 받기", type="primary", key="imp_go", disabled=not (valid and fixed_ok))
    if run:
        with st.spinner("개선형을 탐색하는 중…"):
            result = improve(parent.sequence, class_key, n=n, max_mutations=max_mut, fixed_positions=fixed,
                             seeds=[options[k] for k in refs], use_literature_seeds=False, similarity_weight=weight)
        st.session_state["imp_result"] = {"result": result, "parent": parent_name, "label": cls.label, "parent_rec": parent}
        log_run("imp", st.session_state["imp_result"])
    state = st.session_state.get("imp_result")
    if state and state["parent"] == parent_name and state["result"].parent == parent.sequence:
        _results(state)


def _pick_from_records(recs: dict[str, SeqRecord]) -> SeqRecord | None:
    """검색·채택해 '내 서열'에 있는 서열 중에서 고른다."""
    if not recs:
        st.info("아직 채택한 서열이 없습니다. '검색' 탭에서 서열을 채택하거나, 위에서 '직접 입력'을 고르세요.")
        return None
    names = list(recs)
    cur = current_record()
    name = st.selectbox("개선할 서열 (검색해서 채택한 내 서열)", names, index=names.index(cur.name) if cur else 0, key="imp_parent",
                        format_func=lambda n: ("🧪 " if recs[n].is_designed else "📄 ") + n)
    parent = recs[name]
    st.code(parent.sequence, language=None, wrap_lines=True)
    if not parent.is_designed:
        st.caption(f"출처: {parent.provenance.short()}")
    return parent


def _type_sequence() -> SeqRecord | None:
    """서열을 직접 입력한다. 저장하지 않아도 추천을 받을 수 있고, 개선형을 저장할 때 이 서열도 부모로 함께 저장된다."""
    seq_text = st.text_area("개선할 서열 (1문자 표기, 표준 20 아미노산, 5~40 aa)", height=70, key="imp_direct_seq", placeholder="예: RKLPDAPGMHTW")
    name = st.text_input("이름 (선택, 비워 두면 서열로 정함)", key="imp_direct_name")
    if not seq_text.strip():
        st.info("개선할 서열을 입력하세요. 입력만으로는 '내 서열'에 저장되지 않습니다.")
        return None
    try:
        seq = validate_sequence(seq_text)
    except SequenceError as exc:
        st.error(str(exc))
        return None
    st.caption(f"{len(seq)} aa — 검증 통과 (직접 입력, 출처 미기록)")
    return SeqRecord(name=name.strip() or seq[:14], sequence=seq, provenance=Provenance(SOURCE_MANUAL, "직접 입력"))


def _results(state: dict) -> None:
    result = state["result"]
    parent_rec = records().get(state["parent"])
    if parent_rec is not None and parent_rec.sequence != result.parent:
        parent_rec = None                                                      # 같은 이름의 다른 서열이면 직접 입력한 서열로 본다
    parent_rec = parent_rec or state.get("parent_rec")
    cls = CLASSES[result.class_key]
    m1, m2, m3 = st.columns(3)
    m1.metric("부모 점수", f"{result.parent_score:.1f}", help="표면·유사도 점수에서 제약 위반 벌점을 뺀 값 (벌점이 크면 음수)")
    m2.metric("부모 표면 점수", f"{result.parent_binding:.1f}")
    m3.metric("탐색한 변이체", f"{result.searched:,}개")
    if result.parent_notes:
        st.caption("부모의 제약 위반: " + "; ".join(result.parent_notes))
    if result.fixed_positions:
        st.caption(f"고정한 위치: {', '.join(map(str, result.fixed_positions))}")

    if result.already_good:
        st.info("점수가 부모보다 오르는 변이를 찾지 못했습니다. 이 서열은 **이 점수 기준으로 이미 (국소) 최적**에 가깝거나, 고정한 위치가 많아 "
                "바꿀 수 있는 곳이 부족합니다. 최대 변이 수를 늘리거나 고정 위치를 줄여 보세요. (점수가 낮다고 실제로 나쁜 서열이라는 뜻은 아닙니다.)")
        return
    st.markdown(f"#### 개선형 {len(result.candidates)}개 — {cls.label.split(' (')[0]}")
    rows = []
    for i, c in enumerate(result.candidates):
        rows.append({
            "순위": i + 1, "서열": c.sequence, "변경": " ".join(c.mutations), "변이 수": len(c.mutations),
            "예상 결합 성향": c.score, "부모 대비": f"{c.delta:+.1f}", "표면 점수": c.binding,
            "순전하": c.charge, "GRAVY": c.gravy, "분자량": c.mw, "비고": "; ".join(c.notes),
        })
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch", height=min(420, 36 * len(rows) + 40))
    st.caption("점수에는 변이 1개당 2점의 비용이 반영돼 있습니다. '부모 대비'는 부모 점수와의 차이이고, 같은 변이를 여러 개 겹친 후보가 "
               "순위를 독차지하지 않게 서로 3곳 이상 다른 후보를 골랐습니다. 같은 입력이면 결과가 재현됩니다.")

    excel_download_button("imp_xlsx")
    labels = {i: f"#{i + 1} {' '.join(c.mutations)} → {c.sequence} ({c.delta:+.1f})" for i, c in enumerate(result.candidates)}
    pick = st.selectbox("변이별 표면 가중치 변화 · 잔기 기여 보기", list(labels), format_func=labels.get, key="imp_pick")
    cand = result.candidates[pick]
    st.write("변이별 표면 가중치 변화(+이면 이 표면에 유리한 잔기로 바뀜): " +
             ", ".join(f"{lab} ({chg:+.2f})" for lab, chg in cand.weight_changes.items()))
    st.pyplot(plots.residue_contribution_plot(cand.sequence, cand.contributions, f"{' '.join(cand.mutations)} — {cls.label}"))

    st.markdown("##### 물성 분석 (저장하지 않아도 바로 볼 수 있습니다)")
    items = [("부모", result.parent)] + [(f"#{i + 1} {' '.join(c.mutations)}", c.sequence) for i, c in enumerate(result.candidates)]
    with st.expander("물성 비교 (부모 + 전체 개선형)", expanded=True):
        st.dataframe(properties_table(items, label_column="이름"), hide_index=True, width="stretch")

    chosen = st.multiselect("내 서열에 추가할 개선형 (여러 개 선택 가능)", list(labels), format_func=labels.get, key="imp_add_pick",
                            placeholder="개선형을 선택하세요")
    if st.button(f"선택한 {len(chosen)}개를 내 서열에 추가", disabled=not chosen or parent_rec is None, key="imp_add"):
        first = None
        if parent_rec.name not in records() or records()[parent_rec.name].sequence != parent_rec.sequence:   # 직접 입력한 서열은 부모로 함께 저장
            same = next((r for r in records().values() if r.sequence == parent_rec.sequence), None)
            if same is not None:
                parent_rec = same
            else:
                parent_rec.name = unique_name(parent_rec.name)
                add_record(parent_rec, make_current=False)
        for i in chosen:
            c = result.candidates[i]
            mutation = "+".join(c.mutations)
            var = parent_rec.make_variant(unique_name(f"{parent_rec.name}_{mutation}"), c.sequence, mutation)
            var.reference_info = {"recommended_for": cls.label, "recommend_score": c.score, "improved_from": parent_rec.name}
            add_record(var, make_current=False)
            first = first or var.name
        set_current(first)
        flash(f"{len(chosen)}개를 '{parent_rec.name}'의 개선형(설계된 서열, 비실험)으로 추가했습니다. 서열 분석 탭에서 부모와 정렬·비교해 보세요.")
        st.rerun()

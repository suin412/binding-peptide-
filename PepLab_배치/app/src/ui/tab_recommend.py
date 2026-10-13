"""전체 서열 기반 추천 탭: 알려진 결합 서열 전체를 참고해 결합 성향이 높을 것 같은 새 서열을 자동으로 추천한다."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from ..analysis import plots
from ..analysis.binding_score import CLASSES, class_for_material, score_sequence
from ..design.recommend import CLASS_CONSTRAINTS, LITERATURE_SEEDS, Recommendation, SeedSequence, recommend
from ..provenance import SOURCE_DESIGNED, Provenance, SeqRecord
from ..sequence import SequenceError, validate_sequence
from .common import excel_download_button
from .state import add_record, current_record, flash, log_run, records, set_current, unique_name
from .structure_view import properties_table

DISCLAIMER = ("**예상 결합 성향입니다 — 실제 친화도(Kd)의 예측이 아닙니다.** 표면 화학별 정성적 경향(소수성 고분자 ← 방향족 잔기, "
              "티타늄·실리카 ← 양전하 잔기 등)과, **알려진 결합 서열 전체**와의 유사도(잔기 구성 + 연속한 두 잔기 모티프)로 점수를 매기고, "
              "용해도 등 현실적인 제약을 벌점으로 반영해 후보를 고릅니다. 서열의 입체 구조·pH·이온 세기는 반영하지 못하므로 추천 서열은 "
              "실험 후보를 좁히는 출발점으로만 쓰세요.")


def reference_options(class_key: str, include_unlabeled: bool) -> dict[str, SeedSequence]:
    """참고할 수 있는 '알려진 서열' 전체: 문헌 예시 + 내 서열(같은 표면에 결합한다고 기록) (+ 선택: 결합 대상을 기록하지 않은 내 서열)."""
    out: dict[str, SeedSequence] = {}
    for s in LITERATURE_SEEDS.get(class_key, []):
        out[f"문헌 · {s.name} · {s.sequence}"] = s
    for rec in records().values():
        if rec.is_designed:
            continue
        mapped = class_for_material(rec.binding_target.name) if rec.binding_target.reported else None
        labeled = mapped == class_key
        if not (labeled or (include_unlabeled and not rec.binding_target.reported)):
            continue
        try:
            seq = validate_sequence(rec.sequence, min_len=5)
        except SequenceError:
            continue
        if len(seq) > 40 or any(s.sequence == seq for s in out.values()):
            continue
        tag = "내 서열" if labeled else "내 서열(결합 대상 미기록)"
        out[f"{tag} · {rec.name} · {seq}"] = SeedSequence(rec.name, seq, rec.provenance.short(),
                                                           from_literature=labeled and rec.provenance.source_type != "수기 입력")
    return out


def render() -> None:
    st.markdown("결합할 **대상(표면)** 을 고르면, **알려진 결합 서열 전체**를 참고해 결합 성향이 높을 것 같은 **새 서열을 자동으로 추천**합니다. "
                "특정 서열 하나를 고쳐 쓰고 싶으면 🔧 개선형 추천 탭을 쓰세요.")
    st.warning(DISCLAIMER)

    cur = current_record()
    guess = class_for_material(cur.binding_target.name) if cur and cur.binding_target.reported else None
    keys = list(CLASSES)
    class_key = st.selectbox("결합 대상 (표면 종류)", keys, index=keys.index(guess) if guess else 0,
                             format_func=lambda k: CLASSES[k].label, key=f"rec_class_{cur.name if cur else ''}_{guess}")
    if guess:
        st.caption(f"현재 서열 '{cur.name}'의 결합 대상 '{cur.binding_target.name}' → {CLASSES[guess].label.split(' (')[0]}")
    cls = CLASSES[class_key]
    st.caption(f"가정: {cls.rationale}")

    c1, c2 = st.columns(2)
    length = c1.slider("서열 길이 (aa)", 6, 25, (10, 14), key="rec_len")
    n = c2.slider("추천 개수", 3, 30, 10, key="rec_n")

    st.markdown("##### 참고할 알려진 서열")
    include_unlabeled = st.checkbox("결합 대상을 기록하지 않은 내 서열도 후보로 보기", key=f"rec_unlabeled_{class_key}",
                                    help="결합한다고 확인되지 않은 서열이 섞일 수 있으니, 결합 서열이 맞는 것만 골라 쓰세요.")
    options = reference_options(class_key, include_unlabeled)
    chosen_refs = st.multiselect(
        f"알려진 서열 {len(options)}개 중 추천에 참고할 서열 (기본: 전부)" if options else "참고할 알려진 서열이 없습니다 (표면 가정 점수만 사용)",
        list(options), default=list(options), key=f"rec_refs_{class_key}_{len(options)}_{int(include_unlabeled)}", disabled=not options,
    )
    if options:
        st.caption("문헌 예시는 논문 본문에서 확인한 서열이고, 내 서열은 검색·채택·직접 입력한 서열 중 이 표면에 결합한다고 기록한 것입니다. "
                   "참고 서열이 많을수록 그 서열들의 잔기 구성과 연속 모티프를 닮은 후보가 우선됩니다.")
        with st.expander("참고 서열 목록과 각 서열의 표면 점수"):
            st.dataframe(pd.DataFrame([{
                "구분": k.split(" · ")[0], "이름": s.name, "서열": s.sequence, "길이": len(s.sequence),
                "표면 점수": round(score_sequence(s.sequence, class_key).scaled, 1), "출처": s.source,
            } for k, s in options.items() if k in chosen_refs]), hide_index=True, width="stretch")

    with st.expander("고급 설정"):
        weight = st.slider("알려진 서열과 닮은 정도의 반영 비중", 0.0, 0.8, 0.4, step=0.05, key="rec_sim",
                           help="높일수록 알려진 서열과 비슷한 잔기 구성·모티프가 우선됩니다. 0이면 표면 가중치 점수만 봅니다.")
        kmer = st.slider("유사도에서 연속 모티프(2-mer)가 차지하는 비중", 0.0, 1.0, 0.4, step=0.1, key="rec_kmer",
                         help="0이면 잔기 구성만, 1이면 알려진 서열에 들어 있는 연속한 두 잔기 조합만 봅니다.")
        c = CLASS_CONSTRAINTS[class_key]
        st.caption(f"적용되는 제약: 한 잔기 최대 {c.max_residue_frac:.0%}, 같은 잔기 {c.max_run}연속까지, Cys {c.max_cys}개까지, "
                   f"전하 잔기(D·E·K·R) {c.min_charged}개 이상, 순전하 {c.charge_range[0]:+g}~{c.charge_range[1]:+g}, "
                   f"GRAVY {c.gravy_range[0]:g}~{c.gravy_range[1]:g}")

    rounds = st.session_state.get("rec_round", 0)
    b1, b2 = st.columns([1, 1])
    run = b1.button("서열 추천 받기", type="primary", key="rec_go")
    again = b2.button("다른 후보 더 보기", key="rec_again", help="같은 설정으로 다시 탐색해 다른 후보를 만듭니다.")
    if run or again:
        rounds = rounds + 1 if again else 0
        st.session_state["rec_round"] = rounds
        seeds = [options[k] for k in chosen_refs]
        with st.spinner("후보를 탐색하는 중…"):
            result = recommend(class_key, n=n, min_len=length[0], max_len=length[1], seeds=seeds, use_literature_seeds=False,
                               seed_value=rounds, similarity_weight=weight, kmer_weight=kmer)
        st.session_state["rec_result"] = {"result": result, "label": cls.label}
        log_run("rec", st.session_state["rec_result"])
    state = st.session_state.get("rec_result")
    if state:
        _results(state)


def _results(state: dict) -> None:
    result = state["result"]
    cls = CLASSES[result.class_key]
    if not result.candidates:
        st.info("조건을 만족하는 후보를 찾지 못했습니다. 길이 범위를 넓히거나 닮은 정도의 비중을 낮춰 보세요.")
        return
    st.markdown(f"#### 추천 서열 {len(result.candidates)}개 — {cls.label.split(' (')[0]}")
    if result.seeds:
        st.caption(f"알려진 서열 {len(result.seeds)}개를 참고했습니다.")
    if result.baseline:
        base_rows = [{"알려진 서열": s.name, "서열": s.sequence, "표면 점수": result.baseline.get(s.name), "출처": s.source}
                     for s in result.seeds]
        with st.expander(f"기준선: 알려진 서열 {len(base_rows)}개 자신의 점수 (추천 점수와 비교해 보세요)", expanded=result.weak_baseline):
            st.dataframe(pd.DataFrame(base_rows), hide_index=True, width="stretch")
        if result.weak_baseline:
            st.warning("이 표면에서는 **문헌에서 실제로 결합한다고 알려진 서열의 점수가 낮게** 나옵니다. 이 점수가 이 표면의 실제 결합 서열 특징을 "
                       "충분히 잡지 못한다는 뜻이니, 추천 서열의 점수를 곧이곧대로 믿지 말고 '알려진 서열과 닮은 정도'를 더 중요하게 보세요.")

    rows = []
    for i, c in enumerate(result.candidates):
        rows.append({
            "순위": i + 1, "서열": c.sequence, "길이": c.length,
            "예상 결합 성향": c.score, "표면 점수": c.binding,
            "닮은 알려진 서열": f"{c.nearest} ({c.similarity:.2f})" if c.similarity is not None else "-",
            "신규성": c.novelty, "순전하": c.charge, "GRAVY": c.gravy, "분자량": c.mw,
            "비고": "; ".join(c.notes),
        })
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch", height=min(420, 36 * len(rows) + 40),
                 column_config={"예상 결합 성향": st.column_config.ProgressColumn(
                     "예상 결합 성향", min_value=0, max_value=100, format="%.1f",
                     help="표면 가중치 점수와 알려진 서열과의 유사도를 합친 학습용 점수에서 제약 위반 벌점을 뺀 값 (0~100)")})
    st.caption("신규성: 1이면 알려진 서열과 전혀 겹치지 않음. 추천 후보끼리는 잔기 구성이 서로 3개 이상 다르게 골랐습니다. "
               f"탐색한 서열 {result.evaluations:,}개 · 같은 설정이면 결과가 재현됩니다.")

    excel_download_button("rec_xlsx")
    labels = {i: f"#{i + 1} {c.sequence} (예상 {c.score})" for i, c in enumerate(result.candidates)}
    pick = st.selectbox("잔기별 기여 보기", list(labels), format_func=labels.get, key="rec_pick")
    c: Recommendation = result.candidates[pick]
    st.pyplot(plots.residue_contribution_plot(c.sequence, c.contributions, f"{c.sequence} — {cls.label}"))

    st.markdown("##### 물성 분석 (저장하지 않아도 바로 볼 수 있습니다)")
    items = [(f"#{i + 1}", cand.sequence) for i, cand in enumerate(result.candidates)]
    with st.expander("물성 비교 (전체 후보)", expanded=True):
        st.dataframe(properties_table(items, label_column="순위"), hide_index=True, width="stretch")
        st.caption("불안정성 지수 40 초과는 시험관 내 불안정 가능성을 뜻하는 경험적 지표입니다. 정렬·소수성 윈도우·helical wheel은 내 서열에 추가한 뒤 📊 서열 분석 탭에서 볼 수 있습니다.")

    chosen = st.multiselect("내 서열에 추가할 후보 (여러 개 선택 가능)", list(labels), format_func=labels.get, key="rec_add_pick",
                            placeholder="후보를 선택하세요")
    if st.button(f"선택한 {len(chosen)}개를 내 서열에 추가", disabled=not chosen, key="rec_add"):
        first = None
        for i in chosen:
            cand = result.candidates[i]
            rec = SeqRecord(
                name=unique_name(f"{cand.sequence} · {cls.label.split(' (')[0]} 추천"), sequence=cand.sequence,
                provenance=Provenance(SOURCE_DESIGNED, f"자동 추천 · {cls.label.split(' (')[0]} · 예상 {cand.score}"),
                is_designed=True, mutation="자동 추천",
                reference_info={"recommended_for": cls.label, "recommend_score": cand.score},
            )
            add_record(rec, make_current=False)
            first = first or rec.name
        set_current(first)
        flash(f"{len(chosen)}개를 설계된 서열(비실험)로 추가했습니다. 서열 분석·결합 친화도 탭에서 이어서 확인하세요.")
        st.rerun()

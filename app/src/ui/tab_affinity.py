"""결합 친화도 탭: ① 표면 결합 성향 점수(학습용 예측) ② 문헌에서 보고된 친화도 기록 ③ 예측 vs 보고값 비교."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from ..affinity import (
    AFFINITY_METHODS, AFFINITY_METRICS, CALCULATED_METHOD, METRIC_COMPUTED, METRIC_KD, STATUS_CONFIRMED,
    STATUS_SUGGESTED, UNITS, Affinity, extract_affinity_candidates, rank_correlation,
)
from ..analysis import plots
from ..analysis.binding_score import CLASSES, class_for_material, score_sequence
from ..provenance import BindingTarget, Provenance, SeqRecord, SOURCE_PAPER
from ..sequence import STANDARD_AA
from .common import excel_download_button
from .state import add_record, current_record, flash, records, unique_name

CAVEAT = ("**학습용 경험적 점수입니다.** 결합 친화도(Kd)의 예측이 아니라, 표면 화학별 *정성적* 경향(예: 소수성 고분자 ← 방향족 잔기, "
          "티타늄·실리카 ← 양전하 잔기)을 단순화해 서열끼리 줄 세우는 용도입니다. pH·이온 세기·표면 상태·펩타이드 입체 구조는 "
          "반영하지 않으며, 아래 ②③에서 문헌 값과 비교해 얼마나 맞는지 직접 확인하세요.")

# 문헌 예시: PMID 37735840 (플라스틱 결합 펩타이드의 계산 설계) 표에 실린 상위 설계 서열과 PepBD 계산 점수
LITERATURE_EXAMPLES = [
    ("PET (폴리에틸렌 테레프탈레이트)", "WEWWHVFHHRLR", -63.8),
    ("폴리에틸렌", "HWMWAMKWHMRH", -59.6),
    ("폴리프로필렌", "WWQRHMFHFRTW", -52.3),
    ("폴리스티렌", "WHWQRIIWQQMR", -44.5),
]
EXAMPLE_PMID = "PMID 37735840"


def render() -> None:
    recs = records()
    if not recs:
        st.info("먼저 '검색' 탭에서 서열을 추가하세요. 문헌 예시는 이 탭의 ②에서 한 번에 불러올 수 있습니다.")
        _examples_button(show_only=True)
        return
    st.warning(CAVEAT)
    excel_download_button("aff_xlsx")
    t1, t2, t3 = st.tabs(["① 결합 성향 점수 (예측)", "② 보고된 친화도 기록 (문헌)", "③ 예측 vs 보고값 비교"])
    with t1:
        _scores(recs)
    with t2:
        _records(recs)
    with t3:
        _compare(recs)


# ---------------------------------------------------------------------- ① 예측
def _class_select(default_key: str, key: str, label: str = "표면 분류") -> str:
    keys = list(CLASSES)
    return st.selectbox(label, keys, index=keys.index(default_key) if default_key in keys else 0,
                        format_func=lambda k: CLASSES[k].label, key=key)


def _scores(recs: dict) -> None:
    cur = current_record()
    guess = class_for_material(cur.binding_target.name) if cur else None
    cls_key = _class_select(guess or "polymer", key=f"aff_class_{cur.name if cur else ''}")
    cls = CLASSES[cls_key]
    if guess is None and cur is not None:
        st.caption("현재 서열의 결합 대상에서 표면 분류를 추정하지 못해 기본값을 썼습니다. 직접 골라 주세요.")
    st.caption(f"가정: {cls.rationale}")
    with st.expander("잔기별 가중치 (가정) 보기"):
        weights = pd.DataFrame({"잔기": list(STANDARD_AA), "가중치": [cls.weights.get(a, 0.0) for a in STANDARD_AA]})
        st.dataframe(weights.sort_values("가중치", ascending=False), hide_index=True, width="stretch", height=300)

    names = list(recs)
    default = [cur.name] if cur else names[:1]
    default += [n for n, r in recs.items() if cur and r.parent == cur.name][:4]
    chosen = st.multiselect("점수를 볼 서열", names, default=default, key="aff_pick")
    if not chosen:
        return
    rows = []
    scores = {}
    for n in chosen:
        s = score_sequence(recs[n].sequence, cls_key)
        scores[n] = s
        top = sorted(range(len(s.contributions)), key=lambda i: -s.contributions[i])[:3]
        rows.append({
            "이름": n, "길이": len(recs[n].sequence), "결합 성향 (0–100)": round(s.percentile, 1),
            "가중치 척도 (0–100)": round(s.scaled, 1), "z": round(s.z, 2),
            "잔기당 평균 가중치": round(s.raw, 3),
            "가장 기여한 잔기": ", ".join(f"{recs[n].sequence[i]}{i + 1}" for i in top if s.contributions[i] > 0),
        })
    st.dataframe(pd.DataFrame(rows).sort_values("결합 성향 (0–100)", ascending=False), hide_index=True, width="stretch",
                 column_config={"결합 성향 (0–100)": st.column_config.ProgressColumn(
                     "결합 성향 (0–100)", min_value=0, max_value=100, format="%.1f",
                     help="같은 길이의 무작위 서열보다 점수가 높을 확률(백분위). 50이면 무작위 수준입니다.")})
    st.caption("결합 성향은 절대 친화도가 아니라 **무작위 서열 대비 백분위**입니다(강한 후보끼리는 100 근처로 몰립니다). "
               "그럴 때는 가중치 척도(이 표면의 최소~최대 가중치 구간에 대응)로 후보를 구분하세요. 짧은 서열일수록 우연히 높게 나오기 쉽습니다.")
    pick = st.selectbox("잔기별 기여 보기", chosen, key="aff_residue_pick")
    s = scores[pick]
    st.pyplot(plots.residue_contribution_plot(recs[pick].sequence, s.contributions, f"{pick} — {cls.label}"))


# --------------------------------------------------------------- ② 보고된 친화도
def _examples_button(show_only: bool = False) -> None:
    st.markdown("**문헌 예시로 연습하기**")
    st.caption(f"{EXAMPLE_PMID}(*In Silico Design and Analysis of Plastic-Binding Peptides*)의 표에 실린 PET·PE·PP·PS 결합 설계 펩타이드 4개와 "
               "PepBD **계산** 점수를 불러옵니다. 실험 측정값이 아니라 계산값이며, 단위는 kcal/mol로 가정했으니 원문에서 확인하세요.")
    if st.button("문헌 예시 4개 추가", key="aff_examples"):
        added = 0
        existing = {r.sequence for r in records().values()}
        for target, seq, energy in LITERATURE_EXAMPLES:
            if seq in existing:
                continue
            rec = SeqRecord(
                name=unique_name(f"{seq} · {target.split(' ')[0]} 결합 설계"), sequence=seq,
                provenance=Provenance(SOURCE_PAPER, EXAMPLE_PMID, "https://pubmed.ncbi.nlm.nih.gov/37735840/",
                                      title="In Silico Design and Analysis of Plastic-Binding Peptides"),
                screening_methods=["인실리코 스크리닝"],
                binding_target=BindingTarget(target, "고분자", "계산 설계 대상", STATUS_CONFIRMED),
                affinities=[Affinity(METRIC_COMPUTED, energy, "kcal/mol", CALCULATED_METHOD, target,
                                     "PepBD 계산(논문 표의 PepBD 점수)", EXAMPLE_PMID, STATUS_CONFIRMED,
                                     "측정값이 아닌 계산 결합 에너지. 단위는 원문 확인 필요")],
            )
            add_record(rec, make_current=False)
            added += 1
        flash(f"문헌 예시 {added}개를 추가했습니다." if added else "이미 모두 추가돼 있습니다.", "success" if added else "info")
        st.rerun()


def _fill_from_candidate(cands: list, key: str, rec_name: str) -> None:
    """후보 값을 입력 폼에 채운다 (버튼 콜백에서 호출 — 위젯 생성 전에 상태를 바꿀 수 있다)."""
    idx = st.session_state.get(key)
    if idx is None or idx >= len(cands):
        return
    c = cands[idx]
    st.session_state[f"aff_f_metric_{rec_name}"] = c.metric
    st.session_state[f"aff_f_value_{rec_name}"] = f"{c.value:g}"
    st.session_state[f"aff_f_unit_{rec_name}_{c.metric}"] = c.unit
    if c.methods:
        st.session_state[f"aff_f_method_{rec_name}"] = c.methods[0]
    st.session_state[f"aff_f_status_{rec_name}"] = STATUS_SUGGESTED
    st.session_state[f"aff_f_note_{rec_name}"] = f"본문에서 찾은 후보: “{c.context[:120]}” (어느 펩타이드의 값인지 확인 필요)"


def _records(recs: dict) -> None:
    cur = current_record()
    names = list(recs)
    name = st.selectbox("친화도를 기록할 서열", names, index=names.index(cur.name) if cur else 0, key="aff_rec_pick")
    rec = recs[name]
    if rec.is_designed:
        st.caption("설계된 서열(비실험)입니다. 실제로 측정되지 않았다면 기록하지 마세요 (계산값은 '계산 결합 에너지'로 기록).")

    st.markdown(f"**{name}** 의 보고된 친화도")
    if rec.affinities:
        st.dataframe(pd.DataFrame([{
            "#": i + 1, "지표": a.metric, "값": a.value, "단위": a.unit, "측정 방법": a.method, "대상": a.target,
            "조건": a.conditions, "출처": a.source, "ΔG (kcal/mol)": None if a.delta_g is None else round(a.delta_g, 2),
            "상태": a.status + (" · 계산값" if a.computed else ""),
        } for i, a in enumerate(rec.affinities)]), hide_index=True, width="stretch")
        c1, c2 = st.columns([3, 1])
        which = c1.selectbox("삭제할 항목", range(len(rec.affinities)), format_func=lambda i: f"#{i + 1} {rec.affinities[i].label()}",
                             key=f"aff_del_{name}_{len(rec.affinities)}")
        if c2.button("삭제", key=f"aff_del_btn_{name}"):
            rec.affinities.pop(which)
            st.rerun()
    else:
        st.caption("아직 기록이 없습니다 — 친화도: 미보고. 아래에서 논문에 적힌 값을 추가하세요.")

    with st.expander("논문 본문에서 Kd·ΔG 값 후보 찾기 (붙여넣기)"):
        text = st.text_area("본문 일부를 붙여넣으세요", height=120, key=f"aff_text_{name}",
                            help="예: 'The K_D was 12 nM as measured by QCM-D …' 처럼 값이 있는 문단")
        cands = extract_affinity_candidates(text) if text.strip() else []
        if text.strip() and not cands:
            st.caption("Kd·Ka·ΔG 형태의 값을 찾지 못했습니다. 아래 폼에 직접 입력하세요.")
        if cands:
            st.dataframe(pd.DataFrame([{"#": i + 1, "값": c.label(), "측정 방법 후보": ", ".join(c.methods) or "-", "문맥": c.context}
                                       for i, c in enumerate(cands)]), hide_index=True, width="stretch")
            ckey = f"aff_cand_{name}"
            st.selectbox("폼에 채울 후보", range(len(cands)), format_func=lambda i: f"#{i + 1} {cands[i].label()}", key=ckey)
            st.button("선택한 후보로 아래 폼 채우기", key=f"aff_fill_{name}", on_click=_fill_from_candidate, args=(cands, ckey, name))
            st.caption("자동 추출은 어느 펩타이드·어느 대상의 값인지 알 수 없습니다. 채운 뒤 문맥을 읽고 확정하세요 ('제안' 상태로 들어갑니다).")

    st.markdown("**친화도 추가**")
    c1, c2, c3 = st.columns([2, 2, 2])
    metric = c1.selectbox("지표", AFFINITY_METRICS, key=f"aff_f_metric_{name}")
    value_text = c2.text_input("값 (예: 12, 3.2e-7)", key=f"aff_f_value_{name}")
    unit = c3.selectbox("단위", UNITS[metric], key=f"aff_f_unit_{name}_{metric}")
    default_method = CALCULATED_METHOD if metric == METRIC_COMPUTED else None
    method_options = AFFINITY_METHODS + ["직접 입력"]
    method_index = method_options.index(default_method) if default_method in method_options else 0
    c1, c2 = st.columns(2)
    method_sel = c1.selectbox("측정 방법", method_options, index=method_index, key=f"aff_f_method_{name}")
    method = c2.text_input("측정 방법 이름", key=f"aff_f_method_custom_{name}") if method_sel == "직접 입력" else method_sel
    c1, c2 = st.columns(2)
    target = c1.text_input("결합 대상(표면)", value=rec.binding_target.name, key=f"aff_f_target_{name}")
    conditions = c2.text_input("측정 조건 (pH, 완충액, 이온 세기, 온도)", key=f"aff_f_cond_{name}")
    c1, c2 = st.columns(2)
    source = c1.text_input("출처 (PMID/DOI)", value=rec.provenance.identifier if rec.provenance.source_type == SOURCE_PAPER else "",
                           key=f"aff_f_source_{name}")
    status = c2.radio("상태", [STATUS_CONFIRMED, STATUS_SUGGESTED], horizontal=True, key=f"aff_f_status_{name}")
    note = st.text_input("메모", key=f"aff_f_note_{name}")

    value = None
    if value_text.strip():
        try:
            value = float(value_text.replace(",", ""))
        except ValueError:
            st.error("값은 숫자로 입력하세요 (예: 12, 3.2e-7).")
    if metric == METRIC_KD and value is not None and value <= 0:
        st.error("Kd는 0보다 커야 합니다.")
        value = None
    if value is not None:
        preview = Affinity(metric, value, unit, method, target, conditions, source, status, note)
        if preview.delta_g is not None:
            st.caption(f"ΔG ≈ {preview.delta_g:.2f} kcal/mol (298 K 환산, 더 음수일수록 강한 결합)")
    if st.button("친화도 추가", type="primary", key=f"aff_add_{name}", disabled=value is None):
        rec.affinities.append(Affinity(metric, value, unit, method, target.strip(), conditions.strip(), source.strip(), status, note.strip()))
        flash(f"'{name}'에 친화도 {metric} {value:g} {unit}를 기록했습니다.")
        st.rerun()

    st.divider()
    _examples_button()


# ------------------------------------------------------------------ ③ 비교
def _compare(recs: dict) -> None:
    fallback = _class_select("polymer", key="aff_cmp_class", label="결합 대상에서 표면 분류를 알 수 없을 때 사용할 분류")
    rows, points = [], []
    for name, rec in recs.items():
        usable = [a for a in rec.affinities if a.delta_g is not None]
        if not usable:
            continue
        measured = [a for a in usable if not a.computed]
        aff = (measured or usable)[0]
        cls_key = class_for_material(rec.binding_target.name) or class_for_material(aff.target) or fallback
        s = score_sequence(rec.sequence, cls_key)
        rows.append({
            "이름": name, "표면 분류": CLASSES[cls_key].label.split(" (")[0], "결합 성향(백분위)": round(s.percentile, 1),
            "결합 성향(가중치 척도)": round(s.scaled, 1),
            "보고된 값": aff.label(), "ΔG (kcal/mol)": round(aff.delta_g, 2), "종류": "계산값" if aff.computed else "측정값",
        })
        points.append({"name": name, "x": s.scaled, "dg": aff.delta_g, "computed": aff.computed})
    if not rows:
        st.info("비교하려면 ②에서 Kd·Ka·ΔG(또는 계산 결합 에너지)를 기록한 서열이 필요합니다. ②의 '문헌 예시 4개 추가'로 바로 체험할 수 있습니다.")
        return
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
    rho = rank_correlation([p["x"] for p in points], [-p["dg"] for p in points])
    st.pyplot(plots.affinity_scatter_plot(points, rho))
    if rho is None:
        st.caption("점이 3개 미만이거나 값이 모두 같아 상관을 계산하지 않았습니다.")
    else:
        verdict = "예측 순위와 보고값 순위가 같은 방향" if rho > 0.3 else "방향이 거의 맞지 않음" if rho > -0.3 else "예측과 반대 방향"
        st.metric("순위상관 ρ (점수 vs 결합 세기)", f"{rho:.2f}", help="+1이면 점수가 높을수록 더 강하게 결합한다고 보고된 경우")
        st.caption(f"→ {verdict}. 단, n={len(points)}개로는 통계적 의미가 거의 없습니다.")
    st.info("비교 시 주의: 서로 다른 표면·측정 방법·조건·계산/측정값이 섞여 있으면 직접 비교할 수 없습니다. "
            "이 비교는 점수가 어떤 경향을 놓치는지 *생각해 보는* 용도입니다.")

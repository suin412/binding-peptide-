"""서열 분석 탭: 물성 비교, 서열 정렬, 소수성·조성·helical wheel."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from ..analysis import plots
from ..analysis.align import align_sequences, render_html
from ..analysis.properties import compute_properties
from .common import excel_download_button, quick_add_form
from .state import current_record, records

PROP_COLUMNS = {
    "length": "길이(aa)", "mw": "분자량(Da)", "pI": "pI", "net_charge_pH7.4": "순전하(pH 7.4)",
    "gravy": "GRAVY", "instability_index": "불안정성 지수", "aromaticity": "방향족성",
}


def render() -> None:
    quick_add_form("analysis", expanded=not records())
    recs = records()
    if not recs:
        st.info("위에서 새 서열을 입력하거나, '검색'·'서열 추천' 탭에서 서열을 추가하세요.")
        return
    names = list(recs)
    cur = current_record()
    default = [cur.name] if cur else names[:1]
    children = [n for n, r in recs.items() if cur and r.parent == cur.name][:3]
    chosen = st.multiselect("분석할 서열 (여러 개 선택하면 비교)", names, default=default + children, key="an_pick")
    if not chosen:
        st.info("분석할 서열을 선택하세요.")
        return

    rows = []
    for n in chosen:
        p = compute_properties(recs[n].sequence)
        rows.append({"이름": n, "구분": "설계" if recs[n].is_designed else "원본", **{PROP_COLUMNS[k]: v for k, v in p.items()}})
    st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
    excel_download_button("an_xlsx")
    st.caption("순전하는 Biopython `charge_at_pH(7.4)`(측쇄·말단기 pKa 기반 추정)입니다. 불안정성 지수 40 초과는 시험관 내 불안정 가능성을 뜻하는 경험적 지표입니다.")

    if len(chosen) >= 2:
        st.markdown("##### 서열 정렬")
        c1, c2 = st.columns(2)
        a = c1.selectbox("기준(A)", chosen, key="an_align_a")
        b = c2.selectbox("비교(B)", [n for n in chosen if n != a] or chosen, key="an_align_b")
        res = align_sequences(recs[a].sequence, recs[b].sequence)
        m1, m2, m3 = st.columns(3)
        m1.metric("동일성", f"{res.identity * 100:.1f}%")
        m2.metric("일치 칸", f"{res.matches}/{len(res.midline)}")
        m3.metric("차이 위치(A 기준)", ", ".join(map(str, res.diff_positions_a[:8])) + ("…" if len(res.diff_positions_a) > 8 else "") or "없음")
        st.markdown(render_html(res, a, b), unsafe_allow_html=True)

    series = {n: recs[n].sequence for n in chosen}
    st.markdown("##### 시각화")
    c1, c2 = st.columns(2)
    window = c1.slider("소수성 윈도우", 3, 15, 7, step=2, key="an_window")
    c1.pyplot(plots.hydropathy_plot(series, window))
    c2.pyplot(plots.composition_plot(series))
    wheel_name = st.selectbox("Helical wheel (양친매성 헬릭스 확인)", chosen, key="an_wheel")
    wheel_seq = recs[wheel_name].sequence
    limit = st.slider("표시할 잔기 수", 5, max(5, min(40, len(wheel_seq))), min(18, max(5, len(wheel_seq))), key="an_wheel_n") if len(wheel_seq) > 5 else len(wheel_seq)
    st.pyplot(plots.helical_wheel(wheel_seq[:limit], f"{wheel_name} (1-{limit})"))

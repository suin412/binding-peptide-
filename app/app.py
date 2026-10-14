"""PepLab — 펩타이드/단백질 설계 연습 앱 (Streamlit)."""
from __future__ import annotations

import pandas as pd
import streamlit as st

from src.analysis.properties import compute_properties
from src.sequence import to_fasta
from src.ui import tab_affinity, tab_analysis, tab_improve, tab_recommend, tab_search
from src.ui.common import apply_context, context_editor, excel_download_button, render_header
from src.ui.state import current_record, delete_record, init_state, records, show_flash

st.set_page_config(page_title="PepLab — 펩타이드 설계 연습", page_icon="🧬", layout="wide")
init_state()


def sidebar():
    """사이드바를 그리고, 엑셀 버튼 자리를 돌려준다(버튼은 탭이 모두 그려진 뒤 채운다 — 그래야 방금 실행한 추천 결과가 바로 담긴다)."""
    st.sidebar.title("🧬 PepLab")
    st.sidebar.caption("펩타이드/단백질 검색 → 개선형·전체 서열 기반 추천 → 분석 연습용 앱")
    excel_slot = st.sidebar.container()
    recs = records()
    st.sidebar.subheader(f"내 서열 ({len(recs)})")
    if not recs:
        st.sidebar.info("아직 서열이 없습니다. '검색' 탭에서 프리셋이나 UniProt 항목을 추가하세요.")
        return excel_slot
    names = list(recs)
    cur = current_record()

    def label(n: str) -> str:
        r = recs[n]
        return ("🧪 " if r.is_designed else "📄 ") + n

    selected = st.sidebar.radio("현재 서열", names, index=names.index(cur.name) if cur else 0,
                                format_func=label, key="sidebar_pick")
    st.session_state["current"] = selected
    rec = recs[selected]
    if st.sidebar.button("🗑 현재 서열 삭제", key="sidebar_delete"):
        delete_record(selected)
        st.rerun()

    if not rec.is_designed:
        with st.sidebar.expander("발굴 방법 · 결합 대상 편집"):
            methods, note, target = context_editor(f"side_ctx_{selected}", initial=rec)
            if st.button("저장", key=f"side_ctx_save_{selected}"):
                apply_context(rec, methods, note, target)
                st.rerun()

    st.sidebar.divider()
    st.sidebar.download_button("FASTA 내보내기", to_fasta([(n, r.sequence) for n, r in recs.items()]),
                               file_name="peplab_sequences.fasta", key="dl_fasta")
    rows = []
    for n, r in recs.items():
        p = compute_properties(r.sequence)
        rows.append({
            "name": n, "sequence": r.sequence, "length": p["length"], "parent": r.parent or "", "mutation": r.mutation,
            "designed": r.is_designed, "mw": p["mw"], "pI": p["pI"], "charge_pH7.4": p["net_charge_pH7.4"], "gravy": p["gravy"],
            "source": r.provenance.short(), "binding_target": r.binding_target.label(), "screening_methods": r.methods_label(),
            "reported_affinity": r.affinity_label(),
        })
    st.sidebar.download_button("CSV 내보내기", pd.DataFrame(rows).to_csv(index=False).encode("utf-8-sig"),
                               file_name="peplab_sequences.csv", mime="text/csv", key="dl_csv")
    return excel_slot


excel_slot = sidebar()

# 탭 위쪽은 항상 같은 컨테이너 하나로 둔다: 서열을 추가하면 헤더·안내 메시지가 생기고 사라지는데, 이것이 탭과 같은 층에 있으면
# 탭의 위치가 밀려 화면이 첫 탭(검색)으로 되돌아가기 때문이다.
top = st.container()
with top:
    show_flash()
    st.title("펩타이드 · 단백질 설계 연습")
    cur = current_record()
    if cur is not None:
        render_header(cur)
        st.divider()

t_search, t_imp, t_analysis, t_aff, t_rec = st.tabs(["🔎 검색", "🔧 개선형 추천", "📊 서열 분석", "🧲 결합 친화도", "🎯 전체 서열 기반 추천"])
with t_search:
    tab_search.render()
with t_imp:
    tab_improve.render()
with t_analysis:
    tab_analysis.render()
with t_aff:
    tab_affinity.render()
with t_rec:
    tab_recommend.render()

if records() or st.session_state.get("rec_result") or st.session_state.get("imp_result") or st.session_state.get("lit_result"):      # 저장하지 않은 추천 결과만 있어도 내보낼 수 있다
    excel_download_button("dl_xlsx", target=excel_slot)

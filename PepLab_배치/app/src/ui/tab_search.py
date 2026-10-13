"""검색 탭: UniProt → 논문/특허 → 직접 입력·프리셋 (서열 출처 폴백 체인)."""
from __future__ import annotations

import hashlib
from collections import Counter
from datetime import date

import pandas as pd
import streamlit as st

from ..presets import SAMPLE_TEXTS, get_preset, preset_names, preset_record
from ..provenance import (
    BINDING_CATEGORIES, SCREENING_METHODS, STATUS_CONFIRMED, STATUS_SUGGESTED, BindingTarget,
    SOURCE_MANUAL, SOURCE_PAPER, SOURCE_PATENT, SOURCE_UNIPROT, Provenance, SeqRecord,
)
from ..sequence import SequenceError, parse_fasta, validate_sequence
from ..sources import binding, europepmc, literature, seqsearch, uniprot
from ..sources.http import ApiError
from ..sources.query_assist import expand_query
from ..sources.seqextract import Candidate, extract_candidates
from .common import apply_context, context_editor
from .state import add_record, current_record, flash, set_current, unique_name


def render() -> None:
    sub_db, sub_paper, sub_manual = st.tabs(["① UniProt 검색", "② 논문 · 특허에서 찾기", "③ 직접 입력 · 프리셋"])
    with sub_db:
        _uniprot_section()
    with sub_paper:
        _paper_patent_section()
    with sub_manual:
        _manual_section()


# ---------------------------------------------------------------- UniProt
def _uniprot_section() -> None:
    st.caption("ID(예: P01504), 키워드, 유전자명으로 검색합니다. 펩타이드는 UniProt에 없을 수 있으며, 그 경우 ②에서 논문·특허로 찾습니다.")
    c1, c2, c3 = st.columns([4, 2, 1])
    query = c1.text_input("검색어", placeholder="melittin, LL-37, P01504 …", key="up_query")
    reviewed = c2.checkbox("검토된(Swiss-Prot) 항목만", value=True, key="up_reviewed")
    go = c3.button("검색", key="up_go", width="stretch")
    if go and query.strip():
        try:
            with st.spinner("UniProt 검색 중…"):
                st.session_state["up_hits"] = uniprot.search(query, size=15, reviewed_only=reviewed)
                st.session_state["up_hits_query"] = query
        except ApiError as exc:
            st.error(str(exc))
            return

    hits = st.session_state.get("up_hits")
    if hits is None:
        return
    if not hits:
        st.info(f"'{st.session_state.get('up_hits_query', '')}'에 대한 UniProt 결과가 없습니다. "
                "②번 탭에서 논문·특허로 찾아보거나, '검토된 항목만'을 해제해 보세요.")
        return

    st.dataframe(
        pd.DataFrame([{"AC": h.accession, "ID": h.entry_id, "이름": h.name, "생물종": h.organism, "길이(aa)": h.length}
                      for h in hits]),
        hide_index=True, width="stretch",
    )
    labels = {h.accession: f"{h.accession} — {h.name} ({h.organism}, {h.length} aa)" for h in hits}
    acc = st.selectbox("자세히 볼 항목", list(labels), format_func=labels.get, key="up_selected")
    try:
        with st.spinner("항목 정보를 가져오는 중…"):
            entry = uniprot.fetch(acc)
    except ApiError as exc:
        st.error(str(exc))
        return
    _uniprot_detail(entry)


def _uniprot_detail(entry: uniprot.UniProtEntry) -> None:
    st.markdown(f"#### {entry.name}")
    st.caption(f"{entry.accession} · {entry.entry_id} · {entry.organism} · {entry.length} aa · [UniProt에서 열기]({entry.url})")
    if entry.function:
        with st.expander("기능 요약"):
            st.write(entry.function)
    pe_num = entry.protein_existence[:1]
    pe_text = uniprot.PE_LABELS.get(pe_num, entry.protein_existence)
    extra = " · 키워드 'Direct protein sequencing'" if entry.direct_sequencing else ""
    st.caption(f"참고 정보(발굴 방법이 아님): 단백질 존재 수준 {pe_text}{extra}")
    links = [f"[PDB {p}](https://www.rcsb.org/structure/{p})" for p in entry.pdb_ids[:6]]
    if links:
        st.markdown("실험 구조: " + " · ".join(links))
    refs = [r for r in entry.references if r["pmid"]]
    if refs:
        with st.expander(f"참고문헌 ({len(refs)}건)"):
            for r in refs[:15]:
                st.markdown(f"- [{r['title'] or 'PubMed'}](https://pubmed.ncbi.nlm.nih.gov/{r['pmid']}/) — {r['journal']} {r['year']} (PMID {r['pmid']})")

    region_labels = [f"전체 서열 (1-{entry.length}, {entry.length} aa)"] + [
        f"{r.label} ({r.start}-{r.end}, {r.end - r.start + 1} aa)" for r in entry.regions
    ]
    choice = st.selectbox("사용할 구간 — 성숙 펩타이드/체인을 고를 수 있습니다", region_labels, key=f"up_region_{entry.accession}")
    idx = region_labels.index(choice)
    if idx == 0:
        seq, start, end, full = entry.sequence, 1, entry.length, True
        default_name = f"{entry.name} ({entry.accession})"
    else:
        region = entry.regions[idx - 1]
        seq, start, end, full = entry.sequence[region.start - 1:region.end], region.start, region.end, False
        default_name = f"{region.label} ({entry.accession}:{region.start}-{region.end})"
    st.code(seq, language=None, wrap_lines=True)

    with st.expander("발굴 방법 · 결합 대상 기록 (선택)"):
        methods, note, target = context_editor(f"up_ctx_{entry.accession}")
    name = st.text_input("저장할 이름", value=default_name, key=f"up_name_{entry.accession}_{idx}")
    if st.button("이 서열 사용", type="primary", key=f"up_use_{entry.accession}_{idx}"):
        ident = entry.accession if full else f"{entry.accession}:{start}-{end}"
        rec = SeqRecord(
            name=unique_name(name), sequence=seq,
            provenance=Provenance(SOURCE_UNIPROT, ident, entry.url, title=entry.name),
            reference_info={
                "accession": entry.accession, "full_length": full, "organism": entry.organism,
                "protein_existence": pe_text, "direct_sequencing": entry.direct_sequencing, "pdb_ids": entry.pdb_ids,
            },
        )
        apply_context(rec, methods, note, target)
        add_record(rec)
        flash(f"'{rec.name}'을(를) 추가했습니다. 개선형 추천·서열 분석 탭에서 이어서 사용할 수 있습니다.")
        st.rerun()


# ------------------------------------------------------------ 논문 · 특허
MODE_DIRECT = "서열 바로 찾기 (키워드 → 여러 논문)"
MODE_PASTE = "텍스트 붙여넣기 (논문 본문 · 특허)"
MODE_TITLE = "논문 이름으로 검색 (선택)"
MODE_SEQ = "서열로 논문 찾기 (서열 → 논문)"


def _paper_patent_section() -> None:
    st.caption("UniProt에 없는 펩타이드는 논문 초록(Europe PMC)이나 특허 텍스트에서 서열 **후보**를 추출합니다. "
               "추출은 자동이지만 채택은 항상 직접 확인한 뒤에 합니다.")
    mode = st.radio("방식", [MODE_DIRECT, MODE_PASTE, MODE_TITLE, MODE_SEQ], horizontal=True, key="pp_mode")
    if mode == MODE_DIRECT:
        _literature_search()
    elif mode == MODE_PASTE:
        _paste_text()
    elif mode == MODE_SEQ:
        _sequence_search()
    else:
        _paper_search()


def _paper_provenance(paper) -> Provenance:
    """논문이면 '논문', Europe PMC의 특허 항목이면 '특허'로 출처 유형을 기록한다."""
    return Provenance(SOURCE_PATENT if paper.is_patent else SOURCE_PAPER, paper.identifier, paper.url, title=paper.title)


def _lit_sig(*parts) -> str:
    return hashlib.md5("|".join(map(str, parts)).encode()).hexdigest()[:8]


def _literature_search() -> None:
    """키워드 하나로 여러 논문(전 기간)을 훑어 서열 후보를 서열 단위로 바로 보여준다."""
    st.markdown("키워드를 넣으면 **여러 논문의 초록과 본문**을 읽어서 서열을 바로 보여줍니다. 논문을 하나씩 고를 필요가 없습니다. "
                "본문(전문)까지 읽으므로 시간이 걸립니다.")
    query = st.text_input("키워드 (영어 권장)", placeholder="titanium binding peptide", key="lit_query",
                          help="예: gold binding peptide, PET binding peptide, silica binding peptide phage display. "
                               "따옴표·괄호·AND/OR/NOT를 직접 쓰면 Europe PMC 검색식으로 그대로 사용합니다.")
    preview = expand_query(query) if query.strip() else None
    if preview and preview.changed:
        st.caption("🔎 " + " · ".join(preview.notes) + " — 실제 검색식은 결과 아래에서 확인할 수 있습니다.")
    this_year = date.today().year
    with st.expander("검색 범위 설정 (기본: 1975년~현재, 전 기간)"):
        years = st.slider("출판 연도", europepmc.EARLIEST_YEAR, this_year, (europepmc.EARLIEST_YEAR, this_year), key="lit_years")
        sort = st.selectbox("논문을 고르는 방식", europepmc.SORT_MODES, key="lit_sort",
                            help="Europe PMC 기본 정렬은 최근 논문 쪽으로 치우쳐 있어, 기본값은 연도 구간별로 고르게 가져옵니다.")
        max_papers = st.select_slider("확인할 논문 수 — 전체 일치가 이 수 이하이면 전부 읽습니다 (많을수록 오래 걸립니다)",
                                    [100, 500, 1000, 3000, 5000], value=3000, key="lit_max")
        min_len = st.number_input("서열 최소 길이 (짧을수록 오탐이 늘어납니다)", 5, 30, literature.DEFAULT_MIN_LEN, key="lit_minlen")
        expand = st.checkbox("재료 약어 자동 확장 (예: PET → polyethylene terephthalate, 영상 PET 논문 제외)", value=True, key="lit_expand",
                             help="PET·PS·PE·HA 같은 약어는 의학·생물 분야의 다른 뜻과 겹쳐 엉뚱한 논문이 많이 나옵니다. 끄면 입력한 그대로 검색합니다.")
        use_ft = st.checkbox("초록만이 아니라 **논문 본문(전문)까지 읽기** — 서열은 대부분 본문에 있습니다 (시간이 걸립니다)",
                             value=True, key="lit_ft",
                             help="PMC에 본문이 공개된 논문(오픈액세스 + 저자 원고)만 읽을 수 있습니다. Europe PMC가 막혀 있으면 NCBI PMC로 다시 시도합니다.")
        max_ft = st.number_input("본문을 읽을 최대 논문 수 (제목·초록에 peptide가 나오는 논문부터)", 5, 500, 100, key="lit_maxft") if use_ft else 100
        if use_ft:
            est = max(1, round(int(max_ft) * 3 / 3 / 60))
            st.caption(f"예상 소요: 본문 {int(max_ft)}편 기준 약 {est}~{est * 2 + 1}분 (3편씩 동시에 읽음). 같은 검색을 다시 하면 검색 결과는 캐시되지만 본문은 저작권상 저장하지 않아 다시 읽습니다.")

    if st.button("서열 찾기", type="primary", key="lit_go", disabled=not query.strip()):
        bar = st.progress(0.0, text="시작…")
        try:
            result = literature.find_sequences(
                query, max_papers=int(max_papers), year_from=years[0], year_to=years[1], sort=sort,
                use_fulltext=use_ft, max_fulltext=int(max_ft), min_len=int(min_len), expand_aliases=expand,
                progress=lambda frac, msg: bar.progress(frac, text=msg),
            )
        except ApiError as exc:
            bar.empty()
            st.error(str(exc))
            return
        bar.empty()
        st.session_state["lit_result"] = {
            "sig": _lit_sig(query, years, sort, max_papers, min_len, use_ft, max_ft, expand, len(result.papers)),
            "query": query, "result": result, "fulltext": use_ft,
        }

    state = st.session_state.get("lit_result")
    if state:
        _literature_results(state)


def _literature_results(state: dict) -> None:
    result: literature.LiteratureResult = state["result"]
    sig = state["sig"]
    if not result.papers:
        st.info(f"'{state['query']}'에 해당하는 논문이 없습니다. 영어 키워드로 바꾸거나 연도 범위를 넓혀 보세요.")
        return

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("확인한 논문", f"{len(result.papers):,}편", help=f"Europe PMC 전체 일치 {result.hit_count:,}편 중 선택 방식에 따라 가져온 수")
    c2.metric("서열 후보 (고유)", f"{len(result.groups):,}개")
    c3.metric("후보 출현 (중복 포함)", f"{len(result.hits):,}건")
    c4.metric("본문을 읽은 논문", f"{result.fulltext_read}편" if state["fulltext"] else "사용 안 함",
              help="PMC에 본문이 있어 실제로 전문을 읽은 논문 수")
    st.caption(f"Europe PMC 전체 일치 {result.hit_count:,}편 중 {len(result.papers):,}편을 확인했습니다. "
               "서열이 초록에 없거나 비오픈액세스 논문은 놓칠 수 있습니다.")
    if result.expansion_notes:
        st.info("검색어를 보정했습니다: " + " · ".join(result.expansion_notes))
    if result.effective_query:
        with st.expander("실제로 보낸 Europe PMC 검색식"):
            st.code(result.effective_query, language=None, wrap_lines=True)
    if not state["fulltext"] and len(result.groups) < 10:
        st.info("초록만으로는 서열이 적게 나올 수 있습니다. 서열은 보통 본문(방법·결과)에 있으니 **검색 범위 설정 ▸ 오픈액세스 전문까지 확인**을 켜고 다시 찾아보세요.")
    if state["fulltext"]:
        st.caption(f"본문을 읽은 논문 {result.fulltext_read}편 · 본문을 읽지 못해 초록만 확인한 논문 {len(result.unread)}편 "
                   f"(PMC에 본문이 없거나 접근할 수 없음{f', 전문 요청 실패 {result.fulltext_failed}편 포함' if result.fulltext_failed else ''}).")
        if result.unread:
            with st.expander(f"본문을 읽지 못한 논문 {len(result.unread)}편 — 직접 열어서 서열을 확인할 수 있습니다"):
                for p in result.unread[:200]:
                    st.markdown(f"- {p.year} · {p.title[:100]} — " + (f"[원문 링크]({p.url})" if p.url else p.identifier or "링크 없음"))

    with st.expander("연도별 분포 — 전 기간이 포함됐는지 확인"):
        years_hit = Counter(int(h.paper.year) for h in {h.paper.key: h for h in result.hits}.values() if h.paper.year.isdigit())
        df_years = pd.DataFrame({"논문 수": pd.Series(result.year_counts), "서열이 나온 논문 수": pd.Series(years_hit)}).fillna(0).astype(int).sort_index()
        st.bar_chart(df_years)

    if not result.groups:
        st.info("서열 후보가 없습니다. 최소 길이를 줄이거나, 전문 확인을 켜거나, 키워드를 바꿔 보세요. "
                "서열이 이미지/표로만 있는 논문은 읽을 수 없습니다.")
        return

    n_unc = sum(1 for g in result.groups if g.uncertain)
    h1, h2 = st.columns(2)
    hide_unc = h1.checkbox(f"문맥 불확실 후보 숨기기 ({n_unc}개)", value=True, key=f"lit_hide_unc_{sig}")
    hide_warn = h2.checkbox("경고 있는 후보도 숨기기", key=f"lit_hide_{sig}")
    groups = [g for g in result.groups if not (hide_unc and g.uncertain) and not (hide_warn and g.has_warning)]
    if not groups:
        st.info("조건에 맞는 후보가 없습니다. 위 체크를 해제해 보세요.")
        return

    rows = []
    for i, g in enumerate(groups):
        best = g.best
        yrs = g.years
        rows.append({
            "#": i + 1, "서열": g.sequence, "길이": len(g.sequence), "관련도": round(g.score, 1), "논문 수": g.paper_count,
            "연도": (f"{yrs[0]}" if len(yrs) == 1 else f"{yrs[0]}–{yrs[-1]}") if yrs else "",
            "대표 논문": f"{best.paper.identifier} · {best.paper.title[:60]}",
            "결합 대상 후보": ", ".join(t.name for t in best.targets[:2]),
            "출처": "/".join(sorted({h.source for h in g.hits})),
            "경고": ("? 문맥 불확실" if g.uncertain else "⚠" if g.has_warning else ""),
        })
    table_key = f"lit_table_{sig}_{int(hide_unc)}{int(hide_warn)}"      # 필터가 바뀌면 선택도 초기화
    event = st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch", height=min(420, 36 * len(rows) + 40),
                         on_select="rerun", selection_mode="multi-row", key=table_key)
    clicked = list(event.selection.rows) if event else []
    labels = {i: f"#{i + 1} · {g.sequence[:26]}{'…' if len(g.sequence) > 26 else ''} · 논문 {g.paper_count}편" for i, g in enumerate(groups)}
    flags = f"{int(hide_unc)}{int(hide_warn)}"
    picked = st.multiselect("채택할 서열 — 여러 개를 고를 수 있고, 표의 행을 클릭(여러 행 가능)해도 됩니다", list(labels),
                            format_func=labels.get, key=f"lit_multi_{sig}_{flags}", placeholder="서열을 선택하세요")
    all_on = st.checkbox(f"표시된 후보 {len(groups)}개 전부 선택 (한 번에 최대 {MAX_BULK}개)", key=f"lit_all_{sig}_{flags}")
    selected = list(range(len(groups))) if all_on else sorted(set(clicked) | set(picked))
    if not selected:
        st.caption("표에서 행을 클릭하거나 위 목록에서 서열을 고르세요. 하나를 고르면 자세히 확인·수정하며 채택하고, 여러 개를 고르면 한꺼번에 채택합니다.")
        return
    if len(selected) > 1:
        items = []
        for gi in selected[:MAX_BULK]:
            g = groups[gi]
            h = max(g.hits, key=lambda x: (not x.candidate.uncertain, x.score))
            items.append((h.candidate, _paper_provenance(h.paper),
                          h.targets, h.paper.identifier))
        if len(selected) > MAX_BULK:
            st.warning(f"한 번에 {MAX_BULK}개까지만 채택합니다. 앞의 {MAX_BULK}개만 처리됩니다.")
        _bulk_adopt_form(items, key=f"lit_bulk_{sig}_{flags}")
        return
    gi = selected[0]
    group = groups[gi]
    if len(group.hits) > 1:
        plabels = {i: f"{h.paper.year} · {h.paper.identifier} · {h.paper.title[:70]} ({h.source})" for i, h in enumerate(group.hits)}
        hi = st.selectbox("이 서열이 나온 논문 (출처로 기록됩니다)", list(plabels), format_func=plabels.get, key=f"lit_paper_{sig}_{gi}")
    else:
        hi = 0
    hit = group.hits[hi]
    paper = hit.paper
    st.markdown(f"**{paper.title}**  \n{paper.authors[:150]} — *{paper.journal}* {paper.year}"
                + (f" · [원문 링크]({paper.url})" if paper.url else ""))
    st.caption(f"추출 위치: {hit.source} · 문맥: “{hit.candidate.context}”")
    if hit.source == "전문":
        st.caption("전문에서 추출한 후보입니다. 참고문헌 영역은 제외했지만 표·방법 문장의 오탐 가능성이 있으니 문맥을 꼭 확인하세요.")
    prov = _paper_provenance(paper)
    _adopt_form(hit.candidate, prov, hit.targets, key=f"lit_{sig}_{gi}_{hi}", name_hint=paper.identifier)


def _sequence_search() -> None:
    """서열을 넣으면 그 서열이 언급된 논문과 위치(초록/본문)·문맥을 찾아 준다."""
    st.markdown("서열을 넣으면 **그 서열이 언급된 논문**을 찾고, 초록·본문(전문)에서 **어느 문맥에 나오는지** 확인합니다. "
                "논문 속 결합 대상·발굴 방법을 알아낼 때 쓰세요.")
    cur = current_record()
    seq_in = st.text_area("서열 (1문자 표기)", value=cur.sequence if cur else "", height=70,
                          key=f"sq_seq_{cur.name if cur else ''}",
                          help="현재 선택한 서열이 기본으로 들어갑니다. 정확히 같은 서열이 나온 논문을 찾습니다.")
    with st.expander("검색 옵션"):
        max_papers = st.select_slider("확인할 논문 수 (Europe PMC 관련도 순)", [20, 50, 100, 200], value=50, key="sq_max")
        contained = st.checkbox("더 긴 서열의 일부로 들어 있는 경우도 찾기 (느릴 수 있음, 짧은 서열은 잡음이 많음)", key="sq_contained")
        three = st.checkbox("3문자 표기(Gly-Ile-Gly …)로 적힌 경우도 찾기 (잡음이 많음)", key="sq_three")
        verify = st.checkbox("본문(전문)을 읽어 위치·문맥 확인 (PMC에 공개된 논문만, 시간이 걸림)", value=True, key="sq_verify")
    clean = None
    if seq_in.strip():
        try:
            clean = seqsearch.validate(seq_in)
        except SequenceError as exc:
            st.error(str(exc))
    if st.button("논문 찾기", type="primary", key="sq_go", disabled=clean is None):
        bar = st.progress(0.0, text="시작…")
        try:
            result = seqsearch.find_papers_for_sequence(
                clean, max_papers=int(max_papers), contained=contained, three_letter=three, verify_fulltext=verify,
                progress=lambda frac, msg: bar.progress(frac, text=msg),
            )
        except ApiError as exc:
            bar.empty()
            st.error(str(exc))
            return
        bar.empty()
        st.session_state["sq_result"] = {"result": result, "sig": _lit_sig(clean, max_papers, contained, three, verify)}
    state = st.session_state.get("sq_result")
    if state:
        _sequence_results(state)


def _sequence_results(state: dict) -> None:
    result: seqsearch.SeqSearchResult = state["result"]
    sig = state["sig"]
    located = [m for m in result.matches if m.located]
    c1, c2, c3 = st.columns(3)
    c1.metric("Europe PMC 일치", f"{result.hit_count:,}편")
    c2.metric("서열 위치 확인", f"{len(located)}편", help="초록 또는 본문에서 서열이 실제로 나오는 문맥을 찾은 논문")
    c3.metric("본문을 읽은 논문", f"{result.fulltext_read}편")
    for w in result.warnings:
        st.caption("⚠ " + w)
    with st.expander("실제로 보낸 Europe PMC 검색식"):
        st.code(" ".join(result.queries), language=None, wrap_lines=True)
    if not result.matches:
        st.info("이 서열이 들어 있는 논문을 찾지 못했습니다. 서열이 본문 표·이미지에만 있거나 PMC에 공개되지 않았을 수 있습니다. "
                "'더 긴 서열의 일부' 또는 '3문자 표기' 옵션을 켜 보세요.")
        return
    if result.fulltext_failed:
        st.caption(f"본문을 가져오지 못한 논문 {result.fulltext_failed}편은 초록만 확인했습니다.")

    st.dataframe(pd.DataFrame([{
        "#": i + 1, "연도": m.paper.year, "논문": f"{m.paper.identifier} · {m.paper.title[:70]}", "일치 방식": m.best_kind,
        "위치": m.where,
        "문맥": (m.occurrences[0].context if m.occurrences else "Europe PMC 색인에서 일치 - 본문 위치는 확인하지 못함"),
    } for i, m in enumerate(result.matches)]), hide_index=True, width="stretch", height=min(420, 36 * len(result.matches) + 40))
    labels = {i: f"#{i + 1} · {m.paper.year} · {m.paper.identifier} · {m.paper.title[:60]} ({m.best_kind})"
              for i, m in enumerate(result.matches)}
    idx = st.selectbox("자세히 볼 논문", list(labels), format_func=labels.get, key=f"sq_pick_{sig}")
    match = result.matches[idx]
    paper = match.paper
    st.markdown(f"**{paper.title}**  \n{paper.authors[:150]} — *{paper.journal}* {paper.year}"
                + (f" · [원문 링크]({paper.url})" if paper.url else ""))
    if not match.located:
        st.warning("Europe PMC 검색에서는 일치했지만 읽을 수 있는 초록·본문에서 서열 위치를 찾지 못했습니다 "
                   "(표·그림·보충자료 안에만 있거나 본문이 공개되지 않은 경우). 원문 링크에서 직접 확인하세요.")
    for occ in match.occurrences:
        st.markdown(f"- **{occ.where}** · {occ.kind}: “{occ.context}”")
    context_text = " ".join(o.context for o in match.occurrences)
    suggestions = binding.suggest_targets(f"{paper.title}. {paper.abstract} {context_text}", plastic_context=True)
    if suggestions:
        st.caption("이 논문에서 제안된 결합 대상 후보: " + ", ".join(s.name for s in suggestions[:4]))
    cand = Candidate(result.sequence, "1문자", context=match.occurrences[0].context if match.occurrences else "")
    prov = _paper_provenance(paper)
    with st.expander("이 서열을 이 논문을 출처로 채택"):
        _adopt_form(cand, prov, suggestions, key=f"sq_{sig}_{idx}", name_hint=paper.identifier)


def _paper_search() -> None:
    c1, c2 = st.columns([5, 1])
    query = c1.text_input("논문 검색어", placeholder="titanium binding peptide phage display", key="pp_query")
    go = c2.button("검색", key="pp_go", width="stretch")
    if go and query.strip():
        try:
            with st.spinner("Europe PMC 검색 중…"):
                st.session_state["pp_papers"] = europepmc.search(query, page_size=10)
        except ApiError as exc:
            st.error(str(exc))
            return
    papers = st.session_state.get("pp_papers")
    if papers is None:
        return
    if not papers:
        st.info("검색 결과가 없습니다. 영어 키워드로 다시 검색해 보세요.")
        return
    labels = {i: f"{p.year} · {p.title[:90]} ({p.identifier or 'ID 없음'})" for i, p in enumerate(papers)}
    idx = st.selectbox("논문 선택", list(labels), format_func=labels.get, key="pp_paper_idx")
    paper = papers[idx]
    st.markdown(f"**{paper.title}**  \n{paper.authors[:150]} — *{paper.journal}* {paper.year}"
                + (f" · [원문 링크]({paper.url})" if paper.url else ""))
    if not paper.abstract:
        st.warning("이 논문은 초록이 제공되지 않습니다. 원문에서 서열을 찾아 ③번 탭이나 '텍스트 붙여넣기'로 입력해 주세요.")
        return
    with st.expander("초록", expanded=False):
        st.write(paper.abstract)
    text = f"{paper.title}. {paper.abstract}"
    prov = _paper_provenance(paper)
    _candidates_ui(text, prov, key=f"paper_{idx}_{paper.identifier}", default_name_hint=paper.identifier)


def _paste_text() -> None:
    sample = st.selectbox("연습용 가상 샘플 (실제 논문/특허 아님)", ["(직접 붙여넣기)"] + list(SAMPLE_TEXTS), key="pp_sample")
    if sample != "(직접 붙여넣기)":
        if st.button("샘플 텍스트 불러오기", key="pp_load_sample"):
            st.session_state["pp_text"] = SAMPLE_TEXTS[sample]
            st.session_state["pp_src_type"] = SOURCE_PATENT if "특허" in sample else SOURCE_PAPER
            st.rerun()
    c1, c2 = st.columns([1, 2])
    src_type = c1.radio("출처 유형", [SOURCE_PAPER, SOURCE_PATENT], horizontal=True, key="pp_src_type")
    ident = c2.text_input("식별자 (PMID, DOI, 특허번호 등)", key="pp_ident",
                          placeholder="예: PMID 12345678 / US 1234567 B2")
    url = st.text_input("링크 (선택)", key="pp_url")
    text = st.text_area("본문 / 초록 / 서열목록 붙여넣기", height=200, key="pp_text",
                        help="특허의 서열목록(<400> 블록)이나 논문 본문을 붙여넣으세요. 저장되는 것은 서열과 출처 식별자뿐입니다.")
    if src_type == SOURCE_PATENT:
        st.caption("특허 검색은 [Google Patents](https://patents.google.com/), [Lens](https://www.lens.org/), "
                   "[WIPO PATENTSCOPE](https://patentscope.wipo.int/)에서 직접 찾은 뒤 본문을 붙여넣으세요 (사이트 약관상 자동 수집은 하지 않습니다).")
    if not text.strip():
        return
    prov = Provenance(src_type, ident.strip(), url.strip())
    _candidates_ui(text, prov, key="paste", default_name_hint=ident.strip())


def _candidates_ui(text: str, prov: Provenance, *, key: str, default_name_hint: str = "") -> None:
    candidates = extract_candidates(text)
    suggestions = binding.suggest_targets(text)

    if not candidates:
        st.info("서열 후보를 찾지 못했습니다. 서열이 이미지/표로 되어 있거나 표기가 특이할 수 있어요. ③번 탭에서 직접 입력해 주세요.")
    else:
        st.markdown(f"**서열 후보 {len(candidates)}개** — 자동 추출 결과이므로 오탐이 있을 수 있습니다.")
        st.dataframe(
            pd.DataFrame([{
                "#": i + 1, "서열": c.sequence, "형식": c.kind, "SEQ ID NO": c.seq_id_no, "길이": len(c.sequence),
                "경고": len(c.warnings), "문맥": c.context,
            } for i, c in enumerate(candidates)]),
            hide_index=True, width="stretch",
        )
        labels = {i: f"#{i + 1} · {c.sequence[:24]}{'…' if len(c.sequence) > 24 else ''}"
                     + (f" · SEQ ID NO: {c.seq_id_no}" if c.seq_id_no else "") for i, c in enumerate(candidates)}
        picked = st.multiselect("채택할 후보 — 여러 개를 고를 수 있습니다", list(labels), format_func=labels.get,
                                key=f"{key}_cand", placeholder="후보를 선택하세요")
        all_on = st.checkbox(f"후보 {len(candidates)}개 전부 선택 (한 번에 최대 {MAX_BULK}개)", key=f"{key}_all")
        chosen = list(range(len(candidates)))[:MAX_BULK] if all_on else picked
        if len(chosen) == 1:
            _adopt_form(candidates[chosen[0]], prov, suggestions, key=f"{key}_{chosen[0]}", name_hint=default_name_hint)
        elif len(chosen) > 1:
            _bulk_adopt_form([(candidates[i], prov, suggestions, default_name_hint) for i in chosen], key=f"{key}_bulk")
        else:
            st.caption("하나를 고르면 자세히 확인·수정하며 채택하고, 여러 개를 고르면 한꺼번에 채택합니다.")

    if suggestions:
        with st.expander(f"결합 대상 후보 제안 ({len(suggestions)}) — 소재 키워드 매칭"):
            for s in suggestions:
                st.markdown(f"- **{s.name}** ({s.category}, {s.mentions}회) — “{s.evidence}”")
    else:
        st.caption("텍스트에서 소재 키워드가 발견되지 않아 결합 대상은 '미보고'가 기본값입니다.")


MAX_BULK = 200
BT_AUTO = "각 논문에서 제안된 1순위 후보를 \u0027제안\u0027 상태로 기록"
BT_SAME = "선택한 서열 모두 같은 대상으로 직접 입력"
BT_NONE = "미보고로 두기"


def _bulk_adopt_form(items: list[tuple[Candidate, Provenance, list, str]], *, key: str) -> None:
    """여러 서열 후보를 한꺼번에 채택한다. 서열마다 출처는 각자의 논문, 발굴 방법은 공통, 결합 대상은 방식 선택."""
    n = len(items)
    st.markdown(f"##### 선택한 서열 {n}개를 한꺼번에 채택")
    st.dataframe(pd.DataFrame([{
        "서열": c.sequence, "길이": len(c.sequence), "출처": f"{p.identifier}" + (f" · SEQ ID NO: {c.seq_id_no}" if c.seq_id_no else ""),
        "제안된 결합 대상": ", ".join(t.name for t in tg[:1]) or "미보고",
        "경고": "⚠" if c.warnings else "",
    } for c, p, tg, _ in items]), hide_index=True, width="stretch", height=min(300, 36 * n + 40))
    skipped = [c.sequence for c, *_ in items if not c.valid]
    if skipped:
        st.warning(f"표준 20 아미노산이 아닌 문자(X 등)가 있는 후보 {len(skipped)}개는 건너뜁니다. "
                   "하나씩 선택해 서열을 고친 뒤 채택하세요.")
    st.caption("서열마다 출처는 그 서열이 나온 논문으로 기록되고, 비표준 요소 메모도 서열별로 보존됩니다.")
    with st.expander("발굴 방법 · 결합 대상 (선택한 서열 모두에 적용)", expanded=True):
        methods = st.multiselect("발굴(스크리닝) 방법 — 논문을 읽고 직접 선택", SCREENING_METHODS, key=f"{key}_methods",
                                 placeholder="선택하지 않으면 \u0027미보고\u0027로 기록됩니다")
        note = st.text_input("발굴 방법 메모", key=f"{key}_method_note")
        mode = st.radio("결합 대상", [BT_AUTO, BT_SAME, BT_NONE], key=f"{key}_bt_mode")
        same_name = same_cat = same_note = ""
        same_confirm = False
        if mode == BT_SAME:
            same_name = st.text_input("결합 대상 이름 (예: 폴리에틸렌 테레프탈레이트)", key=f"{key}_bt_name")
            same_cat = st.selectbox("분류", BINDING_CATEGORIES, key=f"{key}_bt_cat")
            same_note = st.text_input("형태·조건 메모", key=f"{key}_bt_note")
            same_confirm = st.checkbox("논문에서 직접 확인했으므로 \u0027확정\u0027으로 기록", key=f"{key}_bt_confirm")
        else:
            st.caption("\u0027제안\u0027은 소재 키워드 매칭 결과라 틀릴 수 있습니다. 채택 후 사이드바에서 서열마다 확정·수정하세요.")
    ready = n - len(skipped)
    if st.button(f"선택한 {ready}개 채택", type="primary", key=f"{key}_go", disabled=ready == 0):
        added = []
        for c, p, tg, _hint in items:
            try:
                seq = validate_sequence(c.sequence)
            except SequenceError:
                continue
            ident = p.identifier
            if c.seq_id_no:
                ident = f"{ident} · SEQ ID NO: {c.seq_id_no}".strip(" ·")
            name = f"{seq[:14]}" + (f" (SEQ ID NO {c.seq_id_no})" if c.seq_id_no else "") + (f" · {p.identifier}" if p.identifier else "")
            rec = SeqRecord(name=unique_name(name), sequence=seq,
                            provenance=Provenance(p.source_type, ident, p.url, title=p.title),
                            nonstandard_notes=c.nonstandard_notes)
            if mode == BT_AUTO and tg:
                target = BindingTarget(tg[0].name, tg[0].category, "", STATUS_SUGGESTED)
            elif mode == BT_SAME and same_name.strip():
                target = BindingTarget(same_name.strip(), same_cat, same_note, STATUS_CONFIRMED if same_confirm else STATUS_SUGGESTED)
            else:
                target = BindingTarget()
            apply_context(rec, methods, note, target)
            add_record(rec, make_current=False)
            added.append(rec)
        if added:
            set_current(added[0].name)
        flash(f"{len(added)}개 서열을 채택했습니다. 사이드바 \u0027내 서열\u0027에서 확인하세요." + (f" (건너뜀 {len(skipped)}개)" if skipped else ""))
        st.rerun()


def _adopt_form(cand: Candidate, prov: Provenance, suggestions, *, key: str, name_hint: str) -> None:
    sig = hashlib.md5(cand.sequence.encode()).hexdigest()[:6]
    k = f"{key}_{sig}"
    for w in cand.warnings:
        st.warning(w)
    seq_text = st.text_area("서열 (필요하면 수정)", value=cand.sequence, key=f"{k}_seq", height=80)
    seq_ok = None
    try:
        seq_ok = validate_sequence(seq_text)
    except SequenceError as exc:
        st.error(str(exc))
    notes = st.text_area("비표준 요소 메모 (원문 표기 보존)", value=cand.nonstandard_notes, key=f"{k}_notes", height=70)
    default_name = f"{cand.sequence[:14]}" + (f" (SEQ ID NO {cand.seq_id_no})" if cand.seq_id_no else "") + (f" · {name_hint}" if name_hint else "")
    name = st.text_input("저장할 이름", value=default_name.strip(), key=f"{k}_name")
    with st.expander("발굴 방법 · 결합 대상 기록", expanded=bool(suggestions)):
        methods, note, target = context_editor(f"{k}_ctx", suggestions=suggestions)
    if st.button("이 후보를 채택", type="primary", key=f"{k}_adopt", disabled=seq_ok is None):
        ident = prov.identifier
        if cand.seq_id_no:
            ident = f"{ident} · SEQ ID NO: {cand.seq_id_no}".strip(" ·")
        rec = SeqRecord(
            name=unique_name(name), sequence=seq_ok,
            provenance=Provenance(prov.source_type, ident, prov.url, title=prov.title),
            nonstandard_notes=notes.strip(),
        )
        apply_context(rec, methods, note, target)
        add_record(rec)
        flash(f"'{rec.name}'을(를) 채택했습니다. 출처: {rec.provenance.short()}")
        st.rerun()


# ------------------------------------------------------------ 직접 입력 · 프리셋
def _manual_section() -> None:
    st.markdown("**예시 프리셋** — 서열은 UniProt에서 직접 확인한 값입니다.")
    c1, c2 = st.columns([4, 1])
    name = c1.selectbox("프리셋", preset_names(), key="preset_name", label_visibility="collapsed")
    p = get_preset(name)
    c1.caption(f"{p['description']} · UniProt {p['accession']} ({p['start']}-{p['end']}, {len(p['sequence'])} aa) — {p['note']}")
    if c2.button("불러오기", key="preset_load", width="stretch"):
        add_record(preset_record(name))
        flash(f"프리셋 '{name}'을(를) 추가했습니다.")
        st.rerun()

    st.divider()
    st.markdown("**직접 입력 / FASTA 업로드**")
    uploaded = st.file_uploader("FASTA 파일 (선택)", type=["fasta", "fa", "faa", "txt"], key="manual_fasta")
    default_name, default_seq = "", ""
    if uploaded is not None:
        entries = parse_fasta(uploaded.getvalue().decode("utf-8", errors="replace"))
        if entries:
            idx = st.selectbox("파일 안의 서열", range(len(entries)), format_func=lambda i: f"{entries[i][0]} ({len(entries[i][1])} aa)",
                               key="manual_fasta_idx") if len(entries) > 1 else 0
            default_name, default_seq = entries[idx]
        else:
            st.warning("FASTA에서 서열을 읽지 못했습니다.")
    sig = hashlib.md5((default_name + default_seq).encode()).hexdigest()[:6]
    seq_text = st.text_area("서열 (1문자 표기, 표준 20 아미노산)", value=default_seq, key=f"manual_seq_{sig}", height=90)
    name_in = st.text_input("이름", value=default_name, key=f"manual_name_{sig}")
    seq_ok = None
    if seq_text.strip():
        try:
            seq_ok = validate_sequence(seq_text)
            st.caption(f"{len(seq_ok)} aa — 검증 통과")
        except SequenceError as exc:
            st.error(str(exc))
    c1, c2 = st.columns(2)
    src = c1.selectbox("출처 유형", [SOURCE_MANUAL, SOURCE_PAPER, SOURCE_PATENT, "기타"], key=f"manual_src_{sig}")
    ident = c2.text_input("식별자 (PMID, DOI, 특허번호 …)", key=f"manual_ident_{sig}")
    notes = st.text_input("비표준 요소 메모 (선택)", key=f"manual_notes_{sig}")
    with st.expander("발굴 방법 · 결합 대상 기록 (선택)"):
        methods, note, target = context_editor(f"manual_ctx_{sig}")
    if st.button("서열 추가", type="primary", disabled=seq_ok is None or not name_in.strip(), key=f"manual_add_{sig}"):
        rec = SeqRecord(
            name=unique_name(name_in), sequence=seq_ok,
            provenance=Provenance(src, ident.strip()), nonstandard_notes=notes.strip(),
        )
        apply_context(rec, methods, note, target)
        add_record(rec)
        flash(f"'{rec.name}'을(를) 추가했습니다.")
        st.rerun()

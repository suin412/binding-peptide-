"""Streamlit 세션 상태: 서열 레코드 저장소와 현재 선택."""
from __future__ import annotations

import streamlit as st

from ..provenance import SeqRecord


def init_state() -> None:
    st.session_state.setdefault("records", {})       # dict[str, SeqRecord]
    st.session_state.setdefault("current", None)     # 현재 선택한 레코드 이름
    st.session_state.setdefault("structures", {})    # (source, sequence) -> dict
    st.session_state.setdefault("flash", [])         # 다음 실행에서 한 번 보여줄 메시지
    # 사이드바 라디오는 위젯 상태로 선택을 기억하므로, 코드에서 현재 서열을 바꾸면 다음 실행의 맨 앞(위젯 생성 전)에서 반영한다
    pending = st.session_state.pop("_pending_current", None)
    if pending in st.session_state["records"]:
        st.session_state["current"] = pending
        st.session_state["sidebar_pick"] = pending


def records() -> dict[str, SeqRecord]:
    return st.session_state["records"]


def current_record() -> SeqRecord | None:
    name = st.session_state.get("current")
    return records().get(name) if name else None


RUN_LOG_LIMIT = 100


def log_run(kind: str, state: dict) -> None:
    """추천 실행 기록을 쌓는다 (kind: 'imp' | 'rec'). 엑셀 내보내기가 지금까지의 모든 실행을 한 파일에 담는 데 쓴다."""
    from datetime import datetime

    log = st.session_state.setdefault("run_log", [])
    log.append({"kind": kind, "time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"), **state})
    del log[:-RUN_LOG_LIMIT]


def unique_name(base: str) -> str:
    base = (base or "서열").strip() or "서열"
    if base not in records():
        return base
    i = 2
    while f"{base} ({i})" in records():
        i += 1
    return f"{base} ({i})"


def add_record(rec: SeqRecord, *, make_current: bool = True) -> SeqRecord:
    """이름 충돌 시 번호를 붙여 저장한다. 변이체의 부모 이름은 호출 전에 확정돼 있어야 한다."""
    rec.name = unique_name(rec.name)
    records()[rec.name] = rec
    if make_current:
        _make_current(rec.name)
    return rec


def _make_current(name: str) -> None:
    st.session_state["current"] = name
    st.session_state["_pending_current"] = name


def set_current(name: str) -> None:
    if name in records():
        _make_current(name)


def delete_record(name: str) -> None:
    records().pop(name, None)
    if st.session_state.get("current") == name:
        names = list(records())
        st.session_state["current"] = names[-1] if names else None


def flash(message: str, kind: str = "success") -> None:
    st.session_state["flash"].append((kind, message))


def show_flash() -> None:
    for kind, message in st.session_state.get("flash", []):
        getattr(st, kind, st.info)(message)
    st.session_state["flash"] = []

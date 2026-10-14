"""전체 결과를 엑셀(.xlsx) 한 파일로 내보낸다. 외부 라이브러리 없이 zip + XML로 직접 만든다(휴대용 배포에서도 그대로 동작)."""
from __future__ import annotations

import io
import math
import numbers
import re
import zipfile
from datetime import datetime
from xml.sax.saxutils import escape

from .affinity import rank_correlation
from .analysis.binding_score import CLASSES, class_for_material, score_sequence
from .analysis.properties import compute_properties
from .provenance import NOT_REPORTED, SeqRecord

Sheet = tuple[list[str], list[list]]      # (머리글, 행들)
_ILLEGAL_XML = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f]")
_BAD_SHEET_CHARS = re.compile(r"[\[\]:*?/\\]")


def _col(n: int) -> str:
    s = ""
    while n:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def _cell(ref: str, value, style: int = 0) -> str:
    if value is None or (isinstance(value, numbers.Real) and not isinstance(value, bool) and not math.isfinite(float(value))):
        return f'<c r="{ref}" s="{style}"/>' if style else ""
    if isinstance(value, bool):
        value = "예" if value else "아니오"
    if isinstance(value, numbers.Integral):
        return f'<c r="{ref}" s="{style}"><v>{int(value)}</v></c>'
    if isinstance(value, numbers.Real):                        # numpy 실수도 순수 숫자로 기록한다
        return f'<c r="{ref}" s="{style}"><v>{float(value)!r}</v></c>'
    text = escape(_ILLEGAL_XML.sub("", str(value)))
    return f'<c r="{ref}" s="{style}" t="inlineStr"><is><t xml:space="preserve">{text}</t></is></c>'


def _sheet_xml(headers: list[str], rows: list[list]) -> str:
    widths = [min(60, max(8, len(h) * 2 + 2, *(min(60, len(str(r[i])) + 2) for r in rows if i < len(r) and r[i] is not None))) for i, h in enumerate(headers)]
    cols = "".join(f'<col min="{i}" max="{i}" width="{w}" customWidth="1"/>' for i, w in enumerate(widths, 1))
    out = [f'<row r="1">{"".join(_cell(_col(i) + "1", h, 1) for i, h in enumerate(headers, 1))}</row>']
    for n, row in enumerate(rows, 2):
        out.append(f'<row r="{n}">{"".join(_cell(_col(i) + str(n), v) for i, v in enumerate(row, 1))}</row>')
    return ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>'
            f'<cols>{cols}</cols><sheetData>{"".join(out)}</sheetData></worksheet>')


def xlsx_bytes(sheets: dict[str, Sheet]) -> bytes:
    """{시트 이름: (머리글, 행들)} → .xlsx 파일 내용. 시트 이름은 엑셀 규칙(31자, 특수문자 제외)에 맞게 다듬는다."""
    names: list[str] = []
    for raw in sheets:
        base = (_BAD_SHEET_CHARS.sub("-", raw).strip() or "시트")[:31]
        name, k = base, 2
        while name.lower() in (n.lower() for n in names):
            name = f"{base[:28]}_{k}"
            k += 1
        names.append(name)
    n = len(names)
    sheet_overrides = "".join(f'<Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' for i in range(1, n + 1))
    content_types = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                     '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
                     '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
                     '<Default Extension="xml" ContentType="application/xml"/>'
                     '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
                     '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
                     f'{sheet_overrides}</Types>')
    root_rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                 '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                 '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>')
    workbook = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'
                + "".join(f'<sheet name="{escape(nm, {chr(34): "&quot;"})}" sheetId="{i}" r:id="rId{i}"/>' for i, nm in enumerate(names, 1))
                + '</sheets></workbook>')
    wb_rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
               '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
               + "".join(f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>' for i in range(1, n + 1))
               + f'<Relationship Id="rId{n + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>')
    styles = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
              '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
              '<fonts count="2"><font><sz val="11"/><name val="Malgun Gothic"/></font><font><b/><sz val="11"/><name val="Malgun Gothic"/></font></fonts>'
              '<fills count="3"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill>'
              '<fill><patternFill patternType="solid"><fgColor rgb="FFDCE6F4"/></patternFill></fill></fills>'
              '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
              '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
              '<cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
              '<xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1"/></cellXfs></styleSheet>')
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("_rels/.rels", root_rels)
        z.writestr("xl/workbook.xml", workbook)
        z.writestr("xl/_rels/workbook.xml.rels", wb_rels)
        z.writestr("xl/styles.xml", styles)
        for i, (headers, rows) in enumerate(sheets.values(), 1):
            z.writestr(f"xl/worksheets/sheet{i}.xml", _sheet_xml(headers, rows))
    return buf.getvalue()


# ---------------------------------------------------------------- 시트 구성
def records_sheet(recs: dict[str, SeqRecord]) -> Sheet:
    headers = ["이름", "서열", "길이(aa)", "구분", "부모", "변이", "출처 유형", "출처 식별자", "출처 URL", "제목", "확인일",
               "발굴 방법", "결합 대상", "결합 대상 분류", "결합 대상 상태", "조건·메모", "보고된 친화도",
               "분자량(Da)", "pI", "순전하(pH 7.4)", "GRAVY", "불안정성 지수", "방향족성", "비표준 요소 메모", "설계 참고(어디서 왔나)"]
    rows = []
    for r in recs.values():
        p = compute_properties(r.sequence)
        ref = r.reference_info or {}
        origin = "; ".join(f"{k}={v}" for k, v in ref.items() if k in ("improved_from", "recommended_for", "recommend_score"))
        rows.append([
            r.name, r.sequence, p["length"], "설계된 서열(비실험)" if r.is_designed else "원본", r.parent or "", r.mutation,
            r.provenance.source_type, r.provenance.identifier, r.provenance.url, r.provenance.title, r.provenance.retrieved,
            r.methods_label(), r.binding_target.label(), r.binding_target.category, r.binding_target.status if r.binding_target.reported else NOT_REPORTED,
            r.binding_target.note, r.affinity_label(),
            p["mw"], p["pI"], p["net_charge_pH7.4"], p["gravy"], p["instability_index"], p["aromaticity"], r.nonstandard_notes, origin,
        ])
    return headers, rows


def affinity_sheet(recs: dict[str, SeqRecord]) -> Sheet:
    headers = ["서열 이름", "서열", "지표", "값", "단위", "구분", "측정 방법", "결합 대상", "측정 조건", "출처", "상태", "ΔG (kcal/mol, 298 K)", "메모"]
    rows = []
    for r in recs.values():
        for a in r.affinities:
            dg = a.delta_g
            rows.append([r.name, r.sequence, a.metric, a.value, a.unit, "계산값(측정 아님)" if a.computed else "보고된 측정값", a.method,
                         a.target, a.conditions, a.source, a.status, None if dg is None else round(dg, 2), a.note])
    return headers, rows


def scores_sheet(recs: dict[str, SeqRecord]) -> Sheet:
    """내 서열 전체의 표면별 결합 성향(가중치 척도 0–100). 서열의 결합 대상에서 추정한 표면은 백분위도 함께 적는다."""
    keys = list(CLASSES)
    headers = (["이름", "서열", "추정한 표면(결합 대상 기준)", "그 표면 백분위", "그 표면 가중치 척도"]
               + [f"{CLASSES[k].label.split(' (')[0]} 가중치 척도" for k in keys])
    rows = []
    for r in recs.values():
        guess = class_for_material(r.binding_target.name) if r.binding_target.reported else None
        all_scores = {k: score_sequence(r.sequence, k) for k in keys}
        g = all_scores.get(guess)
        rows.append([r.name, r.sequence, CLASSES[guess].label.split(" (")[0] if guess else NOT_REPORTED,
                     round(g.percentile, 1) if g else None, round(g.scaled, 1) if g else None] + [round(all_scores[k].scaled, 1) for k in keys])
    return headers, rows


def compare_sheet(recs: dict[str, SeqRecord]) -> tuple[Sheet, float | None]:
    """예측 점수 vs 보고된 친화도(ΔG) 비교 (결합 친화도 탭 ③과 같은 계산). 표면을 모르면 폴리머로 가정한다."""
    headers = ["이름", "서열", "표면 분류", "결합 성향(백분위)", "결합 성향(가중치 척도)", "보고된 값", "ΔG (kcal/mol)", "종류"]
    rows, xs, ys = [], [], []
    for name, rec in recs.items():
        usable = [a for a in rec.affinities if a.delta_g is not None]
        if not usable:
            continue
        aff = ([a for a in usable if not a.computed] or usable)[0]
        key = class_for_material(rec.binding_target.name) or class_for_material(aff.target) or "polymer"
        s = score_sequence(rec.sequence, key)
        rows.append([name, rec.sequence, CLASSES[key].label.split(" (")[0], round(s.percentile, 1), round(s.scaled, 1),
                     aff.label(), round(aff.delta_g, 2), "계산값" if aff.computed else "측정값"])
        xs.append(s.scaled)
        ys.append(-aff.delta_g)
    return (headers, rows), rank_correlation(xs, ys) if rows else None


def improve_sheet(runs: list[dict]) -> Sheet:
    """개선형 추천을 실행한 모든 기록을 한 시트에 쌓는다 (실행마다 '부모' 줄 + 개선형 후보 줄)."""
    headers = ["실행 번호", "실행 시각", "순위", "서열", "변경", "변이 수", "예상 결합 성향", "부모 대비", "표면 점수", "순전하", "GRAVY", "분자량(Da)",
               "표면 가중치 변화", "비고", "부모 이름", "부모 서열", "결합 표면"]
    rows = []
    for n, run in enumerate(runs, 1):
        res = run["result"]
        label = CLASSES[res.class_key].label
        when = run.get("time", "")
        rows.append([n, when, 0, res.parent, "(부모)", 0, res.parent_score, 0.0, res.parent_binding, None, None, None, "", "; ".join(res.parent_notes),
                     run.get("parent", ""), res.parent, label])
        for i, c in enumerate(res.candidates, 1):
            rows.append([n, when, i, c.sequence, " ".join(c.mutations), len(c.mutations), c.score, round(c.delta, 1), c.binding, c.charge, c.gravy, c.mw,
                         ", ".join(f"{k} {v:+.2f}" for k, v in c.weight_changes.items()), "; ".join(c.notes), run.get("parent", ""), res.parent, label])
    return headers, rows


def recommend_sheet(runs: list[dict]) -> Sheet:
    """알려진 전체 서열 기반 추천을 실행한 모든 기록을 한 시트에 쌓는다."""
    headers = ["실행 번호", "실행 시각", "결합 표면", "순위", "서열", "길이(aa)", "예상 결합 성향", "표면 점수", "무작위 대비 백분위", "알려진 서열과의 유사도",
               "가장 닮은 서열", "새로움(1=완전히 새로움)", "순전하", "GRAVY", "분자량(Da)", "비고", "탐색 시드"]
    rows = []
    for n, run in enumerate(runs, 1):
        res = run["result"]
        label = CLASSES[res.class_key].label
        for i, c in enumerate(res.candidates, 1):
            rows.append([n, run.get("time", ""), label, i, c.sequence, c.length, c.score, c.binding, c.percentile, c.similarity, c.nearest, c.novelty,
                         c.charge, c.gravy, c.mw, "; ".join(c.notes), res.seed_value])
    return headers, rows


def references_sheet(runs: list[dict]) -> Sheet:
    headers = ["실행 번호", "결합 표면", "참고한 알려진 서열", "서열", "출처", "구분", "이 점수에서의 기준선 점수"]
    rows = []
    for n, run in enumerate(runs, 1):
        res = run["result"]
        label = CLASSES[res.class_key].label
        rows += [[n, label, s.name, s.sequence, s.source, "문헌" if s.from_literature else "내 서열", res.baseline.get(s.name)] for s in res.seeds]
    return headers, rows


def literature_sheet(state: dict) -> Sheet:
    """검색 탭의 '논문·특허에서 서열 찾기' 결과: 서열 한 줄 + 대표 논문."""
    res = state["result"]
    headers = ["검색어", "서열", "길이(aa)", "나온 논문 수", "연도", "서열 문맥", "관련도", "찾은 곳", "대표 논문 제목", "대표 논문 연도", "논문 식별자", "논문 URL",
               "결합 대상 후보(제안, 미검증)", "경고"]
    rows = []
    for g in res.groups:
        b = g.best
        p = b.paper
        targets = "; ".join(t.name for t in b.targets)
        warns = "; ".join(dict.fromkeys(w for h in g.hits for w in h.candidate.warnings))
        years = g.years
        rows.append([state.get("query", ""), g.sequence, len(g.sequence), g.paper_count, f"{years[0]}-{years[-1]}" if years else "",
                     "불확실" if g.uncertain else "서열 문맥", round(g.score, 2), b.source, p.title, p.year, p.identifier, p.url, targets, warns])
    return headers, rows


def notes_sheet(recs: dict, imp_runs: list[dict], rec_runs: list[dict], lit_state) -> Sheet:
    lines = [
        ["파일 생성", datetime.now().strftime("%Y-%m-%d %H:%M")],
        ["프로그램", "PepLab (펩타이드·단백질 검색 → 서열 추천 → 분석 연습용 앱)"],
        ["중요", "추천 점수는 표면별 잔기 가중치로 만든 학습용 휴리스틱입니다. 결합 친화도(Kd) 예측이 아니며, 점수가 높다고 실제 결합이 좋다는 뜻이 아닙니다."],
        ["중요", "'설계된 서열'과 추천 결과는 실험으로 측정된 서열이 아닙니다. '미보고'는 값이 확인되지 않았다는 뜻입니다."],
        ["중요", "추천 시트에는 이번 접속 중 실행한 모든 추천이 '실행 번호' 순서로 한 시트에 쌓여 있습니다(앱을 닫으면 기록은 사라집니다)."],
        ["시트 '내 서열'", f"{len(recs)}개"],
        ["시트 '결합 친화도'", f"{sum(len(r.affinities) for r in recs.values())}건 (문헌 보고값과 계산값을 구분)"],
        ["시트 '개선형 추천'", f"{len(imp_runs)}회 실행 기록" if imp_runs else "없음 — 개선형 추천 탭에서 먼저 실행하세요"],
        ["시트 '전체 서열 기반 추천'", f"{len(rec_runs)}회 실행 기록" if rec_runs else "없음 — 전체 서열 기반 추천 탭에서 먼저 실행하세요"],
        ["시트 '논문 검색 결과'", f"서열 {len(lit_state['result'].groups)}개" if lit_state else "없음 — 검색 탭의 '논문·특허에서 찾기'를 먼저 실행하세요"],
    ]
    return ["항목", "내용"], lines


def build_export(recs: dict[str, SeqRecord], rec_state: dict | None = None, imp_state: dict | None = None, *,
                 runs: list[dict] | None = None, lit_state: dict | None = None) -> bytes:
    """지금까지의 모든 결과(내 서열·친화도·점수·모든 추천 실행 기록·논문 검색 결과)를 한 엑셀 파일로 묶는다.

    runs가 없으면 마지막 추천 결과(rec_state, imp_state)만 담는다.
    """
    if runs is None:
        runs = ([{"kind": "imp", **imp_state}] if imp_state else []) + ([{"kind": "rec", **rec_state}] if rec_state else [])
    imp_runs = [r for r in runs if r["kind"] == "imp"]
    rec_runs = [r for r in runs if r["kind"] == "rec"]
    compare, rho = compare_sheet(recs)
    notes = notes_sheet(recs, imp_runs, rec_runs, lit_state)
    if rho is not None:
        notes[1].append(["순위상관 ρ (점수 vs 결합 세기)", f"{rho:.2f} (n={len(compare[1])}; 표본이 적어 통계적 의미는 거의 없음)"])
    sheets: dict[str, Sheet] = {"안내": notes, "내 서열": records_sheet(recs), "결합 성향 점수": scores_sheet(recs),
                                "결합 친화도": affinity_sheet(recs), "예측 vs 보고값": compare}
    if imp_runs:
        sheets["개선형 추천"] = improve_sheet(imp_runs)
    if rec_runs:
        sheets["전체 서열 기반 추천"] = recommend_sheet(rec_runs)
        sheets["참고한 알려진 서열"] = references_sheet(rec_runs)
    if lit_state:
        sheets["논문 검색 결과"] = literature_sheet(lit_state)
    return xlsx_bytes(sheets)

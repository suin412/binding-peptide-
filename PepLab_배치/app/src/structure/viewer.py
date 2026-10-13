"""py3Dmol 3D 뷰어: pLDDT 색상, cartoon/stick, 변이 잔기 강조, 나란히/겹쳐 보기용 HTML 생성."""
from __future__ import annotations

import py3Dmol

from .pdbutil import PLDDT_BANDS, band, plddt_per_residue

HIGHLIGHT_COLOR = "#e91e63"


def _style_model(view, model_index: int, pdb_text: str, style: str, color_mode: str,
                 solid_color: str, highlight: list[int]) -> None:
    sel_model = {"model": model_index}
    shape = "cartoon" if style == "cartoon" else "stick"
    if color_mode == "plddt":
        groups: dict[str, list[int]] = {}
        for resi, value in plddt_per_residue(pdb_text):
            groups.setdefault(band(value)[1], []).append(resi)
        for color, resis in groups.items():
            view.setStyle({**sel_model, "resi": resis}, {shape: {"color": color}})
    else:
        view.setStyle(sel_model, {shape: {"color": solid_color}})
    if highlight:
        view.addStyle({**sel_model, "resi": highlight}, {"stick": {"color": HIGHLIGHT_COLOR, "radius": 0.25}})


def build_html(models: list[dict], *, style: str = "cartoon", width: int | str = "100%", height: int = 420) -> str:
    """models: [{'pdb': str, 'color_mode': 'plddt'|'solid', 'solid_color': str, 'highlight': [resi]}].

    width는 기본 '100%'(컬럼 폭에 맞춤)이며, iframe 안의 여백·스크롤바를 없애는 스타일을 함께 넣는다.
    """
    view = py3Dmol.view(width=width, height=height)
    for i, m in enumerate(models):
        view.addModel(m["pdb"], "pdb")
        _style_model(view, i, m["pdb"], style, m.get("color_mode", "plddt"),
                     m.get("solid_color", "#4c78a8"), m.get("highlight", []))
    view.zoomTo()
    return "<style>html,body{margin:0;padding:0;overflow:hidden}</style>" + view._make_html()


def legend_html() -> str:
    items = "".join(
        f'<span style="display:inline-block;margin-right:12px">'
        f'<span style="display:inline-block;width:12px;height:12px;background:{color};border-radius:2px;'
        f'vertical-align:middle;margin-right:4px"></span>{label} (≥{threshold})</span>'
        for threshold, label, color in PLDDT_BANDS
    )
    items += (f'<span style="display:inline-block"><span style="display:inline-block;width:12px;height:12px;'
              f'background:{HIGHLIGHT_COLOR};border-radius:2px;vertical-align:middle;margin-right:4px"></span>'
              f'변이 잔기</span>')
    return f'<div style="font-size:0.85rem">{items}</div>'

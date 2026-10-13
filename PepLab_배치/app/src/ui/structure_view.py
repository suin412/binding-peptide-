"""추천 탭이 공유하는 물성 비교 표. (구조 예측 화면은 없앴다 — src/structure/ 코드는 화면에서 쓰이지 않는다.)"""
from __future__ import annotations

import pandas as pd

from ..analysis.properties import compute_properties


def properties_table(items: list[tuple[str, str]], *, label_column: str = "서열") -> pd.DataFrame:
    rows = []
    for label, seq in items:
        p = compute_properties(seq)
        rows.append({label_column: label, "서열": seq, "분자량(Da)": p["mw"], "pI": p["pI"], "순전하(pH 7.4)": p["net_charge_pH7.4"],
                     "GRAVY": p["gravy"], "불안정성 지수": p["instability_index"], "방향족성": p["aromaticity"]})
    return pd.DataFrame(rows)

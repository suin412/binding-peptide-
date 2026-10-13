"""AlphaFold DB 조회 (UniProt accession 기준, 이미 예측된 구조)."""
from __future__ import annotations

from dataclasses import dataclass

from .http import ApiError, get_json, get_text

API = "https://alphafold.ebi.ac.uk/api/prediction"


@dataclass
class AlphaFoldModel:
    accession: str
    pdb_text: str
    sequence: str          # AlphaFold 모델이 가진 (UniProt 전체) 서열
    mean_plddt: float
    pdb_url: str


def fetch(accession: str) -> AlphaFoldModel:
    data = get_json(f"{API}/{accession.strip().upper()}", service="AlphaFold DB")
    if not data:
        raise ApiError("AlphaFold DB에 해당 항목의 예측 구조가 없습니다.")
    return parse_model(data, get_text)


def parse_model(data: list[dict], text_getter) -> AlphaFoldModel:
    item = data[0]
    pdb_url = item.get("pdbUrl")
    if not pdb_url:
        raise ApiError("AlphaFold DB 응답에 PDB 파일 주소가 없습니다.")
    pdb_text = text_getter(pdb_url, service="AlphaFold DB")
    return AlphaFoldModel(
        accession=item.get("uniprotAccession", ""),
        pdb_text=pdb_text,
        sequence=item.get("sequence") or item.get("uniprotSequence", ""),
        mean_plddt=float(item.get("globalMetricValue") or 0.0),
        pdb_url=pdb_url,
    )

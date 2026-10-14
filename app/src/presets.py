"""프리셋 로더 + 논문/특허 입력 연습용 가상 샘플 텍스트."""
from __future__ import annotations

from .presets_data import PRESETS
from .provenance import Provenance, SeqRecord, SOURCE_PRESET


# 프리셋 서열이 UniProt 항목 전체(성숙 구간이 아님)와 일치하는 accession
FULL_LENGTH_ACCESSIONS = {"P01504", "P42212"}


def preset_names() -> list[str]:
    return [p["name"] for p in PRESETS]


def get_preset(name: str) -> dict:
    for p in PRESETS:
        if p["name"] == name:
            return p
    raise KeyError(name)


def preset_record(name: str) -> SeqRecord:
    p = get_preset(name)
    prov = Provenance(
        source_type=SOURCE_PRESET,
        identifier=f"UniProt {p['accession']} ({p['start']}-{p['end']})",
        url=f"https://www.uniprot.org/uniprotkb/{p['accession']}",
        title=p["description"],
    )
    return SeqRecord(
        name=p["name"],
        sequence=p["sequence"],
        provenance=prov,
        nonstandard_notes=p["nonstandard_notes"],
        # 항목 전체가 프리셋과 같은 경우에만 AlphaFold DB를 쓸 수 있다 (서열 불일치는 구조 탭에서 한 번 더 검사)
        reference_info={"accession": p["accession"], "full_length": p["accession"] in FULL_LENGTH_ACCESSIONS},
    )


# 서열 추출·결합 대상 제안 기능을 연습하기 위한 **가상 텍스트**입니다. 실제 논문/특허가 아니며,
# 포함된 서열도 실존하는 펩타이드가 아닙니다 (기능 시연용).
SAMPLE_TEXTS = {
    "가상 논문 초록 (티타늄 결합 펩타이드)": (
        "In this practice abstract (fictional), we screened a 12-mer phage display library against "
        "titanium oxide surfaces. After three rounds of biopanning, the dominant clone displayed the "
        "peptide NQWDTSAHKVMR (named Ti-P1), which adsorbed onto titanium plates in PBS at pH 7.4. "
        "A control peptide, GSGSGSGSGSGS, showed no preference for polystyrene."
    ),
    "가상 특허 서열목록 (3문자 표기 + SEQ ID NO)": (
        "[Claim 1] An isolated peptide comprising the amino acid sequence of SEQ ID NO: 1. "
        "[Claim 2] The peptide of claim 1 that binds to gold nanoparticles.\n"
        "<210> 1\n<211> 10\n<212> PRT\n<213> Artificial Sequence\n<400> 1\n"
        "Ala Gly Pro Trp His Lys Arg Ser Thr Val 10\n"
        "<210> 2\n<211> 8\n<212> PRT\n<213> Artificial Sequence\n<400> 2\n"
        "Ac-Gly-Aib-Leu-Lys-Orn-Trp-Phe-Ser-NH2"
    ),
}

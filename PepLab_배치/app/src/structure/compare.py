"""구조 중첩 및 Cα RMSD (선택 기능). 원자 수가 같은 경우에만 계산한다."""
from __future__ import annotations

import io
from dataclasses import dataclass

import numpy as np
from Bio.PDB import PDBIO, PDBParser, Superimposer


@dataclass
class RmsdResult:
    rmsd: float | None
    reason: str = ""              # rmsd가 None일 때 이유
    aligned_pdb_b: str = ""       # b를 a에 겹친 PDB (겹쳐 보기용)


def _ca_atoms(structure):
    return [res["CA"] for res in structure.get_residues() if "CA" in res]


def compute_rmsd(pdb_a: str, pdb_b: str) -> RmsdResult:
    parser = PDBParser(QUIET=True)
    sa = parser.get_structure("a", io.StringIO(pdb_a))
    sb = parser.get_structure("b", io.StringIO(pdb_b))
    ca_a, ca_b = _ca_atoms(sa), _ca_atoms(sb)
    if len(ca_a) != len(ca_b):
        return RmsdResult(None, f"Cα 개수가 달라(원본 {len(ca_a)}, 비교 {len(ca_b)}) RMSD를 계산할 수 없습니다. "
                                "길이가 같은 변이체(점돌연변이 등)에서만 사용할 수 있습니다.")
    if not ca_a:
        return RmsdResult(None, "Cα 원자를 찾지 못했습니다.")
    sup = Superimposer()
    sup.set_atoms(ca_a, ca_b)
    sup.apply(list(sb.get_atoms()))
    buf = io.StringIO()
    io_ = PDBIO()
    io_.set_structure(sb)
    io_.save(buf)
    return RmsdResult(float(sup.rms), "", buf.getvalue())


def same_length(pdb_a: str, pdb_b: str) -> bool:
    parser = PDBParser(QUIET=True)
    return len(_ca_atoms(parser.get_structure("a", io.StringIO(pdb_a)))) == \
        len(_ca_atoms(parser.get_structure("b", io.StringIO(pdb_b))))

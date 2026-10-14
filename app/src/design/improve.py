"""채택한 서열 하나를 골라 그 서열의 **개선형**(몇 개 잔기만 바꾼 변이체)을 자동으로 추천한다.

부모 서열에서 최대 k개 위치를 치환한 모든 변이체 중 점수가 부모보다 오르는 것을 빔 탐색으로 찾는다.
점수는 추천(`recommend.py`)과 같은 학습용 점수(표면 가중치 + 기존 결합 서열과의 유사도 − 제약 벌점)이며,
바꾼 잔기 수만큼 작은 비용을 물려 "조금만 바꿔서 많이 오르는" 변이를 우선한다.
**Kd 개선을 예측하지 않는다.** 점수가 오르는 방향의 변이를 제안할 뿐, 실제 결합이 좋아지는지는 실험으로 확인해야 한다.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from ..analysis.binding_score import CLASSES, score_sequence
from ..analysis.properties import compute_properties
from ..sequence import STANDARD_AA, validate_sequence
from .recommend import LITERATURE_SEEDS, SeedSequence, _Scorer, penalties

MUTATION_COST = 2.0          # 바꾼 잔기 하나당 점수에서 빼는 값 (적게 바꾸는 쪽 선호)
MIN_HAMMING_BETWEEN_PICKS = 2      # 변이 2개 이하일 때 / 3개 이상이면 3 (서로 거의 같은 후보가 줄줄이 나오지 않게)


@dataclass
class Improvement:
    sequence: str
    mutations: list[str]          # 예: ["T9W", "P4H"] (1-based, 위치순)
    score: float                  # 최종 적합도 = 표면·유사도 점수 − 제약 벌점 − 변이 비용 (벌점이 크면 음수가 될 수 있다)
    binding: float                # 표면 가중치 척도 점수
    delta: float                  # 부모 대비 점수 변화
    similarity: float | None
    charge: float
    gravy: float
    mw: float
    contributions: list[float]
    weight_changes: dict[str, float]   # 변이별 표면 가중치 변화 (+이면 그 표면에 유리한 쪽으로 바뀜)
    notes: list[str] = field(default_factory=list)


@dataclass
class ImproveResult:
    class_key: str
    parent: str
    parent_score: float
    parent_binding: float
    parent_notes: list[str]
    candidates: list[Improvement]
    seeds: list[SeedSequence]
    evaluations: int
    searched: int                 # 평가한 서로 다른 변이체 수
    fixed_positions: list[int]

    @property
    def already_good(self) -> bool:
        """점수가 오르는 변이를 하나도 못 찾았다 = 이 서열은 이 점수 기준으로 이미 (국소) 최적."""
        return not self.candidates


def parse_positions(text: str, length: int) -> list[int]:
    """'1-3, 7' → [1, 2, 3, 7] (1-based). 범위를 벗어나거나 형식이 틀리면 ValueError."""
    out: set[int] = set()
    for part in (text or "").replace(" ", "").split(","):
        if not part:
            continue
        try:
            if "-" in part:
                lo, hi = (int(x) for x in part.split("-", 1))
                if lo > hi:
                    raise ValueError
                out.update(range(lo, hi + 1))
            else:
                out.add(int(part))
        except ValueError:
            raise ValueError(f"'{part}'을(를) 해석할 수 없습니다. 예: 1-3, 7") from None
    bad = [p for p in out if p < 1 or p > length]
    if bad:
        raise ValueError(f"위치 {sorted(bad)}이(가) 서열 범위(1-{length})를 벗어났습니다.")
    return sorted(out)


def _hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def improve(
    parent: str, class_key: str, *, n: int = 8, max_mutations: int = 3, fixed_positions: list[int] | None = None,
    seeds: list[SeedSequence] | None = None, use_literature_seeds: bool = True, similarity_weight: float = 0.2,
    kmer_weight: float = 0.4, mutation_cost: float = MUTATION_COST, beam: int = 100,
) -> ImproveResult:
    """부모 서열의 개선형 n개를 추천한다 (결정적: 같은 입력이면 같은 결과)."""
    if class_key not in CLASSES:
        raise KeyError(class_key)
    parent = validate_sequence(parent, min_len=5)
    if len(parent) > 40:
        raise ValueError("개선형 추천은 40 aa 이하 서열만 지원합니다.")
    if not 1 <= max_mutations <= 6:
        raise ValueError("최대 변이 수는 1~6이어야 합니다.")
    fixed = set(fixed_positions or [])
    free = [i for i in range(len(parent)) if (i + 1) not in fixed]
    if not free:
        raise ValueError("모든 위치가 고정돼 바꿀 수 있는 위치가 없습니다.")

    all_seeds = list(LITERATURE_SEEDS.get(class_key, [])) if use_literature_seeds else []
    all_seeds += [s for s in (seeds or []) if s.sequence not in {x.sequence for x in all_seeds}]
    scorer = _Scorer(class_key, all_seeds, w_similarity=similarity_weight, kmer_weight=kmer_weight, penalize_copies=False)
    weights = CLASSES[class_key].weights
    parent_fit = scorer.fitness(parent)

    # 빔 탐색: 한 단계마다 아직 안 바꾼 자유 위치를 19가지 잔기로 치환해 보고 상위 beam개만 다음 단계로
    found: dict[str, tuple[float, tuple[tuple[int, str], ...]]] = {}
    frontier: list[tuple[str, tuple[tuple[int, str], ...]]] = [(parent, ())]
    for _ in range(max_mutations):
        scored: dict[str, tuple[float, tuple[tuple[int, str], ...]]] = {}
        for seq, muts in frontier:
            used = {p for p, _ in muts}
            for pos in free:
                if pos in used:
                    continue
                for aa in sorted(STANDARD_AA):
                    if aa == parent[pos]:
                        continue
                    new = seq[:pos] + aa + seq[pos + 1:]
                    if new in found or new in scored:
                        continue
                    new_muts = muts + ((pos, aa),)
                    scored[new] = (scorer.fitness(new) - mutation_cost * len(new_muts), new_muts)
        if not scored:
            break
        ranked = sorted(scored.items(), key=lambda kv: (-kv[1][0], kv[0]))
        found.update(dict(ranked))
        frontier = [(seq, data[1]) for seq, data in ranked[:beam]]

    # 부모보다 점수가 오른 변이체만, 서로 최소 2개 위치가 다른 것부터 n개
    chosen: list[str] = []
    min_hamming = MIN_HAMMING_BETWEEN_PICKS if max_mutations < 3 else MIN_HAMMING_BETWEEN_PICKS + 1
    for seq, (fit, muts) in sorted(found.items(), key=lambda kv: (-kv[1][0], kv[0])):
        if fit <= parent_fit:
            break
        if all(_hamming(seq, c) >= min_hamming for c in chosen):
            chosen.append(seq)
        if len(chosen) >= n:
            break

    out: list[Improvement] = []
    for seq in chosen:
        fit, muts = found[seq]
        sc = score_sequence(seq, class_key)
        props = compute_properties(seq)
        sim, _ = scorer.similarity(seq) if all_seeds else (None, "")
        _, notes = penalties(seq, scorer.constraints)
        ordered = sorted(muts)
        labels = [f"{parent[p]}{p + 1}{aa}" for p, aa in ordered]
        out.append(Improvement(
            sequence=seq, mutations=labels, score=round(fit, 1), binding=round(sc.scaled, 1),
            delta=round(fit - parent_fit, 1), similarity=None if sim is None else round(sim, 3), charge=props["net_charge_pH7.4"],
            gravy=props["gravy"], mw=props["mw"], contributions=sc.contributions,
            weight_changes={lab: round(weights.get(aa, 0.0) - weights.get(parent[p], 0.0), 2) for lab, (p, aa) in zip(labels, ordered)},
            notes=notes,
        ))
    _, parent_notes = penalties(parent, scorer.constraints)
    return ImproveResult(
        class_key=class_key, parent=parent, parent_score=round(parent_fit, 1),
        parent_binding=round(score_sequence(parent, class_key).scaled, 1), parent_notes=parent_notes, candidates=out,
        seeds=all_seeds, evaluations=scorer.evaluations, searched=len(found), fixed_positions=sorted(fixed),
    )

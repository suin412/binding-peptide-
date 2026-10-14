"""대상(표면)을 정하면 결합 성향이 높을 것 같은 서열을 자동으로 추천한다.

사용자가 변이·무작위 서열을 직접 만들 필요 없이, 내부에서 진화적 탐색(유전 알고리즘)으로 후보를 만들고 점수로 골라 준다.

점수는 `analysis.binding_score`의 **학습용 경험적 결합 성향 점수**다 (Kd 예측이 아님).
점수만 최대화하면 W만 12개 같은 비현실적 서열이 나오므로, 아래 제약을 벌점으로 걸어 현실적인 후보만 남긴다.
 - 한 잔기가 25%를 넘지 않게, 같은 잔기 4연속 금지, Cys 개수 제한
 - 용해도를 위해 전하 잔기(D·E·K·R)를 최소 2개, 순전하·GRAVY 범위 제한(표면 종류별로 다름)
같은 표면에 결합한다고 보고된 기존 서열(문헌 예시 + 내가 가진 서열)이 있으면 그 서열들의 잔기 구성과 닮을수록 가산하고,
기존 서열과 거의 같은 복사본은 제외한다(신규 후보만 추천).
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass, field

from Bio.SeqUtils.ProtParamData import kd

from ..analysis.binding_score import BACKGROUND, CLASSES, score_sequence
from ..analysis.properties import compute_properties
from ..sequence import STANDARD_AA

CHARGED = set("DEKR")
HYDROPHOBIC = set("WFYLIVM")
MIN_COMPOSITION_DISTANCE = 3      # 추천 후보끼리 잔기 구성이 최소 이만큼(잔기 수) 달라야 한다


@dataclass(frozen=True)
class Constraints:
    max_residue_frac: float = 0.25
    max_run: int = 3
    max_cys: int = 1
    min_charged: int = 2
    charge_range: tuple[float, float] = (-4.0, 6.0)
    gravy_range: tuple[float, float] = (-1.2, 0.9)


CLASS_CONSTRAINTS: dict[str, Constraints] = {
    "polymer": Constraints(),
    "titanium": Constraints(charge_range=(-2.0, 8.0)),
    "noble": Constraints(max_cys=3),
    "silica": Constraints(charge_range=(-1.0, 8.0)),
    "apatite": Constraints(charge_range=(-10.0, 1.0)),
    "carbon": Constraints(gravy_range=(-1.5, 1.2)),
}


@dataclass(frozen=True)
class SeedSequence:
    name: str
    sequence: str
    source: str            # 출처 설명 (PMID 등)
    from_literature: bool = True


# 문헌 본문에서 직접 확인한 결합 서열 (표면 분류별). 폴리머 4개는 계산으로 설계한 서열이다.
LITERATURE_SEEDS: dict[str, list[SeedSequence]] = {
    "titanium": [
        SeedSequence("TBP-1", "RKLPDAPGMHTW", "PMID 34948009"),
        SeedSequence("TiBP1", "RPRENRGRERGL", "PMID 34948009"),
        SeedSequence("TiBP2", "SRPNGYGGSESS", "PMID 34948009"),
    ],
    "polymer": [
        SeedSequence("PET 설계", "WEWWHVFHHRLR", "PMID 37735840 (계산 설계)"),
        SeedSequence("PE 설계", "HWMWAMKWHMRH", "PMID 37735840 (계산 설계)"),
        SeedSequence("PP 설계", "WWQRHMFHFRTW", "PMID 37735840 (계산 설계)"),
        SeedSequence("PS 설계", "WHWQRIIWQQMR", "PMID 37735840 (계산 설계)"),
    ],
}


@dataclass
class Recommendation:
    sequence: str
    score: float                 # 최종 예상 결합 성향 (0~100, 제약 위반 벌점 반영)
    binding: float               # 표면 가중치 척도 점수 (0~100)
    percentile: float            # 같은 길이 무작위 서열 대비 백분위
    similarity: float | None     # 기존 결합 서열과의 잔기 구성 유사도 (0~1), 기존 서열이 없으면 None
    nearest: str                 # 가장 닮은 기존 서열 이름
    novelty: float               # 1 - 기존 서열과의 최대 위치 동일성 (1이면 완전히 새로움)
    length: int
    charge: float
    gravy: float
    mw: float
    contributions: list[float]
    notes: list[str] = field(default_factory=list)


@dataclass
class RecommendResult:
    class_key: str
    candidates: list[Recommendation]
    seeds: list[SeedSequence]
    evaluations: int
    seed_value: int
    baseline: dict[str, float] = field(default_factory=dict)   # 기존(문헌·내 서열) 결합 서열 자신의 점수 - 점수가 얼마나 믿을 만한지 가늠하는 기준

    @property
    def weak_baseline(self) -> bool:
        """문헌에서 결합한다고 알려진 서열의 점수가 낮으면, 이 점수가 그 표면의 실제 결합 서열 특징을 잘 못 잡는다는 뜻이다."""
        lit = [self.baseline[s.name] for s in self.seeds if s.from_literature and s.name in self.baseline]
        return bool(lit) and sum(lit) / len(lit) < 45.0


# ---------------------------------------------------------------------------------------------
def quick_charge(seq: str) -> float:
    """탐색 중 빠르게 쓰는 순전하 근사 (pH 7.4: K·R +1, D·E -1, H +0.1)."""
    return seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E") + 0.1 * seq.count("H")


def quick_gravy(seq: str) -> float:
    return sum(kd[a] for a in seq) / len(seq)


def _longest_run(seq: str) -> int:
    best = run = 1
    for a, b in zip(seq, seq[1:]):
        run = run + 1 if a == b else 1
        best = max(best, run)
    return best


def penalties(seq: str, c: Constraints) -> tuple[float, list[str]]:
    """제약 위반 벌점과 사유 목록."""
    total, notes = 0.0, []
    n = len(seq)
    over = sum(max(0.0, seq.count(a) / n - c.max_residue_frac) for a in set(seq))
    if over > 0:
        total += 200 * over
        notes.append("특정 잔기가 25%를 넘음")
    run = _longest_run(seq)
    if run > c.max_run:
        total += 10 * (run - c.max_run)
        notes.append(f"같은 잔기 {run}연속")
    # 순서에 대한 약한 규칙: 소수성 잔기나 전하 잔기가 4개 이상 연달아 있으면 응집·용해도 문제 가능
    for label, group in (("소수성", HYDROPHOBIC), ("전하", CHARGED)):
        longest = cur = 0
        for a in seq:
            cur = cur + 1 if a in group else 0
            longest = max(longest, cur)
        if longest >= 4:
            total += 4 * (longest - 3)
            notes.append(f"{label} 잔기 {longest}연속")
    cys = seq.count("C")
    if cys > c.max_cys:
        total += 15 * (cys - c.max_cys)
        notes.append(f"Cys {cys}개")
    charged = sum(seq.count(a) for a in CHARGED)
    if charged < c.min_charged:
        total += 10 * (c.min_charged - charged)
        notes.append("전하 잔기가 적어 용해도가 낮을 수 있음")
    q = quick_charge(seq)
    lo, hi = c.charge_range
    if q < lo or q > hi:
        total += 5 * (lo - q if q < lo else q - hi)
        notes.append(f"순전하 {q:+.1f}가 권장 범위({lo:+g}~{hi:+g}) 밖")
    g = quick_gravy(seq)
    glo, ghi = c.gravy_range
    if g < glo or g > ghi:
        total += 20 * (glo - g if g < glo else g - ghi)
        notes.append(f"GRAVY {g:.2f}가 권장 범위({glo:g}~{ghi:g}) 밖")
    return total, notes


def _composition(seq: str) -> list[float]:
    return [seq.count(a) / len(seq) for a in sorted(STANDARD_AA)]


def _cosine(u: list[float], v: list[float]) -> float:
    dot = sum(a * b for a, b in zip(u, v))
    nu, nv = math.sqrt(sum(a * a for a in u)), math.sqrt(sum(b * b for b in v))
    return dot / (nu * nv) if nu and nv else 0.0


def composition_distance(a: str, b: str) -> int:
    """두 서열의 잔기 구성이 다른 정도 (한쪽을 다른 쪽 구성으로 바꾸는 데 필요한 최소 잔기 교체 수 근사)."""
    return sum(abs(a.count(x) - b.count(x)) for x in STANDARD_AA) // 2 + abs(len(a) - len(b)) // 2


def kmers(seq: str, k: int = 2) -> set[str]:
    """서열 안의 연속한 k개 잔기 조각(모티프) 집합."""
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def identity(a: str, b: str) -> float:
    """갭 없는 정렬에서 가장 높은 위치 동일성 (짧은 쪽 길이 기준, 0~1)."""
    short, long_ = (a, b) if len(a) <= len(b) else (b, a)
    best = 0
    for shift in range(len(long_) - len(short) + 1):
        best = max(best, sum(x == y for x, y in zip(short, long_[shift:])))
    return best / len(short)


class _Scorer:
    """같은 표면 분류에 대한 점수 계산을 반복 호출용으로 미리 준비해 둔다."""

    def __init__(self, class_key: str, seeds: list[SeedSequence], w_similarity: float = 0.25, *,
                 kmer_weight: float = 0.4, penalize_copies: bool = True):
        self.cls = CLASSES[class_key]
        self.constraints = CLASS_CONSTRAINTS.get(class_key, Constraints())
        w = self.cls.weights
        self.lo, self.hi = min(w.values()), max(w.values())
        self.seeds = seeds
        self.seed_comp = [(s, _composition(s.sequence)) for s in seeds]
        self.known_kmers: set[str] = set()
        for s in seeds:
            self.known_kmers |= kmers(s.sequence)
        self.kmer_weight = kmer_weight
        self.penalize_copies = penalize_copies
        self.w_similarity = w_similarity if seeds else 0.0
        self.evaluations = 0

    def scaled(self, seq: str) -> float:
        raw = sum(self.cls.weights.get(a, 0.0) for a in seq) / len(seq)
        return min(max(100.0 * (raw - self.lo) / (self.hi - self.lo), 0.0), 100.0)

    def similarity(self, seq: str) -> tuple[float, str]:
        """기존 결합 서열과의 유사도 (0~1) = 잔기 구성 코사인 유사도 + 2-mer(연속한 두 잔기) 모티프 포함 비율의 가중합.

        구성 유사도는 가장 닮은 기존 서열 하나와, 모티프 비율은 기존 서열 전체와 비교한다. 이름은 구성이 가장 닮은 기존 서열이다.
        """
        if not self.seed_comp:
            return 0.0, ""
        comp = _composition(seq)
        best_sim, best_name = -1.0, ""
        for s, sc in self.seed_comp:
            sim = _cosine(comp, sc)
            if sim > best_sim:
                best_sim, best_name = sim, s.name
        mine = kmers(seq)
        coverage = len(mine & self.known_kmers) / len(mine) if mine else 0.0
        return (1 - self.kmer_weight) * best_sim + self.kmer_weight * coverage, best_name

    def fitness(self, seq: str) -> float:
        self.evaluations += 1
        base = self.scaled(seq)
        if self.w_similarity:
            base = (1 - self.w_similarity) * base + self.w_similarity * 100.0 * self.similarity(seq)[0]
        pen, _ = penalties(seq, self.constraints)
        # 기존 서열과 거의 같은 복사본은 추천하지 않는다 (신규 후보만). 개선형 추천은 부모와 비슷한 것이 목적이라 끈다.
        if self.penalize_copies and self.seeds and max(identity(seq, s.sequence) for s in self.seeds) >= 0.8:
            pen += 50
        return base - pen


def _biased_choices(weights: dict[str, float], temperature: float = 2.5) -> tuple[list[str], list[float]]:
    aas = sorted(STANDARD_AA)
    probs = [BACKGROUND[a] * math.exp(temperature * weights.get(a, 0.0)) for a in aas]
    total = sum(probs)
    return aas, [p / total for p in probs]


def recommend(
    class_key: str, *, n: int = 10, min_len: int = 10, max_len: int = 14, seeds: list[SeedSequence] | None = None,
    use_literature_seeds: bool = True, seed_value: int = 0, population: int = 240, generations: int = 120,
    max_pair_identity: float = 0.65, similarity_weight: float = 0.4, kmer_weight: float = 0.4,
) -> RecommendResult:
    """표면 분류에 대해 결합 성향이 높을 것 같은 새 서열 n개를 추천한다 (재현 가능: seed_value가 같으면 결과가 같다)."""
    if class_key not in CLASSES:
        raise KeyError(class_key)
    if not 5 <= min_len <= max_len <= 40:
        raise ValueError("길이 범위는 5~40 사이이고 최소 ≤ 최대여야 합니다.")
    all_seeds = list(LITERATURE_SEEDS.get(class_key, [])) if use_literature_seeds else []
    all_seeds += [s for s in (seeds or []) if s.sequence not in {x.sequence for x in all_seeds}]
    rng = random.Random(seed_value)
    scorer = _Scorer(class_key, all_seeds, w_similarity=similarity_weight, kmer_weight=kmer_weight)
    aas, probs = _biased_choices(scorer.cls.weights)

    def biased_residue() -> str:
        return rng.choices(aas, probs)[0]

    def random_seq() -> str:
        return "".join(biased_residue() for _ in range(rng.randint(min_len, max_len)))

    def mutate_seq(seq: str, rate: float | None = None) -> str:
        chars = list(seq)
        rate = rate if rate is not None else 1.5 / len(chars)
        for i in range(len(chars)):
            if rng.random() < rate:
                chars[i] = biased_residue() if rng.random() < 0.7 else rng.choice(aas)
        if rng.random() < 0.1 and len(chars) < max_len:
            chars.insert(rng.randrange(len(chars) + 1), biased_residue())
        if rng.random() < 0.1 and len(chars) > min_len:
            del chars[rng.randrange(len(chars))]
        return "".join(chars)

    def crossover(a: str, b: str) -> str:
        cut = rng.randint(1, max(1, min(len(a), len(b)) - 1))
        child = a[:cut] + b[cut:]
        child = child[:max_len]
        return child if len(child) >= min_len else child + "".join(biased_residue() for _ in range(min_len - len(child)))

    # 초기 집단: 편향 무작위 + (기존 서열이 있으면) 그 서열들의 변이체
    pop = [random_seq() for _ in range(population)]
    for i, s in enumerate(all_seeds):
        base = s.sequence[:max_len]
        for j in range(min(15, population // max(len(all_seeds), 1) // 4)):
            pop[(i * 15 + j) % population] = mutate_seq(base, rate=0.25)
    if len(all_seeds) >= 2:                      # 알려진 서열들을 섞은(재조합) 후보도 출발점에 넣어 서로의 모티프를 이어 붙여 본다
        for _ in range(min(40, population // 4)):
            a, b = rng.sample(all_seeds, 2)
            pop[rng.randrange(population)] = crossover(a.sequence[:max_len], b.sequence[:max_len])
    seen: dict[str, float] = {}

    def evaluate(seq: str) -> float:
        if seq not in seen:
            seen[seq] = scorer.fitness(seq)
        return seen[seq]

    elite_n = max(4, population // 20)
    for _ in range(generations):
        scored = sorted(((evaluate(s), s) for s in pop), reverse=True)
        nxt = [s for _, s in scored[:elite_n]]

        def pick() -> str:
            return max(rng.sample(scored, 3))[1]

        while len(nxt) < population:
            child = crossover(pick(), pick()) if rng.random() < 0.6 else pick()
            nxt.append(mutate_seq(child))
        pop = nxt
    for s in pop:
        evaluate(s)

    # 점수가 잔기 구성에만 의존하므로 같은 구성을 순서만 바꾼 후보는 하나로 취급한다:
    # 서로 구성이 최소 MIN_COMPOSITION_DISTANCE 이상 다르고 위치 동일성도 낮은 n개를 점수 순으로 고른다.
    chosen: list[str] = []
    for seq, fit in sorted(seen.items(), key=lambda kv: (-kv[1], kv[0])):
        if fit < 0:
            break
        if all(composition_distance(seq, c) >= MIN_COMPOSITION_DISTANCE and identity(seq, c) < max_pair_identity for c in chosen):
            chosen.append(seq)
        if len(chosen) >= n:
            break

    out: list[Recommendation] = []
    for seq in chosen:
        sc = score_sequence(seq, class_key)
        props = compute_properties(seq)
        sim, nearest = scorer.similarity(seq) if all_seeds else (None, "")
        pen, notes = penalties(seq, scorer.constraints)
        novelty = 1.0 - max((identity(seq, s.sequence) for s in all_seeds), default=0.0)
        out.append(Recommendation(
            sequence=seq, score=round(max(seen[seq], 0.0), 1), binding=round(sc.scaled, 1), percentile=round(sc.percentile, 1),
            similarity=None if sim is None else round(sim, 3), nearest=nearest, novelty=round(novelty, 2), length=len(seq),
            charge=props["net_charge_pH7.4"], gravy=props["gravy"], mw=props["mw"], contributions=sc.contributions, notes=notes,
        ))
    baseline = {s.name: round(score_sequence(s.sequence, class_key).scaled, 1) for s in all_seeds}
    return RecommendResult(class_key, out, all_seeds, scorer.evaluations, seed_value, baseline)

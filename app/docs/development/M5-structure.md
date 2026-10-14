# M5 — 구조 확인 (AlphaFold DB, ESMFold, 3D 뷰어, 선택: RMSD)

## 목표
원본 단백질은 AlphaFold DB에서, 변이체와 DB에 없는 서열은 ESMFold로 구조를 얻어 3D로 나란히 보고, 신뢰도(pLDDT)를 확인한다. (PRD 시나리오 1·3·4)

## 범위
- 포함: AlphaFold DB 조회, ESMFold 예측, py3Dmol 뷰어, pLDDT 색상·경고, 나란히 보기, PDB 다운로드, 구조 기반 재순위화(선택, M4 연결), RMSD(선택)
- 제외: RCSB 실험 구조 로드(P2), 로컬 모델 실행, 삽입·삭제 변이체의 RMSD(원자 매핑 필요 → 향후 확장)

## 작업 목록
- [ ] 시작 전: AlphaFold DB(`alphafold.ebi.ac.uk/api/prediction/{uniprot}`)와 ESMFold(`POST api.esmatlas.com/foldSequence/v1/pdb/`) 실제 응답 확인, 샘플 PDB를 `tests/fixtures/`에 저장. 설치 확인: `py3Dmol`, `stmol` 렌더링
- [ ] `src/sources/alphafold.py`: UniProt AC → `pdbUrl` → PDB 텍스트, pLDDT(B-factor 열) 파싱. 캐시 적용
- [ ] `src/structure/esmfold.py`
  - 입력 검증: 표준 20 아미노산만, 길이 상한(약 400 aa) 초과 시 호출 차단 + 안내 메시지 (시나리오 4)
  - POST → PDB 문자열, 타임아웃·재시도, 비동기 진행 표시(`st.spinner`)
  - 서열 단위 캐시(해시 키)
  - 사용자 서열이 외부 API로 전송됨을 구조 탭에서 **한 번 고지**하고, 미공개 서열 입력 주의 문구 표시
- [ ] `src/structure/viewer.py`: py3Dmol/stmol 3D 뷰어 — pLDDT 색상(90/70/50 구간), cartoon/stick 전환, 변이 잔기 강조(변이 위치 목록에서 계산)
- [ ] 구조 탭 UI: 구조 소스 자동 선택(UniProt AC 있으면 AlphaFold DB, 없으면 ESMFold, 사용자가 전환 가능), 평균 pLDDT 표시, 두 구조 나란히 보기 (S-4)
- [ ] pLDDT 경고: 길이 30 aa 이하이거나 평균 pLDDT가 낮으면 "예측 신뢰도 낮음" 경고 + 학습용 설명 (S-5)
- [ ] 예측 PDB 다운로드 버튼 (S-6, P1)
- [ ] **(선택, 기본 꺼짐) RMSD** (S-7, P2)
  - 사용자가 켠 경우에만 계산, Biopython `Superimposer`, Cα 기준
  - 두 구조의 Cα 개수가 같을 때만 활성화. 다르면 비활성화 + "길이가 달라 계산할 수 없음" 안내
- [ ] **(선택) 구조 기반 재순위화** (A-7): M4 스크리닝 순위표에서 토글을 켜면 상위 N개(기본 3)를 ESMFold로 예측해 부모 대비 pLDDT 비교, 순위에 반영. 끈 상태에서도 기존 흐름이 동일하게 동작
- [ ] ESMFold 호출 실패/제한 시 앱 흐름이 끊기지 않도록 오류 메시지와 재시도 버튼

## 산출물
`src/sources/alphafold.py`, `src/structure/esmfold.py`, `viewer.py`, 구조 탭 UI, fixture(샘플 PDB), 테스트

## 완료 기준
- 원본(AlphaFold DB 또는 ESMFold)과 변이체(ESMFold)를 3D로 나란히 보고 pLDDT가 색으로 표시된다 (시나리오 3)
- 400 aa 초과 서열 입력 시 호출 없이 제한 안내가 표시된다 (시나리오 4)
- 짧은 펩타이드에서 낮은 신뢰도 경고가 표시된다
- RMSD를 켜면 길이가 같은 변이체에서만 계산되고, 길이가 다르면 비활성화 이유가 표시된다
- 구조 기반 재순위화를 켜고 상위 3개만 예측해 pLDDT 비교가 표시되며, 꺼도 정상 동작한다 (시나리오 3b)

## 테스트
- PDB 파싱: fixture에서 pLDDT 추출, 잔기 수 확인
- `esmfold.py`: 길이 상한·비표준 문자 차단, HTTP 오류·타임아웃 mock
- RMSD: 동일 구조 → 0, 원자 수 다른 구조 → 비활성(예외가 아닌 상태 반환)
- 뷰어 렌더링은 수동 점검 (브라우저에서 시각 확인)

## 리스크
- ESMFold 공개 API 장애·속도·길이 제한 → 캐시, 재시도, 안내 메시지, 필요하면 호출 지연 표시
- 짧은 펩타이드는 신뢰도가 낮은 것이 정상 → 경고를 학습 포인트로 설명
- Windows에서 3D 뷰어 컴포넌트 호환 문제 → 문제 시 `st.components.v1.html`로 py3Dmol HTML을 직접 임베드
- 사용자 서열의 외부 전송 → 고지 문구

# M2 — UniProt 검색 + 출처 기록

## 목표
UniProt에서 단백질/펩타이드를 검색해 서열과 메타데이터를 가져오고, 출처를 `SeqRecord`에 기록해 서열 헤더와 물성 표시까지 연결한다. (PRD 시나리오 1)

## 범위
- 포함: UniProt 검색·상세·참고문헌, 결과 표 UI, 출처 기록, 기본 물성 표시
- 제외: 논문·특허 검색(M3), 변이 설계(M4), 구조(M5)

## 작업 목록
- [ ] 시작 전: UniProt REST 응답을 실제로 한 번 호출해 JSON 필드 확인 (`rest.uniprot.org/uniprotkb/search`, `/uniprotkb/{acc}.json`, `.fasta`) 후 샘플을 `tests/fixtures/`에 저장
- [ ] `src/sources/uniprot.py`
  - `search(query, size)`: ID·키워드·유전자명 검색 → 이름, 생물종, 길이, AC 목록
  - `fetch(accession)`: 서열, 기능 요약, 관련 PDB/AlphaFold 링크, 참고문헌(PMID·제목)
  - 참고 정보로 `proteinExistence`, `Direct protein sequencing` 키워드 추출 (표시는 "참고 정보"로만, P2)
  - 타임아웃·재시도(2회)·`st.cache_data`(TTL 1일) 적용, 실패 시 사용자 친화적 메시지
- [ ] 검색 탭 UI: 검색창 → 결과 표 → 항목 선택 → 상세 패널(FASTA, 기능 요약, 링크, 참고문헌 PMID 링크) → "이 서열 사용" 버튼
- [ ] 선택 시 `SeqRecord` 생성: `provenance`에 유형=UniProt, 식별자=AC, URL, 확인 일자 기록
- [ ] **발굴 방법 선택 UI**: 목록 복수 선택 + 메모 칸, 기본 "미보고" (공통 컴포넌트로 만들어 M3에서도 재사용)
- [ ] `src/analysis/properties.py`: 길이, 분자량, pI, 순전하(pH 7.4), GRAVY — Biopython `ProteinAnalysis`. 순전하는 직접 계산(구현 시 방식을 코드 주석에 명시)
- [ ] 선택한 서열의 물성 표 표시
- [ ] 검색 결과 없음 안내: "논문에서 찾아보기" 버튼(M3 연결 자리, 이 마일스톤에서는 비활성 + 안내)

## 산출물
`src/sources/uniprot.py`, `src/analysis/properties.py`, 검색 탭 UI, 발굴 방법 선택 컴포넌트, `tests/fixtures/uniprot_*.json`, 관련 테스트

## 완료 기준
- "melittin" 검색 → 결과 표 → 선택 → 서열(26 aa 예상, 실제 값은 UniProt 기준) 확인 → 물성 표시
- 출처(UniProt AC)가 서열 헤더에 표시되고, 발굴 방법을 고르면 헤더에 반영된다
- API 실패 시 앱이 죽지 않고 안내 메시지가 표시된다
- 같은 검색을 반복하면 캐시로 즉시 응답한다

## 테스트
- `uniprot.py`: fixture 기반 파싱(필드 누락 시 안전하게 처리), 검색 결과 0건, HTTP 오류 mock
- `properties.py`: 알려진 서열의 분자량·pI를 기준값과 비교(허용 오차 명시), 순전하 부호 확인

## 리스크
- UniProt 응답 필드 변경 → fixture와 실제 응답을 구현 시작 시 대조
- 매우 긴 단백질의 물성 계산은 문제없으나 화면 표시가 길어짐 → 서열은 접기 처리

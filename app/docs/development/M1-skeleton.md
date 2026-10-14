# M1 — 골격 + 프리셋 + SeqRecord

## 목표
`streamlit run app.py`로 4개 탭이 뜨고, 프리셋 서열을 불러와 세션에 담을 수 있는 최소 뼈대를 만든다. 이후 모든 마일스톤이 공유하는 데이터 모델(`SeqRecord`)을 확정한다.

## 범위
- 포함: 가상환경, 의존성, 앱 진입점, 탭 4개(빈 화면 + 안내 문구), 프리셋 로드, 서열 직접 입력·FASTA 업로드와 검증, `SeqRecord`, 세션 상태 관리
- 제외: 외부 API 호출, 물성 계산, 변이, 구조 (M2 이후)

## 작업 목록
- [ ] 프로젝트 폴더 구조 생성 (`app.py`, `src/`, `tests/`, `data/cache/`) — [PLAN 프로젝트 구조](../PLAN_1test.md) 참고, 처음에는 모듈을 잘게 쪼개지 않고 필요할 때 분리
- [ ] `requirements.txt`: `streamlit`, `biopython`, `requests`, `pandas`, `matplotlib`, `py3Dmol`, `stmol`, `pytest`
- [ ] `.gitignore` (`data/cache/`, `.venv/`, `__pycache__/`) — git 저장소를 만들지는 사용자가 결정
- [ ] `src/provenance.py`: `SeqRecord` 데이터 클래스
  - 필드: `name`, `sequence`, `provenance`(유형, 식별자, URL, 확인 일자), `screening_method`(선택 목록 값 리스트 + 메모, 기본 "미보고"), `binding_target`(명칭, 분류, 메모, 상태: 제안/확정/미보고), `nonstandard_notes`, `parent`, `is_designed`
  - 변이체 생성 시 부모 참조를 설정하고 `is_designed=True`로 만든다
- [ ] 서열 검증 함수 `validate_sequence()`: 공백·줄바꿈 제거, 대문자화, 표준 20 아미노산 외 문자 거부 + 문제 위치 안내 (D-1)
- [ ] FASTA 파서 (Biopython `SeqIO`) + 업로드 UI
- [ ] `src/presets.py`: 프리셋 목록 (인슐린 B/A쇄, GLP-1(7-36), 멜리틴, 마가이닌, LL-37, 유비퀴틴, GFP). **서열은 구현 시 UniProt에서 직접 확인한 값만 사용**하고 출처(UniProt AC)를 함께 기록한다
- [ ] `app.py`: 탭 4개(검색, 설계, 분석·스크리닝, 구조), 사이드바에 "현재 서열 / 저장된 변이체 목록"
- [ ] 서열 헤더 컴포넌트: 이름, 길이, 출처 요약, 발굴 방법/결합 대상 요약 문장(없으면 "미보고"), 설계 서열 배지
- [ ] `st.session_state` 초기화 (`records: dict[str, SeqRecord]`, `current`)

## 산출물
`app.py`, `src/provenance.py`, `src/presets.py`, `requirements.txt`, `tests/test_validate.py`, `tests/test_seqrecord.py`

## 완료 기준
- `streamlit run app.py` 실행 시 오류 없이 4개 탭이 표시된다.
- 프리셋 선택 → 서열 헤더에 이름·길이·출처가 표시되고, 발굴 방법/결합 대상은 "미보고"로 표시된다.
- 잘못된 문자(`B`, `X`, 숫자 등)가 포함된 입력은 문제 위치와 함께 거부된다.
- 새로고침 전까지 저장한 서열이 사이드바에 유지된다.

## 테스트
- `validate_sequence`: 소문자/공백/줄바꿈 정리, 비표준 문자 거부, 빈 입력
- `SeqRecord`: 기본값 "미보고", 변이체의 `parent`·`is_designed`
- 프리셋: 모든 서열이 검증 함수를 통과

## 리스크
- 프리셋 서열을 기억에 의존해 쓰면 오타 위험 → UniProt에서 가져온 값만 사용
- Windows에서 `stmol`/`py3Dmol` 설치 문제 → M5 전에 설치만 미리 확인

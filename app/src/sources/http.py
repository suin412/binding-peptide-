"""공통 HTTP 헬퍼: 타임아웃·재시도·디스크 캐시·사용자 친화적 오류."""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import requests

CACHE_DIR = Path(__file__).resolve().parents[2] / "data" / "cache"
DEFAULT_TTL = 24 * 3600
USER_AGENT = "PepLab-practice-app/0.1 (educational)"


class ApiError(RuntimeError):
    """외부 API 호출 실패. 메시지는 사용자에게 그대로 보여줄 수 있다."""


def _key(method: str, url: str, params: dict | None, body: str | None) -> str:
    raw = json.dumps([method, url, sorted((params or {}).items()), body], ensure_ascii=False)
    return hashlib.sha1(raw.encode("utf-8")).hexdigest()


def _cache_path(key: str, suffix: str) -> Path:
    return CACHE_DIR / f"{key}.{suffix}"


def _read_cache(path: Path, ttl: int) -> str | None:
    try:
        if path.exists() and time.time() - path.stat().st_mtime < ttl:
            return path.read_text(encoding="utf-8")
    except OSError:
        pass
    return None


def _write_cache(path: Path, text: str) -> None:
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    except OSError:
        pass  # 캐시 실패는 치명적이지 않다


def _request(method: str, url: str, *, params=None, data=None, headers=None,
             timeout: int = 20, retries: int = 2, service: str = "외부 서비스") -> requests.Response:
    last: Exception | None = None
    for attempt in range(retries + 1):
        try:
            resp = requests.request(
                method, url, params=params, data=data,
                headers={"User-Agent": USER_AGENT, **(headers or {})}, timeout=timeout,
            )
            if resp.status_code == 404:
                raise ApiError(f"{service}에서 결과를 찾지 못했습니다 (404).")
            if resp.status_code >= 500 or resp.status_code == 429:
                raise requests.HTTPError(f"HTTP {resp.status_code}")
            resp.raise_for_status()
            return resp
        except ApiError:
            raise
        except requests.RequestException as exc:
            last = exc
            if attempt < retries:
                time.sleep(1.0 * (attempt + 1))
    raise ApiError(f"{service} 호출에 실패했습니다 (네트워크 또는 서비스 장애일 수 있습니다): {last}")


def get_json(url: str, params: dict | None = None, *, ttl: int = DEFAULT_TTL,
             service: str = "외부 서비스", timeout: int = 20) -> dict | list:
    key = _key("GET", url, params, None)
    path = _cache_path(key, "json")
    cached = _read_cache(path, ttl)
    if cached is not None:
        return json.loads(cached)
    resp = _request("GET", url, params=params, service=service, timeout=timeout)
    try:
        payload = resp.json()
    except ValueError as exc:
        raise ApiError(f"{service} 응답을 해석하지 못했습니다.") from exc
    _write_cache(path, json.dumps(payload, ensure_ascii=False))
    return payload


def get_text(url: str, params: dict | None = None, *, ttl: int = DEFAULT_TTL,
             service: str = "외부 서비스", timeout: int = 30) -> str:
    """ttl=0이면 캐시를 읽지도 쓰지도 않는다 (저작권이 있는 본문 등)."""
    key = _key("GET", url, params, None)
    path = _cache_path(key, "txt")
    if ttl > 0:
        cached = _read_cache(path, ttl)
        if cached is not None:
            return cached
    resp = _request("GET", url, params=params, service=service, timeout=timeout)
    if ttl > 0:
        _write_cache(path, resp.text)
    return resp.text


def post_text(url: str, body: str, *, ttl: int = 7 * DEFAULT_TTL, service: str = "외부 서비스",
              timeout: int = 120, retries: int = 1) -> str:
    key = _key("POST", url, None, body)
    path = _cache_path(key, "txt")
    cached = _read_cache(path, ttl)
    if cached is not None:
        return cached
    resp = _request("POST", url, data=body, headers={"Content-Type": "text/plain"},
                    service=service, timeout=timeout, retries=retries)
    _write_cache(path, resp.text)
    return resp.text

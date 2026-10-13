"""PepLab 휴대용 실행기: 빈 포트를 찾아 앱 서버를 켜고, 준비되면 브라우저를 연다."""
import socket
import subprocess
import sys
import time
import urllib.request
import webbrowser
from pathlib import Path

HERE = Path(__file__).resolve().parent
WAIT_SECONDS = 240


def port_free(port):
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(("127.0.0.1", port)) != 0


def healthy(port):
    try:
        return urllib.request.urlopen(f"http://127.0.0.1:{port}/_stcore/health", timeout=2).read() == b"ok"
    except Exception:
        return False


def main():
    port = next((p for p in range(8501, 8521) if port_free(p)), None)
    if port is None:
        print("[오류] 8501~8520 포트를 모두 다른 프로그램이 쓰고 있습니다. 열려 있는 PepLab 창을 닫고 다시 실행하세요.")
        return 1
    url = f"http://127.0.0.1:{port}"
    print(f"PepLab 서버를 켜는 중입니다... (처음 켤 때는 1~2분 걸릴 수 있습니다)\n서버가 준비되면 브라우저가 자동으로 열립니다. 주소: {url}")
    proc = subprocess.Popen(
        [sys.executable, "-m", "streamlit", "run", "app.py", "--server.port", str(port), "--server.address", "127.0.0.1",
         "--server.headless", "true", "--browser.gatherUsageStats", "false"], cwd=str(HERE))
    start = time.time()
    opened = False
    try:
        while proc.poll() is None:
            if not opened and healthy(port):
                webbrowser.open(url)
                opened = True
                print(f"\n브라우저가 열렸습니다. 안 보이면 주소창에 {url} 을 입력하세요.\n끝낼 때는 이 검은 창을 닫으세요.")
            if not opened and time.time() - start > WAIT_SECONDS:
                print(f"\n[오류] {WAIT_SECONDS}초가 지나도 서버가 준비되지 않았습니다. 위에 나온 오류 메시지를 확인하세요.")
                proc.terminate()
                return 1
            time.sleep(1)
    except KeyboardInterrupt:
        proc.terminate()
        return 0
    print(f"\n[오류] 서버가 종료되었습니다(코드 {proc.returncode}). 위에 나온 오류 메시지를 확인하세요.")
    return 1


if __name__ == "__main__":
    sys.exit(main())

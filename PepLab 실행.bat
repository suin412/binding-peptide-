@echo off
chcp 65001 >nul
setlocal
set "ROOT=%~dp0"
if not exist "%ROOT%runtime\python.exe" goto missing
if not exist "%ROOT%app\launch_peplab.py" goto missing
set "PYTHONHOME="
set "PYTHONPATH="
set "PYTHONNOUSERSITE=1"
set "MPLBACKEND=Agg"
set "PYTHONUTF8=1"
cd /d "%ROOT%app"
"%ROOT%runtime\python.exe" "%ROOT%app\launch_peplab.py"
if errorlevel 1 goto failed
goto end

:missing
echo.
echo [오류] 이 파일은 PepLab_portable 폴더 안에서만 실행됩니다.
echo        이 파일만 따로 복사하거나 다운로드하면 동작하지 않습니다.
echo        zip 파일 전체를 먼저 "압축 풀기" 한 뒤, 풀린 폴더 안의 이 파일을 실행하세요.
echo        (찾는 위치: %ROOT%runtime  /  %ROOT%app)
echo.
pause
goto end

:failed
echo.
echo 실행에 실패했습니다. 위의 오류 메시지를 확인하세요.
pause

:end

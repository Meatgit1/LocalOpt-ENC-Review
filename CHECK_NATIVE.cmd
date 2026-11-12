@echo off
setlocal EnableExtensions DisableDelayedExpansion
set "NATIVE_SETUP=%~dp0workflow\check_native.ps1"
if not exist "%NATIVE_SETUP%" (
  echo Keep CHECK_NATIVE.cmd beside the extracted workflow folder.
  pause
  exit /b 1
)
echo One-time native setup. No search, local solve or trajectory replay.
echo RemoteSigned applies only to this PowerShell process.
powershell.exe -NoLogo -NoProfile -NonInteractive -ExecutionPolicy RemoteSigned -File "%NATIVE_SETUP%" -BundleRoot "%~dp0."
set "SETUP_STATUS=%ERRORLEVEL%"
echo Exit code: %SETUP_STATUS%
pause
exit /b %SETUP_STATUS%

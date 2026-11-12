@echo off
setlocal EnableExtensions DisableDelayedExpansion
set "THESIS_BUNDLE_ROOT=%~dp0"
echo Start START_MAC.command first and wait for its Windows-worker message.
echo This window starts only the shared run selected by the Mac launcher.
powershell.exe -NoLogo -NoProfile -NonInteractive -Command "$ErrorActionPreference='Stop'; try { $base=$env:THESIS_BUNDLE_ROOT; $pointer=Join-Path $base 'settings/shared_worker.txt'; if (-not (Test-Path -LiteralPath $pointer -PathType Leaf)) { throw 'No shared run selected. Start the Mac launcher first.' }; $name=([IO.File]::ReadAllText($pointer)).Trim(); if ($name -cnotmatch '\A[A-Za-z0-9][A-Za-z0-9_-]{0,63}\z') { throw 'Invalid shared-run name' }; $worker=Join-Path $base ('runs/shared/'+$name+'/native/RUN_NATIVE_WORKER.cmd'); if (-not (Test-Path -LiteralPath $worker -PathType Leaf)) { throw 'Worker is not ready. Wait for its path in the Mac window.' }; & $worker; exit $LASTEXITCODE } catch { [Console]::Error.WriteLine($_.Exception.Message); exit 1 }"
set "BUNDLE_STATUS=%ERRORLEVEL%"
echo Exit code: %BUNDLE_STATUS%
pause
exit /b %BUNDLE_STATUS%

param([Parameter(Mandatory=$true)][string]$BundleRoot)
$ErrorActionPreference = 'Stop'
Set-StrictMode -Version 2.0

# This records the installed EXE's observed interface, not its build history.
function Check-Capabilities($value) {
    $keys = @('schema','modes','profiles','observer_schema','native_build','pointer_bits','fresh_process_per_stage')
    if ($null -eq $value -or
        (@($value.PSObject.Properties.Name | Sort-Object) -join ',') -cne (($keys | Sort-Object) -join ',')) {
        throw 'InTrance returned incomplete interface information.'
    }
    if ($value.schema -isnot [string] -or $value.schema -cne 'intrance_unified_v1' -or
        $value.native_build -isnot [string] -or $value.native_build -cne 'x86_pe32' -or
        $value.observer_schema -isnot [int] -or $value.observer_schema -ne 1 -or
        $value.pointer_bits -isnot [int] -or $value.pointer_bits -ne 32 -or
        $value.fresh_process_per_stage -isnot [bool] -or -not $value.fresh_process_per_stage) {
        throw 'This InTrance executable does not support the local workflow.'
    }
    $modes = @('search','export','schedule')
    $profiles = @('mercury_body_rendezvous_v1','mercury_body_rendezvous_precision_v3')
    if (@($value.modes).Count -ne $modes.Count -or @($value.profiles).Count -ne $profiles.Count) {
        throw 'InTrance supports different run modes.'
    }
    for ($i=0; $i -lt $modes.Count; $i++) {
        if ($value.modes[$i] -isnot [string] -or $value.modes[$i] -cne $modes[$i]) {
            throw 'InTrance supports different run modes.'
        }
    }
    for ($i=0; $i -lt $profiles.Count; $i++) {
        if ($value.profiles[$i] -isnot [string] -or $value.profiles[$i] -cne $profiles[$i]) {
            throw 'InTrance supports different replay settings.'
        }
    }
}

$process = $null
$started = $false
$temporary = $null
try {
    if ([Environment]::OSVersion.Platform -ne 'Win32NT' -or
        $PSVersionTable.PSVersion.Major -ne 5 -or $PSVersionTable.PSVersion -lt [version]'5.1') {
        throw 'Run CHECK_NATIVE.cmd with Windows PowerShell 5.1.'
    }
    if ((Get-ExecutionPolicy -Scope Process) -ne 'RemoteSigned' -or
        (Get-ExecutionPolicy) -ne 'RemoteSigned') {
        throw 'This setup requires RemoteSigned for this PowerShell process only.'
    }
    $bundle = [IO.Path]::GetFullPath($BundleRoot).TrimEnd('\')
    $script = Join-Path $bundle 'workflow/check_native.ps1'
    if (-not [string]::Equals([IO.Path]::GetFullPath($PSCommandPath),
            [IO.Path]::GetFullPath($script),[StringComparison]::OrdinalIgnoreCase)) {
        throw 'Keep CHECK_NATIVE.cmd and workflow/ in the same extracted bundle.'
    }
    $executable = Join-Path $bundle 'workflow/runtime/InTrance.exe'
    if (-not (Test-Path -LiteralPath $executable -PathType Leaf) -or
        ((Get-Item -LiteralPath $executable -Force).Attributes -band [IO.FileAttributes]::ReparsePoint)) {
        throw 'The bundle is missing its InTrance executable.'
    }
    $before = (Get-FileHash -LiteralPath $executable -Algorithm SHA256).Hash.ToLowerInvariant()
    Write-Host 'Checking InTrance. No trajectory is calculated.'
    $start = [Diagnostics.ProcessStartInfo]::new()
    $start.FileName = $executable
    $start.Arguments = '--capabilities'
    $start.WorkingDirectory = Split-Path $executable -Parent
    $start.UseShellExecute = $false
    $start.CreateNoWindow = $true
    $start.RedirectStandardOutput = $true
    $start.RedirectStandardError = $true
    foreach ($name in @($start.EnvironmentVariables.Keys)) {
        if ($name.StartsWith('INTRANCE_',[StringComparison]::OrdinalIgnoreCase)) {
            $start.EnvironmentVariables.Remove($name)
        }
    }
    $process = [Diagnostics.Process]::new()
    $process.StartInfo = $start
    if (-not $process.Start()) { throw 'InTrance could not start.' }
    $started = $true
    $stdout = $process.StandardOutput.ReadToEndAsync()
    $stderr = $process.StandardError.ReadToEndAsync()
    if (-not $process.WaitForExit(10000)) {
        throw 'The native interface check exceeded 10 seconds.'
    }
    if (-not $stdout.Wait(1000) -or -not $stderr.Wait(1000)) {
        throw 'The native interface output did not close.'
    }
    $output = $stdout.Result.Trim()
    if ($process.ExitCode -ne 0 -or -not [string]::IsNullOrWhiteSpace($stderr.Result)) {
        throw 'InTrance could not report its interface.'
    }
    if ($output.Length -lt 2 -or $output.Length -gt 4096 -or $output[0] -cne '{' -or $output[$output.Length-1] -cne '}') {
        throw 'InTrance returned invalid interface information.'
    }
    # The native response is a flat JSON object with seven named fields.
    $names = @([regex]::Matches($output,'"([^"\\]+)"\s*:') | ForEach-Object { $_.Groups[1].Value })
    if ($names.Count -ne 7 -or @($names | Select-Object -Unique).Count -ne 7) {
        throw 'The native interface response contains repeated or extra fields.'
    }
    $capabilities = $output | ConvertFrom-Json -ErrorAction Stop
    Check-Capabilities $capabilities
    $after = (Get-FileHash -LiteralPath $executable -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($before -cne $after) { throw 'InTrance changed during setup.' }

    $directory = Join-Path $bundle '_internal/runtime'
    [IO.Directory]::CreateDirectory($directory) | Out-Null
    $destination = Join-Path $directory 'InTrance.capabilities.json'
    $record = [ordered]@{
        schema='intrance_native_capabilities_receipt_v1'
        executable_sha256=$after
        capabilities=$capabilities
    }
    if (Test-Path -LiteralPath $destination) {
        $item = Get-Item -LiteralPath $destination -Force
        if ($item.Length -gt 8192 -or ($item.Attributes -band [IO.FileAttributes]::ReparsePoint)) {
            throw 'The existing private native setup record is not valid.'
        }
        $existing = Get-Content -LiteralPath $destination -Raw | ConvertFrom-Json -ErrorAction Stop
        if ((@($existing.PSObject.Properties.Name | Sort-Object) -join ',') -cne 'capabilities,executable_sha256,schema' -or
            $existing.schema -cne $record.schema -or $existing.executable_sha256 -cne $after) {
            throw 'InTrance differs from this folder''s earlier setup. Use a new extracted folder.'
        }
        Check-Capabilities $existing.capabilities
    } else {
        $temporary = Join-Path $directory ([Guid]::NewGuid().ToString('N')+'.tmp')
        $stream = [IO.File]::Open($temporary,[IO.FileMode]::CreateNew,[IO.FileAccess]::Write,[IO.FileShare]::None)
        try {
            $bytes = [Text.UTF8Encoding]::new($false).GetBytes(($record | ConvertTo-Json -Depth 6)+"`n")
            $stream.Write($bytes,0,$bytes.Length)
            $stream.Flush()
        } finally { $stream.Dispose() }
        [IO.File]::Move($temporary,$destination)
        $temporary = $null
    }
    Write-Host 'Native interface checked.'
    Write-Host 'You can now check the settings and start a new run.'
    exit 0
} catch {
    [Console]::Error.WriteLine('STOP: '+$_.Exception.Message)
    exit 1
} finally {
    if ($null -ne $process) {
        if ($started -and -not $process.HasExited) {
            $process.Kill()
            $process.WaitForExit(2000) | Out-Null
        }
        $process.Dispose()
    }
    if ($null -ne $temporary -and (Test-Path -LiteralPath $temporary)) {
        [IO.File]::Delete($temporary)
    }
}

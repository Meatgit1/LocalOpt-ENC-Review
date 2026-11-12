param([switch]$Headless)
$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$root=$PSScriptRoot; $expectedConfig='@CONFIG_SHA@'
$stages=@('cold_search','export','precise_baseline','precise_candidate')
$clock=[Diagnostics.Stopwatch]::StartNew()
$sequence=0; $statusSequence=0; $session=$null; $challenge=$null
$state='WAITING'; $failure=$null; $termination=$null; $work=$null; $ownedChildCleanupVerified=$true
function Hash([string]$path) { (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant() }
function Under([string]$base,[string]$relative) {
  if ($relative -notmatch '^[A-Za-z0-9_./-]+$' -or [IO.Path]::IsPathRooted($relative) -or $relative.Split('/') -contains '..') { throw 'Invalid fixed relative path' }
  $full=[IO.Path]::GetFullPath((Join-Path $base $relative))
  if (-not $full.StartsWith([IO.Path]::GetFullPath($base).TrimEnd('\')+'\',[StringComparison]::OrdinalIgnoreCase)) { throw 'Path escape' }
  $part=$base
  foreach ($item in $relative.Split('/')) {
    $part=Join-Path $part $item
    if ((Test-Path -LiteralPath $part) -and (((Get-Item -LiteralPath $part -Force).Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0)) { throw 'Linked/reparse payload refused' }
  }
  return $full
}
function Publish([string]$relative,$value) {
  $to=Under $root $relative
  if (Test-Path -LiteralPath $to) { throw 'Publication already exists' }
  $tmp=$to+'.'+[Guid]::NewGuid().ToString('N')+'.tmp'
  [IO.File]::WriteAllText($tmp,($value | ConvertTo-Json -Depth 15 -Compress),[Text.UTF8Encoding]::new($false))
  [IO.File]::Move($tmp,$to)
}
function Exact-Keys($value,[string[]]$names) {
  if ((@($value.PSObject.Properties.Name | Sort-Object) -join ',') -cne (($names | Sort-Object) -join ',')) { throw 'Wrong JSON field set' }
}
function Heartbeat {
  $script:statusSequence++
  Publish ('_transport/status_{0:D6}.json' -f $statusSequence) ([ordered]@{schema=5;session_id=$session;state=$state;
    sequence=$sequence;challenge=$challenge;elapsed_seconds=$clock.Elapsed.TotalSeconds;failure=$failure;termination=$termination;
    owned_child_cleanup_verified=$ownedChildCleanupVerified})
}
# Progress text is display-only; receipts still determine process and scientific status.
$hostProgress=@{offset=[int64]0;pending=''}
function Read-ProgressLines([string]$path,$cursor) {
  if (-not [IO.File]::Exists($path)) { return }
  $stream=$null
  try {
    $stream=[IO.File]::Open($path,[IO.FileMode]::Open,[IO.FileAccess]::Read,[IO.FileShare]::ReadWrite)
    if ($stream.Length -lt $cursor.offset) { $cursor.offset=[int64]0; $cursor.pending='' }
    $stream.Position=$cursor.offset
    $buffer=New-Object byte[] 65536
    $count=$stream.Read($buffer,0,$buffer.Length); $cursor.offset=$stream.Position
    $text=$cursor.pending+[Text.Encoding]::UTF8.GetString($buffer,0,$count)
    $lines=$text.Split([char]10); $cursor.pending=$lines[-1]
    if ($cursor.pending.Length -gt 4096) { $cursor.pending='' }
    for ($i=0; $i -lt $lines.Length-1; $i++) {
      $line=($lines[$i] -replace '[^\x20-\x7e]',' ').Trim()
      if ($line.Length -gt 0 -and $line.Length -le 1024) { $line }
    }
  } catch [IO.IOException] {
    # The shared writer may briefly hold the file; the next poll tries again.
  } finally { if ($null -ne $stream) { $stream.Dispose() } }
}
function Show-HostProgress {
  foreach ($line in @(Read-ProgressLines (Join-Path (Split-Path $root -Parent) 'progress.log') $hostProgress)) {
    Write-Host $line
  }
}
function Show-NativeProgress([string]$path,$cursor,[bool]$errorStream=$false) {
  foreach ($line in @(Read-ProgressLines $path $cursor)) {
    if ($line -match '^EXIT_TRACE' -or $line -match '^NATIVE_(BRIDGE|SCHEDULE)_(COMPLETE|PROFILE_COMPLETE|REACHED_END)') { continue }
    if ($errorStream -or $line -match '^(NATIVE_BRIDGE:|input:|NC:|SC:|SIM:|EA:|\.\.\.|[+-][0-9]|CSV\s*:|DAT\s*:|GESOP\s*:|VRML\s*:|ANGLS\s*:|total running time:|InTrance run finished!)') {
      Write-Host $line
      try { [IO.File]::AppendAllText((Join-Path $root '_transport/native_console.log'),($line+"`n"),[Text.UTF8Encoding]::new($false)) }
      catch [IO.IOException] { } # A missed display update does not alter the saved native streams.
    }
  }
}
function Check-Budget {
  if ($clock.Elapsed.TotalSeconds -ge 7200) { throw 'Session deadline reached' }
  $stop=Under $root '_transport/stop.json'
  if (Test-Path -LiteralPath $stop) {
    $raw=[IO.File]::ReadAllText($stop)
    $pattern='\A\{"schema":5,"session_id":"(?<session>[a-f0-9]{32})","reason":"(?<reason>[a-z][a-z0-9_]{0,63})"\}\z'
    if ($raw -cnotmatch $pattern) { throw 'Invalid stop record' }
    if ($Matches.session -cne $session) { throw 'Invalid stop nonce' }
    $script:termination=[ordered]@{kind='controlled_host_stop';code=$Matches.reason}
    throw 'HOST_CANCELLED'
  }
}
function Check-Config {
  if ((Hash (Under $root 'transport.json')) -cne $expectedConfig -or (Hash (Under $root 'source.json')) -cne $config.source_sha256) { throw 'Transport configuration/provenance changed' }
  foreach ($item in $config.assets.PSObject.Properties) {
    Check-Budget
    if ((Hash (Under $root $item.Name)) -cne $item.Value) { throw ('Frozen asset changed: '+$item.Name) }
  }
  if (Test-Path -LiteralPath (Under $root 'assets/reference/mercury.eac')) { throw 'Historical controller in fresh inputs' }
}
function Runtime-Preflight {
  # Read-only runtime checks; the new private directory and session publication prove write access.
  $tokens=$null; $errors=$null
  [System.Management.Automation.Language.Parser]::ParseFile($PSCommandPath,[ref]$tokens,[ref]$errors) | Out-Null
  if ($errors.Count -ne 0) { throw 'Native worker syntax preflight failed' }
  if (((Get-Item -LiteralPath $root -Force).Attributes -band [IO.FileAttributes]::ReparsePoint) -ne 0) { throw 'Linked shared workspace root refused' }
  $machines=[ordered]@{}
  foreach ($role in @('search','export','precision')) {
    $path=Under $root ('assets/'+$role+'.exe')
    $stream=[IO.File]::Open($path,[IO.FileMode]::Open,[IO.FileAccess]::Read,[IO.FileShare]::Read)
    $reader=[IO.BinaryReader]::new($stream)
    try {
      if ($reader.ReadUInt16() -ne 0x5a4d -or $stream.Length -lt 128) { throw 'Native executable is not a PE image' }
      $stream.Position=0x3c; $offset=$reader.ReadInt32()
      if ($offset -lt 64 -or $offset -gt ($stream.Length - 26)) { throw 'Invalid PE header offset' }
      $stream.Position=$offset
      if ($reader.ReadUInt32() -ne 0x00004550) { throw 'Invalid PE signature' }
      $machine=$reader.ReadUInt16(); $stream.Position=$offset+24; $magic=$reader.ReadUInt16()
      if ($machine -ne 0x014c -or $magic -ne 0x010b) { throw 'This release requires the pinned x86 PE32 native executables' }
      $machines[$role]='x86_pe32'
    } finally { $reader.Dispose(); $stream.Dispose() }
  }
  $assetBytes=[int64]0
  foreach ($name in $config.assets.PSObject.Properties.Name) { $assetBytes+=(Get-Item -LiteralPath (Under $root $name)).Length }
  $required=[int64](4*$assetBytes+64MB)
  $drive=[IO.DriveInfo]::new([IO.Path]::GetPathRoot($work))
  if ($drive.AvailableFreeSpace -lt $required) { throw ('Insufficient local staging space: require '+$required+' bytes') }
  if (($work.Length + 80) -ge 240) { throw 'Private staging path is too long for native filenames' }
  Publish '_transport/preflight.json' ([ordered]@{schema=5;run_id=$config.run_id;session_id=$session;
    powershell_version=$PSVersionTable.PSVersion.ToString();os_version=[Environment]::OSVersion.VersionString;
    shared_directory=$root;local_directory=$work;asset_bytes=$assetBytes;required_local_bytes=$required;
    available_local_bytes=$drive.AvailableFreeSpace;native_images=$machines;syntax_checked=$true;
    shared_session_published=$true;ready=$true})
}
function Check-Lineage([string]$name,[string]$base,$manifest) {
  $path=Under $base '_audit/source.json'
  if ((Hash $path) -cne $manifest.source_sha256) { throw 'Stage source record changed' }
  $source=Get-Content -LiteralPath $path -Raw | ConvertFrom-Json
  if ($source.run_id -cne $config.run_id -or $source.seed -ne $config.search.seed) { throw 'Wrong stage source run/seed' }
  if ($name -ceq 'cold_search') {
    Exact-Keys $source @('schema','run_id','kind','controller_input','seed','evaluation_limit')
    if ($source.kind -cne 'fresh_cold_search' -or $null -ne $source.controller_input -or $source.evaluation_limit -ne $config.search.evaluation_limit) { throw 'Fresh cold-search source differs' }
    if (Test-Path -LiteralPath (Under $base 'fixture/mercury.eac')) { throw 'Fresh search cannot receive a controller' }
  } else {
    Exact-Keys $source @('schema','run_id','kind','controller_source','controller_sha256','search_receipt_sha256','observer_summary_sha256','observer_metadata_sha256','session_id','seed','completed_evaluations')
    if ($source.kind -cne 'same_run_observer_final' -or $source.controller_source -cne 'stages/cold_search/result/enc_observer/final.eac' -or
        $source.session_id -cne $session -or $source.completed_evaluations -lt 40 -or
        $source.completed_evaluations -gt $config.search.evaluation_limit -or ($source.completed_evaluations % 2) -ne 0 -or
        ($config.search.stop_policy -ceq 'exact_completed_count' -and $source.completed_evaluations -ne $config.search.evaluation_limit)) { throw 'Not this session fresh controller' }
    if ((Hash (Under $root $source.controller_source)) -cne $source.controller_sha256 -or
        (Hash (Under $base 'fixture/mercury.eac')) -cne $source.controller_sha256 -or
        (Hash (Under $root 'stages/cold_search/_audit/receipt.json')) -cne $source.search_receipt_sha256 -or
        (Hash (Under $root 'stages/cold_search/result/enc_observer/summary.json')) -cne $source.observer_summary_sha256 -or
        (Hash (Under $root 'stages/cold_search/result/enc_observer/final.json')) -cne $source.observer_metadata_sha256) { throw 'Fresh controller lineage changed' }
    $summary=Get-Content -LiteralPath (Under $root 'stages/cold_search/result/enc_observer/summary.json') -Raw | ConvertFrom-Json
    if ($summary.completed_evaluations -ne $source.completed_evaluations -or $summary.seed -ne $source.seed -or
        $summary.evaluation_limit -ne $config.search.evaluation_limit -or $summary.status -cne 'complete') { throw 'Attained native count lineage changed' }
    $searchReceipt=Get-Content -LiteralPath (Under $root 'stages/cold_search/_audit/receipt.json') -Raw | ConvertFrom-Json
    if ($searchReceipt.session_id -cne $session -or $searchReceipt.stage_id -cne 'cold_search' -or $searchReceipt.success -ne $true -or
        $searchReceipt.owned_child_cleanup_verified -ne $true -or $searchReceipt.exit_code -ne 0) { throw 'Fresh search receipt not accepted' }
  }
}
function Check-Stage([string]$name,[string]$manifestHash) {
  $base=Under $root ('stages/'+$name); $path=Under $base '_audit/manifest.json'
  if ((Hash $path) -cne $manifestHash) { throw 'Stage manifest changed' }
  $m=Get-Content -LiteralPath $path -Raw | ConvertFrom-Json
  Exact-Keys $m @('schema','stage_id','mode','profile','tolerance','command','executable_sha256','child_seconds','run_id','config_sha256','source_sha256','hashes')
  $search=$name -ceq 'cold_search'; $precise=$name.StartsWith('precise_')
  $mode='export'; $profile='mercury_body_rendezvous_v1'; $tol='1e-6'; $command='evaluate NC'; $exeHash=$config.executables.export; $deadline=90
  if ($search) { $mode='cold_search'; $profile='enc_count_observer_v1'; $command='optimize NC'; $exeHash=$config.executables.search; $deadline=$config.search.child_seconds }
  elseif ($precise) { $mode='schedule'; $profile='mercury_body_rendezvous_precision_v3'; $tol='1e-14'; $exeHash=$config.executables.precision }
  if ($m.schema -ne 5 -or $m.stage_id -cne $name -or $m.mode -cne $mode -or $m.profile -cne $profile -or $m.tolerance -cne $tol -or
      $m.command -cne $command -or $m.executable_sha256 -cne $exeHash -or $m.child_seconds -ne $deadline -or $m.run_id -cne $config.run_id -or
      $m.config_sha256 -cne $expectedConfig) { throw 'Wrong fresh stage contract' }
  $fixed=$config.fixtures.$name
  $fixtureNames=@($fixed.PSObject.Properties.Name)
  if (-not $search) { $fixtureNames+=@('mercury.eac') }
  if ($precise) { $fixtureNames+=$scheduleNames }
  $names=@('InTrance.exe')+@($fixtureNames | ForEach-Object { 'fixture/'+$_ })+@($dataNames | ForEach-Object { 'data/'+$_ })
  Exact-Keys $m.hashes $names
  $actual=@(Get-ChildItem -LiteralPath (Under $base 'fixture') -Force | ForEach-Object { $_.Name } | Sort-Object)
  if (($actual -join ',') -cne (($fixtureNames | Sort-Object) -join ',')) { throw 'Unexpected fixture file' }
  foreach ($item in $m.hashes.PSObject.Properties) {
    Check-Budget
    if ((Hash (Under $base $item.Name)) -cne $item.Value) { throw ('Payload changed: '+$item.Name) }
  }
  if ($m.hashes.'InTrance.exe' -cne $exeHash) { throw 'Executable changed' }
  foreach ($n in $fixed.PSObject.Properties.Name) { if ($m.hashes.('fixture/'+$n) -cne $fixed.$n) { throw 'Fixed mission changed' } }
  foreach ($n in $dataNames) { if ($m.hashes.('data/'+$n) -cne $config.assets.('assets/data/'+$n)) { throw 'Fixed data changed' } }
  Check-Lineage $name $base $m
  return $m
}
function Payload-Hashes([string]$base,$manifest) {
  $values=[ordered]@{}
  foreach ($item in $manifest.hashes.PSObject.Properties) {
    $name=$item.Name.Replace('fixture/','result/')
    $values[$item.Name]=Hash (Under $base $name)
    if ($values[$item.Name] -cne $item.Value) { throw ('Guest input changed: '+$name) }
  }
  return $values
}
function Run-Native($request) {
  $name=$request.stage_id; $base=Under $root ('stages/'+$name); $m=Check-Stage $name $request.manifest_sha256
  if ((Test-Path -LiteralPath (Join-Path $base 'attempt.claim')) -or (Test-Path -LiteralPath (Join-Path $base 'result'))) { throw 'Native stage already attempted' }
  $claim=[IO.File]::Open((Join-Path $base 'attempt.claim'),[IO.FileMode]::CreateNew); $claim.Dispose()
  $guest=Join-Path $work $name; [IO.Directory]::CreateDirectory($guest) | Out-Null
  $p=$null; $out=$null; $err=$null; $outTask=$null; $errTask=$null; $started=$false; $childClock=$null
  $outCursor=@{offset=[int64]0;pending=''}; $errCursor=@{offset=[int64]0;pending=''}
  $record=[ordered]@{schema=5;session_id=$session;sequence=$sequence;challenge=$request.challenge;stage_id=$name;manifest_sha256=$request.manifest_sha256;
    success=$false;failure=$null;exit_code=$null;timed_out=$false;process_wall_seconds=$null;cpu_seconds=$null;pid=$null;
    inputs_before=$null;inputs_after=$null;output_hashes=$null;guest_directory=$guest;owned_child_cleanup_verified=$false;timing_source='Stopwatch';termination=$null}
  try {
    Copy-Item -LiteralPath (Under $base 'InTrance.exe') -Destination $guest
    Copy-Item -LiteralPath (Under $base 'data') -Destination $guest -Recurse
    Copy-Item -LiteralPath (Under $base 'fixture') -Destination (Join-Path $guest 'result') -Recurse
    $record.inputs_before=Payload-Hashes $guest $m
    $si=New-Object Diagnostics.ProcessStartInfo
    $si.FileName=Join-Path $guest 'InTrance.exe'; $si.Arguments='mercury.inp'; $si.WorkingDirectory=Join-Path $guest 'result'
    $si.UseShellExecute=$false; $si.CreateNoWindow=$true; $si.RedirectStandardOutput=$true; $si.RedirectStandardError=$true
    foreach ($key in @($si.EnvironmentVariables.Keys)) { if ([string]$key -like 'INTRANCE_*') { $si.EnvironmentVariables.Remove([string]$key) } }
    if ($name -ceq 'cold_search') {
      [IO.Directory]::CreateDirectory((Join-Path $guest 'result/enc_observer')) | Out-Null
      $si.EnvironmentVariables['INTRANCE_ENC_OBSERVER_CONFIG']='enc_observer.cfg'
    } else {
      $si.EnvironmentVariables['INTRANCE_BRIDGE_MODE']=$m.mode; $si.EnvironmentVariables['INTRANCE_BRIDGE_PROFILE']=$m.profile
    }
    $out=[IO.FileStream]::new((Join-Path $guest 'result/stdout.txt'),[IO.FileMode]::CreateNew,[IO.FileAccess]::Write,[IO.FileShare]::Read,1,$true)
    $err=[IO.FileStream]::new((Join-Path $guest 'result/stderr.txt'),[IO.FileMode]::CreateNew,[IO.FileAccess]::Write,[IO.FileShare]::Read,1,$true)
    Check-Budget
    $p=New-Object Diagnostics.Process; $p.StartInfo=$si
    $childClock=[Diagnostics.Stopwatch]::StartNew()
    if (-not $p.Start()) { throw 'Native start failed' }; $started=$true; $handle=$p.Handle; $record.pid=$p.Id
    $outTask=$p.StandardOutput.BaseStream.CopyToAsync($out); $errTask=$p.StandardError.BaseStream.CopyToAsync($err)
    $titles=@{cold_search='[1/4] ENC search';export='[2/4] Export states and sail commands';precise_baseline='[3/4] Replay starting schedule';precise_candidate='[4/4] Replay optimised schedule'}
    Write-Host $titles[$name]
    if ($name -ceq 'cold_search') { Write-Host 'Native search output: fitness, relative improvement and native t/r/v terms. Final arrival errors are reported in km and m/s after replay.' }
    while (-not $p.HasExited) {
      if ($childClock.Elapsed.TotalSeconds -ge $m.child_seconds) { $record.timed_out=$true; throw 'Native child deadline reached; no partial handover' }
      Show-NativeProgress (Join-Path $guest 'result/stdout.txt') $outCursor
      Show-NativeProgress (Join-Path $guest 'result/stderr.txt') $errCursor $true
      Check-Budget; Start-Sleep -Milliseconds 100
    }
    $childClock.Stop(); $record.process_wall_seconds=$childClock.Elapsed.TotalSeconds
    if (-not $outTask.Wait(5000) -or -not $errTask.Wait(5000)) { throw 'Native stream EOF absent' }
    $out.Flush(); $err.Flush()
    Show-NativeProgress (Join-Path $guest 'result/stdout.txt') $outCursor
    Show-NativeProgress (Join-Path $guest 'result/stderr.txt') $errCursor $true
    Write-Host ('InTrance calculation finished in {0:F2} s.' -f $childClock.Elapsed.TotalSeconds)
    $record.exit_code=$p.ExitCode; $record.cpu_seconds=$p.TotalProcessorTime.TotalSeconds
    if ($record.exit_code -ne 0 -or $record.process_wall_seconds -gt $m.child_seconds) { throw 'Native exit/deadline failed' }
    $record.inputs_after=Payload-Hashes $guest $m
    Check-Stage $name $request.manifest_sha256 | Out-Null
    [IO.File]::WriteAllText((Join-Path $guest 'result/native_exit_code.txt'),"0`n")
    [IO.File]::WriteAllText((Join-Path $guest 'result/guest_identity_verified.txt'),"Fresh workflow inputs verified before and after native child.`n")
    $record.success=$true
  } catch {
    $record.failure=$_.Exception.Message
    if ($record.failure -ceq 'HOST_CANCELLED' -and $null -ne $termination) { $record.termination=$termination }
    else {
      $code='native_stage_error'
      if ($record.timed_out) { $code='native_child_deadline' }
      elseif ($record.failure -ceq 'Native stream EOF absent') { $code='native_stream_eof_missing' }
      $record.termination=[ordered]@{kind='native_or_transport_error';code=$code}
    }
  }
  finally {
    $cleanupOK=$true
    if ($started) {
      try {
        if (-not $p.HasExited) { $p.Kill(); if (-not $p.WaitForExit(5000)) { throw 'Owned child did not exit' } }
        if ($childClock.IsRunning) { $childClock.Stop() }
        $record.exit_code=$p.ExitCode; $record.process_wall_seconds=$childClock.Elapsed.TotalSeconds; $record.cpu_seconds=$p.TotalProcessorTime.TotalSeconds
        if ($null -ne $outTask -and -not $outTask.Wait(5000)) { throw 'stdout did not drain' }
        if ($null -ne $errTask -and -not $errTask.Wait(5000)) { throw 'stderr did not drain' }
      } catch { $cleanupOK=$false; $record.success=$false; $record.failure=([string]$record.failure+'; cleanup: '+$_.Exception.Message) }
    }
    try {
      if ($null -ne $out) { $out.Dispose() }; if ($null -ne $err) { $err.Dispose() }; if ($null -ne $p) { $p.Dispose() }
    } catch { $cleanupOK=$false; $record.success=$false; $record.failure=([string]$record.failure+'; dispose: '+$_.Exception.Message) }
    if (-not $cleanupOK) { $record.termination=[ordered]@{kind='native_or_transport_error';code='owned_cleanup_error'} }
    $record.owned_child_cleanup_verified=$cleanupOK; $script:ownedChildCleanupVerified=$cleanupOK
    try {
      $result=Join-Path $guest 'result'; $destination=Join-Path $base 'result'
      if (Test-Path -LiteralPath $destination) { throw 'Copy-back result exists' }
      if (-not (Test-Path -LiteralPath $result)) { throw 'Guest evidence absent' }
      $publication=Join-Path $base ('result_publish_'+[Guid]::NewGuid().ToString('N'))
      Copy-Item -LiteralPath $result -Destination $publication -Recurse
      [IO.Directory]::Move($publication,$destination)
      $hashes=[ordered]@{}
      foreach ($file in @(Get-ChildItem -LiteralPath $destination -File -Recurse)) {
        $relative=$file.FullName.Substring($destination.TrimEnd('\').Length+1).Replace('\','/')
        $hashes[$relative]=Hash (Under $destination $relative)
      }
      $record.output_hashes=$hashes
    } catch {
      $record.success=$false; $record.failure=([string]$record.failure+'; collection: '+$_.Exception.Message)
      $record.termination=[ordered]@{kind='native_or_transport_error';code='collection_error'}
    }
    Publish ('stages/'+$name+'/_audit/receipt.json') $record
  }
  if (-not $record.success) {
    if ($null -ne $record.termination -and $record.termination.kind -ceq 'controlled_host_stop' -and $record.owned_child_cleanup_verified) { throw 'HOST_CANCELLED' }
    $script:termination=$record.termination
    throw ('Native stage failed: '+$record.failure)
  }
}
try {
  if ([Environment]::OSVersion.Platform -ne 'Win32NT' -or $PSVersionTable.PSVersion.Major -ne 5 -or $PSVersionTable.PSVersion -lt [version]'5.1') { throw 'Windows PowerShell 5.1 required' }
  if ((Get-ExecutionPolicy -Scope MachinePolicy) -ne 'Undefined' -or (Get-ExecutionPolicy -Scope UserPolicy) -ne 'Undefined' -or (Get-ExecutionPolicy) -ne 'RemoteSigned') { throw 'Only process-scoped RemoteSigned permitted; managed policies are not overridden' }
  if ((Hash (Under $root 'transport.json')) -cne $expectedConfig) { throw 'Transport configuration changed' }
  $config=Get-Content -LiteralPath (Under $root 'transport.json') -Raw | ConvertFrom-Json
  if ($config.schema -ne 5 -or $config.kind -cne 'fresh_enc_workflow_transport_v3' -or ($config.stage_ids -join ',') -cne ($stages -join ',') -or
      $config.session_seconds -ne 7200 -or $config.child_seconds -ne 90 -or $config.max_jobs -ne 4 -or $config.search.warm_start -ne $false -or
      $config.search.stop_policy -cne 'exact_completed_count' -or $config.search.deadline_policy -cne 'hard_timeout_reject_not_orderly_cutoff' -or
      $config.search.child_seconds -lt 1 -or $config.search.child_seconds -gt 1200 -or $config.search.evaluation_limit -lt 40 -or
      $config.search.evaluation_limit -gt 1000000 -or ($config.search.evaluation_limit % 2) -ne 0) { throw 'Wrong fresh worker contract' }
  $dataNames=@('bin1600_600.405','ELEMENTS.UNNUM','ELEMENTS.NUMBR','ELEMENTS.COMET')
  $scheduleNames=@('bridge_input.csv','bridge_initial_input.csv','bridge_schedule_limits.csv')
  $claim=[IO.File]::Open((Under $root '_transport/worker.claim'),[IO.FileMode]::CreateNew); $claim.Dispose()
  $session=[Guid]::NewGuid().ToString('N'); $challenge=[Guid]::NewGuid().ToString('N')
  $work=Join-Path ([Environment]::GetFolderPath('LocalApplicationData')) ('InTrance/mercury_'+$session)
  if (Test-Path -LiteralPath $work) { throw 'Guest workspace already exists' }; [IO.Directory]::CreateDirectory($work) | Out-Null
  Publish '_transport/session.json' ([ordered]@{schema=5;session_id=$session;run_id=$config.run_id;config_sha256=$expectedConfig;worker_pid=$PID;session_seconds=7200})
  Check-Config; Runtime-Preflight
  Write-Host 'EARTH -> MERCURY | IDEAL SOLAR SAIL'
  Write-Host ('ENC seed: '+$config.search.seed+' | evaluation limit: '+$config.search.evaluation_limit)
  Write-Host 'ENC search -> export -> local optimisation -> InTrance replay -> results'
  Write-Host 'Keep this window open. Local progress will appear here while the solver runs.'
  Show-HostProgress
  while ($sequence -lt 4) {
    $state='WAITING'; Heartbeat; $heartbeatAt=$clock.Elapsed.TotalSeconds
    if ($sequence -eq 3) { Write-Host 'LOCAL OPTIMISATION | progress below comes from the Python solver' }
    $requestPath=Under $root ('_transport/requests/'+$sequence+'.json')
    while ($true) {
      Show-HostProgress
      Check-Budget
      $pending=@(Get-ChildItem -LiteralPath (Under $root '_transport/requests') -File -Filter '*.json')
      if ($pending.Count -gt 1) { throw 'Multiple pending requests' }
      if ($pending.Count -eq 1) { if ($pending[0].Name -cne ($sequence.ToString()+'.json')) { throw 'Out-of-order request' }; break }
      if ($clock.Elapsed.TotalSeconds-$heartbeatAt -ge 10) { Heartbeat; $heartbeatAt=$clock.Elapsed.TotalSeconds }
      Start-Sleep -Milliseconds 200
    }
    Check-Budget
    if ((Get-Item -LiteralPath $requestPath).Length -gt 1024) { throw 'Oversized request' }
    $raw=[IO.File]::ReadAllText($requestPath)
    $pattern='\A\{"schema":5,"session_id":"(?<session>[a-f0-9]{32})","sequence":(?<sequence>[0-3]),"challenge":"(?<challenge>[a-f0-9]{32})","stage_id":"(?<stage>cold_search|export|precise_baseline|precise_candidate)","manifest_sha256":"(?<hash>[a-f0-9]{64})"\}\z'
    if ($raw -cnotmatch $pattern) { throw 'Invalid canonical request' }
    if ($Matches.session -cne $session -or [int]$Matches.sequence -ne $sequence -or $Matches.challenge -cne $challenge -or $Matches.stage -cne $stages[$sequence]) { throw 'Session/challenge/order mismatch' }
    $request=$raw | ConvertFrom-Json
    Check-Config; Check-Stage $request.stage_id $request.manifest_sha256 | Out-Null
    [IO.File]::Move($requestPath,(Under $root ('_transport/accepted/'+$sequence+'.json')))
    $state='RUNNING'; Heartbeat; Run-Native $request
    $sequence++; $challenge=[Guid]::NewGuid().ToString('N')
  }
  Check-Config; $state='COMPLETED'; $termination=[ordered]@{kind='completed';code='all_native_stages_completed'}; Heartbeat
  Write-Host 'Native replay complete. Checking arrival and numerical agreement, then writing plots.'
  if (-not $Headless) {
    $displayClock=[Diagnostics.Stopwatch]::StartNew()
    $done=Join-Path (Split-Path $root -Parent) 'progress.done'
    while ($displayClock.Elapsed.TotalSeconds -lt 120 -and $clock.Elapsed.TotalSeconds -lt 7200 -and -not (Test-Path -LiteralPath $done)) {
      Show-HostProgress; Start-Sleep -Milliseconds 200
    }
    Show-HostProgress
    if (-not (Test-Path -LiteralPath $done)) { Write-Host 'The host report is still pending. Read result_summary.json when the Python driver finishes.' }
  }
} catch {
  $failure=$_.Exception.Message; $state='FAILED'
  if ($failure -ceq 'HOST_CANCELLED' -and $null -ne $termination -and $termination.kind -ceq 'controlled_host_stop' -and $ownedChildCleanupVerified) {
    $state='CANCELLED'; $failure=$null
  } elseif ($null -eq $termination -or $termination.kind -cne 'native_or_transport_error') {
    $code='worker_error'
    if ($failure -ceq 'Session deadline reached') { $code='worker_session_deadline' }
    $termination=[ordered]@{kind='native_or_transport_error';code=$code}
  }
  if ($null -ne $session) { try { Heartbeat } catch { Write-Host 'Unable to publish terminal status; preserve guest evidence.' } }
  if ($state -ceq 'CANCELLED') {
    Write-Host ('Host requested a controlled stop: '+$termination.code)
    Write-Host 'Transport stopped cleanly. The scientific workflow did not complete or pass acceptance.'
  } else { Write-Host ('STOP: '+$failure) }
  if ($null -ne $work) { Write-Host ('Guest evidence: '+$work) }
  if ($state -ceq 'CANCELLED') { exit 0 }
  exit 1
}
exit 0

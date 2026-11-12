#ifndef INTRANCE_EXIT_DIAGNOSTICS_H
#define INTRANCE_EXIT_DIAGNOSTICS_H

// Diagnostic working copy only. A CRT assertion is logged and terminates this
// process with failure (134); it is never ignored or counted as success.
#if defined(_WIN32) && defined(_DEBUG)
#include <windows.h>
#include <dbghelp.h>
#include <crtdbg.h>
#include <cstdio>
#include <cstring>
#include <process.h>
#pragma comment(lib, "dbghelp.lib")

inline BOOL& InTranceSymbolsReady()
{
  static BOOL ready = FALSE;
  return ready;
}

inline int __cdecl InTranceHeapReport(int kind, char* message, int*)
{
  if (kind != _CRT_ASSERT && kind != _CRT_ERROR) return FALSE;
  static volatile LONG handling = 0;
  if (InterlockedExchange(&handling, 1) == 0)
  {
    std::fprintf(stderr, "\nINTRANCE_FATAL_CRT_REPORT type=%d\n%s\n", kind,
                 message ? message : "(no report text)");
    std::fflush(stderr);
    void* frames[48];
    const USHORT count = CaptureStackBackTrace(0, 48, frames, NULL);
    const HANDLE process = GetCurrentProcess();
    // Preserve raw evidence first; symbol lookup can itself allocate memory.
    for (USHORT i = 0; i < count; ++i)
      std::fprintf(stderr, "RAW_STACK %02u address=%p\n", static_cast<unsigned>(i), frames[i]);
    std::fflush(stderr);
    const BOOL symbols = InTranceSymbolsReady();
    for (USHORT i = 0; i < count; ++i)
    {
      std::fprintf(stderr, "STACK %02u address=%p", static_cast<unsigned>(i), frames[i]);
      if (symbols)
      {
        struct SymbolBuffer { SYMBOL_INFO info; char name[MAX_SYM_NAME]; } buffer = {};
        buffer.info.SizeOfStruct = sizeof(SYMBOL_INFO);
        buffer.info.MaxNameLen = MAX_SYM_NAME;
        DWORD64 offset = 0;
        const DWORD64 address = reinterpret_cast<DWORD64>(frames[i]);
        if (SymFromAddr(process, address, &offset, &buffer.info))
          std::fprintf(stderr, " %s+0x%llx", buffer.info.Name, static_cast<unsigned long long>(offset));
        IMAGEHLP_LINE64 line = {};
        line.SizeOfStruct = sizeof(line);
        DWORD lineOffset = 0;
        if (SymGetLineFromAddr64(process, address, &lineOffset, &line))
          std::fprintf(stderr, " %s:%lu", line.FileName, static_cast<unsigned long>(line.LineNumber));
      }
      std::fprintf(stderr, "\n");
      std::fflush(stderr);
    }
    std::fprintf(stderr, "INTRANCE_FATAL_CRT_REPORT: terminating with exit code 134\n");
    std::fflush(stderr);
  }
  // Do not attempt normal C++ cleanup again after a failed heap assertion.
  TerminateProcess(GetCurrentProcess(), 134);
  _exit(134);
  return TRUE;
}

inline void InstallInTranceExitDiagnostics()
{
  _CrtSetReportHook(InTranceHeapReport);
  char exeDirectory[MAX_PATH] = {};
  const DWORD length = GetModuleFileNameA(NULL, exeDirectory, MAX_PATH);
  if (length > 0 && length < MAX_PATH)
  {
    char* slash = std::strrchr(exeDirectory, '\\');
    if (slash) *slash = '\0';
  }
  else exeDirectory[0] = '\0';
  SymSetOptions(SymGetOptions() | SYMOPT_LOAD_LINES | SYMOPT_UNDNAME | SYMOPT_DEFERRED_LOADS);
  InTranceSymbolsReady() = SymInitialize(GetCurrentProcess(), exeDirectory[0] ? exeDirectory : NULL, TRUE);
}
#else
inline void InstallInTranceExitDiagnostics() {}
#endif
#endif

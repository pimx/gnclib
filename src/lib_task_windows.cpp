
// lib_task_windows.cpp -- реализация TaskOs для Windows:
// CreateThread + Sleep + MemoryBarrier.
#ifdef _WIN32

#include "lib_task.h"

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <mmsystem.h>
#pragma comment(lib, "winmm.lib")

// Обёртка точки входа потока
struct WinThreadArg
{
    void (*fn)(void*);
    void* arg;
};

static DWORD WINAPI winThreadEntry(LPVOID p)
{
    WinThreadArg* a = (WinThreadArg*)p;
    a->fn(a->arg);
    return 0;
}

bool TaskOs::threadStart(void (*fn)(void*), void* arg, long long* handle)
{
    // статический пул аргументов: без динамических аллокаций
    static WinThreadArg s_args[16];
    static int s_next = 0;
    if (s_next >= 16)
    {
        return false;
    }
    WinThreadArg* a = &s_args[s_next];
    s_next = s_next + 1;
    a->fn = fn;
    a->arg = arg;
    HANDLE th = CreateThread(0, 0, winThreadEntry, a, 0, 0);
    if (th == 0)
    {
        return false;
    }
    *handle = (long long)th;
    return true;
}

void TaskOs::threadJoin(long long handle)
{
    WaitForSingleObject((HANDLE)handle, INFINITE);
    CloseHandle((HANDLE)handle);
}

void TaskOs::sleepMs(int ms)
{
    Sleep((DWORD)ms);
}

void TaskOs::memoryBarrier()
{
    MemoryBarrier();
}

//----------------------------------------------------------------------------
// Примитивы реального масштаба времени (Windows)
//----------------------------------------------------------------------------
double TaskOs::wallNow()
{
    LARGE_INTEGER f;
    LARGE_INTEGER c;
    QueryPerformanceFrequency(&f);
    QueryPerformanceCounter(&c);
    return (double)c.QuadPart / (double)f.QuadPart;
}

void TaskOs::sleepSec(double s)
{
    if (s <= 0.0)
    {
        return;
    }
    Sleep((DWORD)(s * 1000.0));
}

void TaskOs::yieldSpin()
{
    YieldProcessor(); // доспин последних долей миллисекунды
}

// Системный таймер по умолчанию 15.6 мс: одиночный Sleep на такте в
// единицы миллисекунд обрушил бы темп в разы
void TaskOs::timerResolutionBegin()
{
    timeBeginPeriod(1);
}

void TaskOs::timerResolutionEnd()
{
    timeEndPeriod(1);
}

#endif // _WIN32

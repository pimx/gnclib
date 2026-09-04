
// lib_task_posix.cpp -- реализация TaskOs для POSIX (Linux):
// pthread + nanosleep + __sync_synchronize.
#ifndef _WIN32

#include "lib_task.h"

#include <pthread.h>
#include <time.h>
#include <sched.h>

// Обёртка точки входа потока
struct PosixThreadArg
{
    void (*fn)(void*);
    void* arg;
};

static void* posixThreadEntry(void* p)
{
    PosixThreadArg* a = (PosixThreadArg*)p;
    void (*fn)(void*) = a->fn;
    void* arg = a->arg;
    a->fn = 0; // подтверждение копирования (аргумент на стеке стартера)
    fn(arg);
    return 0;
}

bool TaskOs::threadStart(void (*fn)(void*), void* arg, long long* handle)
{
    // статический пул аргументов: без динамических аллокаций
    static PosixThreadArg s_args[16];
    static int s_next = 0;
    if (s_next >= 16)
    {
        return false;
    }
    PosixThreadArg* a = &s_args[s_next];
    s_next = s_next + 1;
    a->fn = fn;
    a->arg = arg;
    pthread_t th;
    if (pthread_create(&th, 0, posixThreadEntry, a) != 0)
    {
        return false;
    }
    *handle = (long long)th;
    return true;
}

void TaskOs::threadJoin(long long handle)
{
    pthread_join((pthread_t)handle, 0);
}

void TaskOs::sleepMs(int ms)
{
    struct timespec ts;
    ts.tv_sec = ms / 1000;
    ts.tv_nsec = (long)(ms % 1000) * 1000000L;
    nanosleep(&ts, 0);
}

void TaskOs::memoryBarrier()
{
    __sync_synchronize();
}


//----------------------------------------------------------------------------
// Примитивы реального масштаба времени (POSIX)
//----------------------------------------------------------------------------
double TaskOs::wallNow()
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

void TaskOs::sleepSec(double s)
{
    if (s <= 0.0)
    {
        return;
    }
    struct timespec ts;
    ts.tv_sec = (time_t)s;
    ts.tv_nsec = (long)((s - (double)ts.tv_sec) * 1e9);
    nanosleep(&ts, 0);
}

void TaskOs::yieldSpin()
{
    // короткий сон вместо sched_yield: в цикле точного доспина при
    // высоких множителях темпа (такт десятки микросекунд) sched_yield
    // на нагруженных ядрах дорог -- конвейер отстаёт; nanosleep ~20 мкс
    // держит точность и освобождает процессор
    struct timespec ts;
    ts.tv_sec = 0;
    ts.tv_nsec = 20000;
    nanosleep(&ts, 0);
}

// nanosleep точен сам по себе (~50 мкс), настройка таймера не нужна
void TaskOs::timerResolutionBegin()
{
}

void TaskOs::timerResolutionEnd()
{
}

#endif // !_WIN32

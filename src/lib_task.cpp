
// lib_task.cpp -- переносимая часть ModTaskBase: конфигурация, цикл
// потока, конверты, рассылка, затравка. Обмен -- UdpReader/UdpSender
// (lib_udp); ОС-примитивы потоков -- TaskOs (lib_task_windows.cpp /
// lib_task_posix.cpp).
#include "lib_task.h"

#include <stdio.h>

ModTaskBase::ModTaskBase()
    : m_tick(0),
      m_port(0),
      m_nDest(0),
      m_bindAny(false),
      m_timeMaster(false),
      m_t0(0.0),
      m_dt(1.0 / 400.0),
      m_phase(0.0),
      m_tBase(0.0),
      m_thread(0),
      m_run(false),
      m_failed(false),
      m_idleMs(0),
      m_timeouts(0),
      m_stopTick(-1),
      m_pacer(0)
{
    // m_dest: конструкторы TaskDest дают порт 0 / localhost
}

ModTaskBase::~ModTaskBase()
{
    stop();
}

void ModTaskBase::configure(int port, const TaskDest* dest, int nDest,
                            double t0, double dt, double phase)
{
    m_port = port;
    m_nDest = (nDest > TASK_MAX_DEST) ? TASK_MAX_DEST : nDest;
    for (int i = 0; i < m_nDest; ++i)
    {
        m_dest[i] = dest[i];
    }
    m_t0 = t0;
    m_dt = dt;
    m_phase = phase;
    m_tBase = t0;
    m_tick = 0;
}

void ModTaskBase::setBindAny(bool bindAny)
{
    m_bindAny = bindAny;
}

void ModTaskBase::setTimeMaster(bool master)
{
    m_timeMaster = master;
}

void ModTaskBase::setPacer(TaskPacer* pacer)
{
    m_pacer = pacer;
}

void ModTaskBase::setStopTick(int stopTick)
{
    m_stopTick = stopTick;
}

bool ModTaskBase::start()
{
    // приём: только loopback, таймаут 50 мс (цикл проверяет флаг останова)
    if (!m_reader.Open(m_port, !m_bindAny, 50))
    {
        printf("task %d: recv socket open failed\n", m_port);
        m_failed = true;
        return false;
    }
    if (!m_sender.Open())
    {
        printf("task %d: send socket open failed\n", m_port);
        m_reader.Close();
        m_failed = true;
        return false;
    }
    m_run = true;
    if (!TaskOs::threadStart(&ModTaskBase::threadEntry, this, &m_thread))
    {
        printf("task %d: thread start failed\n", m_port);
        m_run = false;
        m_failed = true;
        return false;
    }
    return true;
}

void ModTaskBase::stop()
{
    if (m_thread != 0)
    {
        m_run = false;
        TaskOs::threadJoin(m_thread);
        m_thread = 0;
        m_reader.Close();
        m_sender.Close();
    }
}

void ModTaskBase::threadEntry(void* self)
{
    ((ModTaskBase*)self)->runLoop();
}

// Цикл задачи: такт за тактом до останова или предела тактов
void ModTaskBase::runLoop()
{
    while (m_run)
    {
        cycle();
        if (!m_run)
        {
            break;
        }
        if (m_pacer != 0)
        {
            m_pacer->pace(tickTime());
        }
        if (m_stopTick >= 0 && m_tick >= m_stopTick)
        {
            break;
        }
        m_tick = m_tick + 1;
        m_tBase = m_tBase + m_dt; // накопление (тождественно стенду)
    }
    m_run = false;
}

// Приём одной датаграммы с проверкой конверта; false -- таймаут/мусор
bool ModTaskBase::recvMsg(MsgEnvelope* env, void* payload, int maxSize)
{
    unsigned char buf[sizeof(MsgEnvelope) + MSG_MAX_PAYLOAD];
    int n = m_reader.Read(buf, (int)sizeof(buf));
    if (n < (int)sizeof(MsgEnvelope))
    {
        // таймаут приёма: учёт простоя (только после первого такта) и
        // периодическая ретрансляция последнего выхода (~200 мс) --
        // устойчивость к порядку запуска процессов
        if (m_tick > 0)
        {
            m_idleMs = m_idleMs + 50;
        }
        m_timeouts = m_timeouts + 1;
        if ((m_timeouts % 4) == 0)
        {
            resendLast();
        }
        return false;
    }
    m_idleMs = 0;
    m_timeouts = 0;
    memcpy(env, buf, sizeof(MsgEnvelope));
    if (env->magic != (uint32_t)MSG_MAGIC)
    {
        return false;
    }
    // size -- беззнаковый (uint16): проверка на отрицательность не нужна
    if ((int)env->size > maxSize ||
        n != (int)sizeof(MsgEnvelope) + (int)env->size)
    {
        return false;
    }
    memcpy(payload, buf + sizeof(MsgEnvelope), env->size);
    return true;
}

// Рассылка полезной нагрузки всем адресатам
void ModTaskBase::sendToAll(const void* payload, int size, int tk, double t)
{
    unsigned char buf[sizeof(MsgEnvelope) + MSG_MAX_PAYLOAD];
    MsgEnvelope env;
    env.magic = (uint32_t)MSG_MAGIC;
    env.srcPort = (uint16_t)m_port;
    env.tick = (int32_t)tk;
    env.size = (uint16_t)size;
    env.time = t;
    memcpy(buf, &env, sizeof(env));
    memcpy(buf + sizeof(env), payload, size);
    for (int i = 0; i < m_nDest; ++i)
    {
        m_sender.Send(m_dest[i].host, m_dest[i].port, buf,
                      (int)sizeof(env) + size);
    }
}

// Подстройка времени задачи под сообщения партнёра (см. lib_task.h)
bool ModTaskBase::maybeResync(const MsgEnvelope& env, int tickOffset)
{
    int need = m_tick + tickOffset;
    int gap = env.tick - need;
    if (gap > -TASK_RESYNC_GAP && gap < TASK_RESYNC_GAP)
    {
        return false; // расхождение в норме
    }
    // принять время партнёра: такт -- по конверту; базовое время --
    // произведением такта на шаг (env.time содержит фазу ОТПРАВИТЕЛЯ,
    // получателю неизвестную; после подстройки битовая тождественность
    // накоплению всё равно нарушена -- событие исключительное)
    m_tick = env.tick - tickOffset;
    m_tBase = m_t0 + (double)m_tick * m_dt;
    printf("task %d: time resync to tick %d (t=%.3f s, gap %d ticks)\n", m_port,
           m_tick, tickTime(), gap);
    return true;
}

// Затравка цикла (tick = -1) от основного потока
bool taskSendPrime(gnc::UdpSender& sender, int destPort, int srcPort,
                   const void* payload, int size, double t)
{
    unsigned char buf[sizeof(MsgEnvelope) + MSG_MAX_PAYLOAD];
    if (size > MSG_MAX_PAYLOAD)
    {
        return false;
    }
    MsgEnvelope env;
    env.magic = (uint32_t)MSG_MAGIC;
    env.srcPort = (uint16_t)srcPort;
    env.tick = -1;
    env.size = (uint16_t)size;
    env.time = t;
    memcpy(buf, &env, sizeof(env));
    memcpy(buf + sizeof(env), payload, size);
    return sender.Send(0, destPort, buf, (int)sizeof(env) + size);
}

//----------------------------------------------------------------------------
// RtPacer -- платформенно-нейтральная логика темпа; все обращения к ОС
// идут через примитивы TaskOs (реализованы в lib_task_windows.cpp /
// lib_task_posix.cpp).
//----------------------------------------------------------------------------
RtPacer::RtPacer()
    : m_speed(1.0),
      m_t0(0.0),
      m_maxLag(0.020),
      m_slip(0.0),
      m_slipMax(0.0),
      m_slipN(0),
      m_started(false),
      m_anchored(false)
{
}

RtPacer::~RtPacer()
{
    if (m_started)
    {
        TaskOs::timerResolutionEnd(); // парный к Begin в start()
    }
}

void RtPacer::start(double speed)
{
    m_speed = (speed > 1e-6) ? speed : 1.0;
    if (!m_started)
    {
        TaskOs::timerResolutionBegin();
    }
    m_started = true;
    m_anchored = false;
    m_slip = 0.0;
    m_slipMax = 0.0;
    m_slipN = 0;
    m_t0 = TaskOs::wallNow();
}

void RtPacer::pace(double t)
{
    if (!m_anchored)
    {
        // якорь: текущий модельный момент соответствует "сейчас"
        m_t0 = TaskOs::wallNow() - t / m_speed;
        m_anchored = true;
    }
    double target = t / m_speed; // настенная цель от якоря, с
    // ПРОЩЕНИЕ КРУПНОГО ОТСТАВАНИЯ. Если такт вышел за реальное время
    // больше чем на m_maxLag, долг НЕ отрабатывается: якорь сдвигается
    // на величину просрочки, и следующий такт идёт в нормальном темпе.
    // Иначе цикл гнал бы такты свободным темпом до погашения долга, а
    // потребитель телеметрии получал бы пачку кадров с интервалом в
    // десятки микросекунд. Потерянное время не скрывается: оно
    // накапливается в slip() и печатается отчётом стенда.
    double late = (TaskOs::wallNow() - m_t0) - target;
    if (late > m_maxLag)
    {
        m_t0 = m_t0 + late;
        m_slip = m_slip + late;
        m_slipN = m_slipN + 1;
        if (late > m_slipMax)
        {
            m_slipMax = late;
        }
        return; // цель уже наступила, спать нечего
    }
    for (;;)
    {
        double ahead = target - (TaskOs::wallNow() - m_t0);
        if (ahead <= 0.0005)
        {
            break;
        }
        if (ahead > 0.0025)
        {
            TaskOs::sleepSec(ahead - 0.002); // грубый сон с запасом
        }
        else
        {
            TaskOs::yieldSpin(); // доспин последних долей мс
        }
    }
}

double RtPacer::lag(double t) const
{
    if (!m_anchored)
    {
        return 0.0;
    }
    return (TaskOs::wallNow() - m_t0) - t / m_speed;
}

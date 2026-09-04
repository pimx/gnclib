//============================================================================
// lib_task.h -- распределённое исполнение моделей: каждая модель работает
// в отдельном потоке и обменивается сообщениями по UDP (localhost).
//
// ИДЕНТИФИКАЦИЯ: у каждого потока свой порт приёма (TaskPort) -- порт
// фактически является идентификатором модели (задачи/потока).
//
// ДЕТЕРМИНИЗМ: конвейер тактов воспроизводит последовательный цикл
// стенда точно. Задача такта k блокируется до прихода ВСЕХ входных
// сообщений с требуемыми номерами тактов (вход имеет смещение такта:
// 0 -- сообщение текущего такта, -1 -- предыдущего), затем выполняет
// update и рассылает выход адресатам с номером такта k. На каждой
// связи в полёте не более одного сообщения -- переупорядочивание
// исключено, результаты бит-в-бит совпадают с последовательным стендом.
//
// Конвейер (смещения входов и фазы времени):
//   ModOperator   (9905): вход MsgController тактом k-1; t = k*dt
//   ModSensors    (9902): вход MsgDynamics тактом k-1;   t = k*dt
//   ModController (9903): входы MsgSensors и MsgOperator тактом k; t = k*dt
//   ModActuators  (9904): вход MsgController тактом k;   t = (k+1)*dt
//   ModDynamics   (9901): вход MsgActuators тактом k;    t = (k+1)*dt
// Затравка цикла (tick = -1): начальное MsgDynamics -> 9902 и
// MsgController по умолчанию -> 9905 (отправляет основной поток).
//
// ПАМЯТЬ: и принимаемые, и отправляемое сообщения -- в ЛОКАЛЬНОЙ памяти
// задачи (поля объекта, доступ только из её потока) -- параллельного
// доступа на запись/чтение нет. Наблюдение со стороны -- ЧЕРЕЗ МОНИТОР:
// стенд включает порт монитора в список адресатов каждой задачи, и
// монитор собирает копии всех выходных сообщений из датаграмм
// (код автора -- поле srcPort конверта). Последний выход задачи
// доступен через lastOut() ПОСЛЕ завершения её потока (для проверок).
//
// УСТОЙЧИВОСТЬ К ПОРЯДКУ ЗАПУСКА ПРОЦЕССОВ: датаграмма, отправленная на
// ещё не привязанный порт, теряется. Задача, простаивающая в ожидании
// входов, периодически (каждые ~200 мс) ПОВТОРЯЕТ отправку последнего
// выхода адресатам; получатели игнорируют такты, отличные от ожидаемых,
// поэтому дубликаты безопасны, а поздно стартовавший процесс получает
// недостающие сообщения.
//
// ПЕРЕНОСИМОСТЬ: интерфейс классов инвариантен к операционной системе;
// обмен по UDP -- библиотечные классы UdpReader/UdpSender (lib_udp.h,
// реализации lib_udp_windows.cpp / lib_udp_posix.cpp); ОС-примитивы
// потоков (запуск, ожидание, сон, барьер памяти) -- класс TaskOs,
// реализации в файлах с суффиксами ОС:
//   lib_task_windows.cpp (CreateThread + Sleep + MemoryBarrier),
//   lib_task_posix.cpp   (pthread + nanosleep + __sync_synchronize).
// Переносимая часть -- lib_task.cpp. Библиотека не зависит от типов
// сообщений проекта: полезная нагрузка -- произвольные POD-структуры
// (параметры шаблонов); порты задач назначает приложение.
//============================================================================
#ifndef LIB_TASK_H
#define LIB_TASK_H

#include <stdint.h>
#include <string.h>
#include "lib_udp.h"

//----------------------------------------------------------------------------
// TaskPacer -- интерфейс темпогенератора для замыкающей задачи кольца
// (реализуется приложением, например RtPacer)
//----------------------------------------------------------------------------
class TaskPacer
{
public:
    virtual void pace(double t) = 0;

protected:
    ~TaskPacer()
    {
    }
};

//----------------------------------------------------------------------------
// TaskOs -- ОС-примитивы потоков (реализации в lib_task_windows.cpp /
// lib_task_posix.cpp; интерфейс инвариантен к ОС)
//----------------------------------------------------------------------------
class TaskOs
{
public:
    // Потоки
    static bool threadStart(void (*fn)(void*), void* arg, long long* handle);
    static void threadJoin(long long handle);

    // Сон и полный барьер памяти
    static void sleepMs(int ms);
    static void memoryBarrier();

    // --- примитивы реального масштаба времени ---
    // Монотонные настенные часы, с (POSIX: CLOCK_MONOTONIC,
    // Windows: QueryPerformanceCounter)
    static double wallNow();
    // Сон с субмиллисекундным разрешением, с
    static void sleepSec(double s);
    // Уступить квант в цикле точного доспина
    static void yieldSpin();
    // Разрешение системного таймера: на Windows timeBeginPeriod(1) /
    // timeEndPeriod(1), на POSIX -- пустые (nanosleep точен сам по себе).
    // Вызовы парные.
    static void timerResolutionBegin();
    static void timerResolutionEnd();
};

//----------------------------------------------------------------------------
// RtPacer -- темпогенератор реального масштаба времени.
//
// Привязывает модельное время к настенным часам: после каждого такта
// pace(t) досыпает так, чтобы модельный момент t наступал не раньше
// t/speed секунд настенного времени от якоря.
//
// ОТСТАВАНИЕ. Мелкую просрочку (джиттер ОС) темпогенератор отрабатывает
// естественным образом: не спит, пока модель не догонит настенную цель.
// Просрочка СВЫШЕ порога setMaxLag() не отрабатывается вовсе -- якорь
// сдвигается вперёд, время списывается в slip(), и следующий такт идёт
// в нормальном темпе. Без этого один длинный такт (например,
// перепланирование MPC на 0.6 с) заставлял бы цикл гнать сотни тактов
// свободным темпом, и потребитель телеметрии получал бы паузу, а затем
// пачку кадров с интервалом в десятки микросекунд.
//
// Коэффициент speed -- модельных секунд на настенную секунду: 1.0
// реальное время, 2.0 вдвое быстрее, 0.5 вдвое медленнее.
//
// Якорь ЛЕНИВЫЙ, ставится первым вызовом pace(): между start() и первым
// тактом может пройти произвольное время (ожидание партнёра в
// раздельных режимах), иначе конвейер после старта мчался бы свободным
// темпом, догоняя настенную цель.
//
// На Windows одиночный Sleep на такте в единицы миллисекунд обрушивает
// темп (системный таймер 15.6 мс), поэтому pace() -- гибрид: грубый сон
// с запасом 2 мс плюс точный доспин.
//----------------------------------------------------------------------------
class RtPacer : public TaskPacer
{
public:
    RtPacer();
    ~RtPacer();

    // Начало отсчёта; speed -- модельных секунд на настенную секунду
    void start(double speed);

    double speed() const
    {
        return m_speed;
    }

    // Дождаться наступления модельного момента t (интерфейс TaskPacer)
    void pace(double t);

    // Отставание модели от реального темпа, с (> 0 -- не успеваем)
    double lag(double t) const;

    // Порог прощения отставания, с: просрочка такта СВЫШЕ порога не
    // отрабатывается ускоренным темпом, а списывается в slip().
    // Умолчание 0.020 с (8 тактов при 400 Гц) -- выше обычного джиттера
    // ОС и ниже стоимости одиночного перепланирования.
    void setMaxLag(double s)
    {
        m_maxLag = (s > 0.0) ? s : 0.0;
    }

    double maxLag() const
    {
        return m_maxLag;
    }

    // Списанное (не отработанное) время, с -- суммарное и наибольшее
    // единичное; число событий списания
    double slip() const
    {
        return m_slip;
    }

    double slipMax() const
    {
        return m_slipMax;
    }

    int slipEvents() const
    {
        return m_slipN;
    }

private:
    double m_speed; // модельных секунд на настенную секунду
    double m_t0; // настенный якорь отсчёта, с
    double m_maxLag; // порог прощения отставания, с
    double m_slip; // суммарное списанное время, с
    double m_slipMax; // наибольшее единичное списание, с
    int m_slipN; // число событий списания
    bool m_started; // start() вызывался (парность timerResolutionEnd)
    bool m_anchored; // якорь установлен первым pace()
};

//----------------------------------------------------------------------------
// Конверт сообщения (перед полезной нагрузкой в датаграмме)
//----------------------------------------------------------------------------
struct MsgEnvelope
{
    uint16_t srcPort; // код автора (типа) сообщения = порт задачи
    uint16_t size; // размер полезной нагрузки, байт
    int32_t tick; // номер такта сообщения (-1 -- затравка)
    double time;    // модельное время сообщения, с
    uint32_t magic; // MSG_MAGIC
};

enum
{
    MSG_MAGIC = 0x33594C46, // историческая сигнатура конверта
                            // (ASCII "FLY3"); менять нельзя --
        // совместимость всех потребителей протокола
    MSG_MAX_PAYLOAD = 1024, // предел полезной нагрузки, байт
    TASK_MAX_DEST = 4, // предел числа адресатов задачи
    TASK_HOST_LEN = 40, // предел длины адреса хоста
    TASK_RESYNC_GAP = 400 // порог подстройки времени, тактов (1 c)
};

//----------------------------------------------------------------------------
// Адресат выхода задачи: порт + хост (для запуска частей конвейера на
// разных вычислителях; умолчание -- localhost)
//----------------------------------------------------------------------------
struct TaskDest
{
    int port; // порт-идентификатор адресата
    char host[TASK_HOST_LEN]; // IP-адрес адресата

    TaskDest() : port(0)
    {
        host[0] = '1';
        host[1] = '2';
        host[2] = '7';
        host[3] = '.';
        host[4] = '0';
        host[5] = '.';
        host[6] = '0';
        host[7] = '.';
        host[8] = '1';
        host[9] = 0;
    }

    void set(int p, const char* h)
    {
        port = p;
        int i = 0;
        if (h == 0)
        {
            h = "127.0.0.1";
        }
        while (h[i] != 0 && i < TASK_HOST_LEN - 1)
        {
            host[i] = h[i];
            i = i + 1;
        }
        host[i] = 0;
    }
};

//----------------------------------------------------------------------------
// ModTaskBase -- общая (нешаблонная) часть задачи: сокеты, поток, цикл,
// конверты, адресаты, темп, останов
//----------------------------------------------------------------------------
class ModTaskBase
{
public:
    ModTaskBase();
    virtual ~ModTaskBase();

    // Конфигурация: порт приёма (идентификатор), адресаты выхода
    // (порт + хост), старт времени, шаг такта, фаза времени (0 или dt)
    void configure(int port, const TaskDest* dest, int nDest, double t0,
                   double dt, double phase);

    // Принимать датаграммы со всех интерфейсов (для конвейера на разных
    // вычислителях); по умолчанию -- только localhost
    void setBindAny(bool bindAny);

    // Темпогенератор (задаётся замыкающей задаче кольца -- динамике)
    void setPacer(TaskPacer* pacer);

    // ИСТОЧНИК ВРЕМЕНИ: задача-мастер (динамика) ведёт время сама
    // (накопление шага) и отправляет его в сообщениях; остальные задачи
    // берут время из входного сообщения (вход 1 -- носитель времени),
    // их результаты соответствуют этому времени
    void setTimeMaster(bool master);

    // Остановиться после публикации такта stopTick (< 0 -- без предела)
    void setStopTick(int stopTick);

    // Запуск потока задачи (сокет должен открыться до затравки цикла)
    bool start();
    // Запрос останова (только флаг, без ожидания): для группового
    // останова -- сначала requestStop всем задачам, затем stop
    void requestStop()
    {
        m_run = false;
    }
    // Останов: флаг + ожидание завершения потока
    void stop();

    int port() const
    {
        return m_port;
    }
    int tick() const
    {
        return m_tick;
    }
    // Время непрерывного простоя ожидания входов, с (сбрасывается любым
    // принятым конвертом); растёт только после первого выполненного
    // такта -- ожидание старта партнёра простоем не считается
    double idleSeconds() const
    {
        return (double)m_idleMs / 1000.0;
    }
    // Цикл задачи работает (false -- завершился по stopTick или stop)
    bool active() const
    {
        return m_run;
    }
    bool failed() const
    {
        return m_failed;
    }

protected:
    // Один такт задачи (реализуется шаблоном): ожидание входов,
    // update модели, публикация и рассылка
    virtual void cycle() = 0;

    // Приём одной датаграммы с проверкой конверта; false -- таймаут
    // (внутри -- учёт простоя и периодическая ретрансляция resendLast)
    bool recvMsg(MsgEnvelope* env, void* payload, int maxSize);
    // Повтор отправки последнего выхода (реализует шаблон; вызывается
    // из потока задачи при простое ожидания)
    virtual void resendLast()
    {
    }

    // ПОДСТРОЙКА ВРЕМЕНИ ПОД СООБЩЕНИЯ: если такт входа отличается от
    // ожидаемого больше чем на TASK_RESYNC_GAP (партнёр стартовал в
    // другой момент или был перезапущен -- в т.ч. на другом
    // вычислителе), задача принимает время партнёра: её такт и
    // накопленное время выставляются по конверту. Возврат true --
    // подстройка выполнена (ожидаемые такты пересчитать).
    bool maybeResync(const MsgEnvelope& env, int tickOffset);
    // Рассылка полезной нагрузки всем адресатам с конвертом
    void sendToAll(const void* payload, int size, int tk, double t);

    // Модельное время такта с учётом фазы. НАКОПЛЕНИЕ (m_tBase += dt),
    // а не произведение tick*dt: битовая тождественность
    // последовательному стенду, который накапливает t += dt.
    double tickTime() const
    {
        return m_tBase + m_phase;
    }
    bool running() const
    {
        return m_run;
    }
    bool timeMaster() const
    {
        return m_timeMaster;
    }

    // --- такт ---
    int m_tick; // номер текущего такта задачи

private:
    static void threadEntry(void* self);
    void runLoop();

    // --- идентификация и адресация ---
    int m_port; // порт приёма (идентификатор задачи)
    TaskDest m_dest[TASK_MAX_DEST]; // адресаты выхода (порт + хост)
    int m_nDest;                    // число адресатов
    bool m_bindAny;    // приём со всех интерфейсов
    bool m_timeMaster; // задача -- источник времени

    // --- время ---
    double m_t0;    // модельное время такта 0, с
    double m_dt;    // шаг такта, с
    double m_phase; // фаза времени задачи (0 или dt), с
    double m_tBase; // накопленное время начала текущего такта, с

    // --- исполнение ---
    gnc::UdpReader m_reader; // приёмный сокет задачи (порт-идентификатор)
    gnc::UdpSender m_sender; // сокет отправки адресатам
    long long m_thread; // описатель потока
    volatile bool m_run; // флаг работы цикла
    volatile bool m_failed; // фатальная ошибка задачи (сокет)
    volatile long m_idleMs; // простой ожидания входов, мс
    int m_timeouts; // счётчик таймаутов подряд (ретрансляция)
    int m_stopTick; // такт останова (< 0 -- нет)
    TaskPacer* m_pacer; // темпогенератор (0 -- свободный темп)
};

//----------------------------------------------------------------------------
// ModTask1 -- задача модели с одним входом:
//   Model::update(t, OutMsg&, const In1&)
//----------------------------------------------------------------------------
template <class Model, class OutMsg, class In1>
class ModTask1 : public ModTaskBase
{
public:
    ModTask1()
        : m_model(0),
          m_lastTick(-1),
          m_lastTime(0.0),
          m_hasOut(false),
          m_in1Port(0),
          m_in1Off(0)
    {
    }

    // Привязка: модель, порт-источник входа и смещение его такта (0/-1)
    void bind(Model* model, int in1Port, int in1TickOffset)
    {
        m_model = model;
        m_in1Port = in1Port;
        m_in1Off = in1TickOffset;
    }

    // Последний выход задачи (читать ПОСЛЕ завершения потока)
    const OutMsg& lastOut() const
    {
        return m_out;
    }

protected:
    void cycle()
    {
        // ожидание входа текущего такта
        int need1 = m_tick + m_in1Off;
        bool got1 = false;
        double tIn1 = 0.0;
        MsgEnvelope env;
        unsigned char buf[MSG_MAX_PAYLOAD];
        while (running() && !got1)
        {
            if (!recvMsg(&env, buf, sizeof(buf)))
            {
                continue; // таймаут (учёт простоя/ретрансляция в базе)
            }
            if (env.srcPort == m_in1Port && maybeResync(env, m_in1Off))
            {
                need1 = m_tick + m_in1Off; // подстройка времени выполнена
            }
            if (env.srcPort == m_in1Port && env.tick == need1 &&
                env.size == (int)sizeof(In1))
            {
                memcpy(&m_in1, buf, sizeof(In1));
                tIn1 = env.time; // носитель времени -- вход 1
                got1 = true;
            }
        }
        if (!running())
        {
            return;
        }

        // update модели в локальный выход; рассылка адресатам
        // (включая порт монитора, если он в списке адресатов).
        // Время: мастер (динамика) ведёт его сам; остальные берут из
        // входного сообщения -- результат соответствует этому времени.
        double t = timeMaster() ? tickTime() : tIn1;
        m_model->update(t, m_out, m_in1);
        m_lastTick = m_tick;
        m_lastTime = t;
        m_hasOut = true;
        sendToAll(&m_out, (int)sizeof(OutMsg), m_tick, t);
    }

    void resendLast()
    {
        if (m_hasOut)
        {
            sendToAll(&m_out, (int)sizeof(OutMsg), m_lastTick, m_lastTime);
        }
    }

private:
    // --- привязка ---
    Model* m_model; // модель (владелец -- стенд)

    // --- выход (локальная память задачи) ---
    OutMsg m_out;      // последнее выходное сообщение
    int m_lastTick;    // такт последнего выхода
    double m_lastTime; // время последнего выхода, с
    bool m_hasOut; // выход публиковался (есть что ретранслировать)

    // --- вход (локальная память задачи) ---
    In1 m_in1; // последнее принятое входное сообщение
    int m_in1Port; // порт-идентификатор источника
    int m_in1Off;  // смещение такта входа (0 или -1)
};

//----------------------------------------------------------------------------
// ModTask2 -- задача модели с двумя входами:
//   Model::update(t, OutMsg&, const In1&, const In2&)
//----------------------------------------------------------------------------
template <class Model, class OutMsg, class In1, class In2>
class ModTask2 : public ModTaskBase
{
public:
    ModTask2()
        : m_model(0),
          m_lastTick(-1),
          m_lastTime(0.0),
          m_hasOut(false),
          m_in1Port(0),
          m_in1Off(0),
          m_in2Port(0),
          m_in2Off(0)
    {
    }

    void bind(Model* model, int in1Port, int in1TickOffset, int in2Port,
              int in2TickOffset)
    {
        m_model = model;
        m_in1Port = in1Port;
        m_in1Off = in1TickOffset;
        m_in2Port = in2Port;
        m_in2Off = in2TickOffset;
    }

    // Последний выход задачи (читать ПОСЛЕ завершения потока)
    const OutMsg& lastOut() const
    {
        return m_out;
    }

protected:
    void cycle()
    {
        int need1 = m_tick + m_in1Off;
        int need2 = m_tick + m_in2Off;
        bool got1 = false;
        bool got2 = false;
        double tIn1 = 0.0;
        MsgEnvelope env;
        unsigned char buf[MSG_MAX_PAYLOAD];
        while (running() && !(got1 && got2))
        {
            if (!recvMsg(&env, buf, sizeof(buf)))
            {
                continue; // таймаут (учёт простоя/ретрансляция в базе)
            }
            if ((env.srcPort == m_in1Port && maybeResync(env, m_in1Off)) ||
                (env.srcPort == m_in2Port && maybeResync(env, m_in2Off)))
            {
                // подстройка времени: пересчёт ожиданий, сбор заново
                need1 = m_tick + m_in1Off;
                need2 = m_tick + m_in2Off;
                got1 = false;
                got2 = false;
            }
            if (env.srcPort == m_in1Port && env.tick == need1 &&
                env.size == (int)sizeof(In1))
            {
                memcpy(&m_in1, buf, sizeof(In1));
                tIn1 = env.time; // носитель времени -- вход 1
                got1 = true;
            }
            else if (env.srcPort == m_in2Port && env.tick == need2 &&
                     env.size == (int)sizeof(In2))
            {
                memcpy(&m_in2, buf, sizeof(In2));
                got2 = true;
            }
        }
        if (!running())
        {
            return;
        }

        double t = timeMaster() ? tickTime() : tIn1;
        m_model->update(t, m_out, m_in1, m_in2);
        m_lastTick = m_tick;
        m_lastTime = t;
        m_hasOut = true;
        sendToAll(&m_out, (int)sizeof(OutMsg), m_tick, t);
    }

    void resendLast()
    {
        if (m_hasOut)
        {
            sendToAll(&m_out, (int)sizeof(OutMsg), m_lastTick, m_lastTime);
        }
    }

private:
    // --- привязка ---
    Model* m_model; // модель (владелец -- стенд)

    // --- выход (локальная память задачи) ---
    OutMsg m_out;      // последнее выходное сообщение
    int m_lastTick;    // такт последнего выхода
    double m_lastTime; // время последнего выхода, с
    bool m_hasOut; // выход публиковался (есть что ретранслировать)

    // --- входы (локальная память задачи) ---
    In1 m_in1;     // вход 1
    int m_in1Port; // порт-идентификатор источника 1
    int m_in1Off;  // смещение такта входа 1
    In2 m_in2;     // вход 2
    int m_in2Port; // порт-идентификатор источника 2
    int m_in2Off;  // смещение такта входа 2
};

//----------------------------------------------------------------------------
// Отправка затравки цикла (tick = -1) от основного потока
//----------------------------------------------------------------------------
bool taskSendPrime(gnc::UdpSender& sender, int destPort, int srcPort,
                   const void* payload, int size, double t);

#endif // LIB_TASK_H

//============================================================================
// lib_telemetry.h -- библиотечный класс отправки UDP-телеметрии.
// ПОЛНАЯ СПЕЦИФИКАЦИЯ ПАКЕТА -- в TELEMETRY.md проекта (контракт
// приёмника-визуализатора).
//
// Пакет v2 -- double[17] в родном порядке байт (136 байт):
//   [0] время, с; [1..3] pos ENU, м; [4..7] кватернион BODY->ENU
//   (Гамильтон, скаляр первым); [8..10] вектор тяги ENU, Н;
//   [11] масса, кг; [12..14] vel ENU, м/с; [15] режим (0..7);
//   [16] rcsMask -- маска сопел РСУ, ИЛИ-защёлка между отправками.
//
// Класс библиотечный: не зависит от типов сообщений проекта -- поля
// пакета передаются явно (Vec3/Quat из lib_linalg). Прореживание и
// ИЛИ-защёлка маски РСУ -- внутри. Отправка -- UdpSender (lib_udp).
//============================================================================
#ifndef LIB_TELEMETRY_H
#define LIB_TELEMETRY_H

#include "lib_frames.h"
#include "lib_udp.h"

class UdpTelemetry
{
public:
    enum
    {
        PKT_DOUBLES = 17, // размер пакета v2, чисел double
        PKT3_DOUBLES = 18, // размер пакета v3 (этап 217): + dFeet
        PKT4_DOUBLES = 20  // v4 (этап 238): + n2Left, propLeft
    };

    struct Config
    {
        bool enable;      // телеметрия включена
        const char* host; // адрес приёмника (0 -> 127.0.0.1)
        int port;         // порт приёмника (<= 0 -> 9000)
        int decimation;   // слать каждый N-й вызов (>= 1)

        Config() : enable(false), host(0), port(0), decimation(1)
        {
        }
    };

    UdpTelemetry();
    ~UdpTelemetry();

    bool setup(const Config& c);
    void close();
    bool isOpen() const
    {
        return m_sender.isOpen();
    }

    // Сбор и отправка пакета v2. Маска РСУ копится ИЛИ-защёлкой между
    // фактическими отправками (прореживание не теряет короткие импульсы).
    //   t            -- модельное время, с
    //   posEnu       -- положение ЦМ в ENU, м
    //   qBodyToEnu   -- ориентация BODY->ENU (Гамильтон, скаляр первым)
    //   thrustEnu    -- вектор тяги в ENU, Н
    //   mass         -- масса, кг
    //   velEnu       -- скорость ЦМ в ENU, м/с
    //   mode         -- режим контроллера (0..7)
    //   rcsMask      -- маска работающих сопел РСУ (биты 0..3)
    void send(double t, VecEnu posEnu, Quat qBodyToEnu, VecEnu thrustEnu,
              double mass, VecEnu velEnu, int mode, int rcsMask);

    // Пакет ФОРМАТА 3 (этап 217): те же 17 даблов + [17] dFeet --
    // расстояние от ЦМ до опор вдоль оси аппарата, м (x_cm + клиренс).
    // Опоры на поверхности при posEnu.u == dFeet на вертикальной стойке;
    // касание отслеживается по опорам. Формат v2 не менялся.
    void sendV3(double t, VecEnu posEnu, Quat qBodyToEnu, VecEnu thrustEnu,
                double mass, VecEnu velEnu, int mode, int rcsMask,
                double dFeet);

    /* пакет v4 (этап 238): v3 + остаток азота и остаток топлива */
    void sendV4(double t, VecEnu posEnu, Quat qBodyToEnu, VecEnu thrustEnu,
                double mass, VecEnu velEnu, int mode, int rcsMask,
                double dFeet, double n2Left, double propLeft);

private:
    // --- конфигурация ---
    Config m_cfg; // адрес, порт, прореживание

    // --- отправка ---
    UdpSender m_sender; // сокет отправки (lib_udp)
    long m_counter;     // счётчик вызовов (прореживание)
    int m_maskAccum;    // ИЛИ-защёлка маски РСУ между пакетами
};

#endif // LIB_TELEMETRY_H

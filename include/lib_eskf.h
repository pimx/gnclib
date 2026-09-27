// gnclib 1.1.0: navigation ESKF (moved from FLY90 / HOP90 V90.43, with
// the Earth rate, IMU-lever noise and gate recovery of HOP90).
//============================================================================
// eskf.h -- HOP90 (перенос FLY90 V90.4 models/eskf.h; добавлен
// радиовысотомер) -- комплексирование ИНС / СНС / высотомеров: фильтр Калмана в
// пространстве ошибок (ESKF, 15 состояний):
//   [ dp(3) dv(3) dTheta(3) dba(3) dbg(3) ]
// dp, dv -- ошибки положения и скорости на площадке; dTheta -- малый
// поворот ошибки ориентации в осях BODY (правое умножение:
// R_ист = R_оц * exp([dTheta]x)); dba, dbg -- ошибки оценок смещений.
//
// Механизация БИНС на каждом такте ИНС (400 Гц); ковариация --
// дискретизация 1-го порядка Fd = I + F*dt. Коррекции ПОСЛЕДОВАТЕЛЬНО
// ПО СКАЛЯРАМ (обращения матриц нет; после каждой ошибка вводится в
// номинал). СНС -- позиция и скорость ENU, 6 скаляров; высотомер --
// коррекция ОГРАНИЧЕНА вертикальным каналом (dp_u, dv_u, dba_bu, форма
// Джозефа): один луч не разделяет высоту и наклон. Адаптивный
// процессный шум при манёврах (пик-удержание с распадом ~1.5 с).
//
// Реализация: eskf_predict.cpp (настройка, механизация, ковариация),
// eskf_update.cpp (скалярные коррекции). Матричная форма -- eskf_math.h.
//============================================================================
#ifndef GNCLIB_LIB_ESKF_H
#define GNCLIB_LIB_ESKF_H

#include "lib_nav_types.h"
#include "lib_eskf_math.h"
#include "lib.h"

#include <math.h>

namespace gnc
{

class Eskf
{
public:
    static const int N = 15;

    // Параметры шумов фильтра.
    struct Config
    {
        // gnclib 1.1.0 (HOP90): lever of the IMU -- smoothing of the gyro
        // derivative (s; 0 -- the raw derivative, as FLY90) and the adaptive
        // lever noise in the velocity process noise (0 -- off, as FLY90)
        double wDotTau = 0.0;
        int leverNoise = 0;
        double leverNoiseK = 1.0; // scale of the lever noise sigma (FLY90 tuning)
        double gyroNoiseStd;  // СКО гироскопов на такте, рад/с
        double accelNoiseStd; // СКО акселерометров на такте, м/с^2
        double gyroBiasRw;    // блуждание смещения гироскопов, рад/с/sqrt(такт)
        double accelBiasRw;   // блуждание смещения акселерометров
        double gnssPosStd;    // СКО позиции СНС, м
        double gnssVelStd;    // СКО скорости СНС, м/с
        double altNoiseStd;   // СКО лазерного высотомера, м
        double raltNoiseStd;  // СКО радиовысотомера: постоянная часть, м
        double raltNoiseRel;  // СКО радиовысотомера: доля дальности

        // Конструктор: значения по умолчанию (переопределяет модель
        // датчиков от паспортов).
        Config()
            : gyroNoiseStd(3e-4), accelNoiseStd(8e-3), gyroBiasRw(5e-7),
              accelBiasRw(5e-6), gnssPosStd(0.1), gnssVelStd(0.05),
              altNoiseStd(0.02), raltNoiseStd(0.6), raltNoiseRel(0.01)
        {
        }
    };

    // Конструктор: фильтр не настроен, тяжесть стандартная.
    Eskf();

    // Рычаг ИНС относительно ЦМ: акселерометры измеряют добавку
    // omgDot x r + omg x (omg x r); без её вычитания угловые манёвры
    // подмешиваются в линейный канал.
    //   r -- точка установки ИНС (связанная), м
    // Установка рычага высотомера относительно ЦМ (V90.14: ЦМ смещается при
    // выработке топлива).
    //   r -- высотомер относительно ЦМ (связанная), м
    void setAltLever(const gnc::CfBodyVector& r)
    {
        m_lever = r;
    }

    void setImuLever(const gnc::CfBodyVector& r)
    {
        m_imuLever = r;
    }

    // Настройка фильтра.
    //   c -- шумы; frame -- площадка (ENU <-> ECEF); gravity -- тяжесть
    //   площадки, м/с^2; leverBody -- высотомер относительно ЦМ, м
    void setup(const Config& c, const gnc::GEO* frame, double gravity,
               gnc::CfBodyVector leverBody);

    // Инициализация точно известным начальным состоянием.
    //   p0, v0 -- положение и скорость (площадка); att0 -- ориентация
    //   t0 -- время, с
    void init(gnc::CfLocalVector p0, gnc::CfLocalVector v0,
              const gnc::CfRotation& att0, Time t0);

    // Один такт комплексирования (на каждом такте ИНС).
    //   outNav -- навигационное решение (выход)
    //   imu, gnss, alt -- измерения такта; t -- время, с
    void update(NavData& outNav, const ImuData& imu, const GnssData& gnss,
                const AltimeterData& alt, const AltimeterData& ralt, Time t);

    // Оценка смещения нуля акселерометров.
    //   возврат: смещение (связанная), м/с^2
    gnc::CfBodyVector accelBias() const
    {
        return m_ba;
    }

    // Оценка смещения нуля гироскопов.
    //   возврат: смещение (связанная), рад/с
    gnc::CfBodyVector gyroBias() const
    {
        return m_bg;
    }

private:
    // Заполнение решения для контроллера.
    //   out -- решение (выход)
    void fillOut(NavData& out) const;

    // Приведение акселерометров к ЦМ по рычагу ИНС (eskf_predict.cpp).
    //   w -- скомпенсированная угловая скорость; a -- удельная сила
    //   (вход-выход); dt -- такт, с
    void leverCompensate(const gnc::CfBodyVector& w, gnc::CfBodyVector* a,
                         double dt);

    // Распространение ковариации P = Fd P Fd^T + Qd (eskf_predict.cpp).
    //   Reb -- ориентация; a, w -- удельная сила и угловая скорость
    //   dt -- такт, с; qScale -- манёвренный масштаб Q
    void propagate(const Mat3& Reb, const gnc::CfBodyVector& a,
                   const gnc::CfBodyVector& w, double dt, double qScale);

    // Коррекция СНС: 6 скаляров позиции и скорости (eskf_update.cpp).
    //   gnss -- измерение
    void gnssCorrect(const GnssData& gnss);

    // Коррекция высотомера: вертикальный канал (eskf_update.cpp).
    //   alt -- измерение
    void altCorrect(const AltimeterData& alt, double std0, double stdRel);

    // Скалярная коррекция z = c*x[idx] + шум(R), строб 6 сигм, ввод
    // ошибки в номинал (eskf_update.cpp).
    //   idx -- индекс состояния; c -- коэффициент; innov -- невязка
    //   R -- дисперсия измерения
    void scalarUpdate(int idx, double c, double innov, double R);

    // Скалярная коррекция с усилением, ограниченным вертикальным каналом
    // {dp_u, dv_u, dba_bu}; ковариация -- форма Джозефа (eskf_update.cpp).
    //   idx, c, innov, R -- как scalarUpdate
    void scalarUpdateMasked(int idx, double c, double innov, double R);

    // Симметризация ковариации (численная гигиена).
    void symmetrize();

    // Установка диагонального блока 3x3 ковариации.
    //   i0 -- первый индекс блока; v -- дисперсия
    void setDiag(int i0, double v);

    Config m_cfg;            // шумы
    const gnc::GEO* m_frame; // площадка
    double m_g;              // тяжесть площадки, м/с^2

    double m_qScale;  // манёвренный масштаб Q (пик-удержание)
    double m_gnssAge; // время с последней коррекции СНС, с
    Time m_time;      // время последнего такта, с
    gnc::CfBodyVector m_wLast;      // угловая скорость последнего такта
    gnc::CfBodyVector m_lever;      // высотомер относительно ЦМ
    gnc::CfBodyVector m_imuLever;   // ИНС относительно ЦМ
    gnc::CfBodyVector m_wPrevLever; // угловая скорость прошлого такта
    bool m_leverInit;               // разность гироскопа определена
    Mat3 m_R;                       // ориентация (матрица, eskf_math.h)
    gnc::CfLocalVector m_v;         // скорость (площадка)
    gnc::CfLocalVector m_p;         // положение (площадка)
    gnc::CfBodyVector m_ba;         // смещение акселерометров (связанная)
    gnc::CfBodyVector m_bg;         // смещение гироскопов (связанная)
    double m_P[N][N];               // ковариация ошибок
    int m_rejectRun;                // V90.35: подряд отброшенных стробом измерений
    gnc::CfBodyVector m_wDotF;      // V90.35: сглаженное угловое ускорение (рычаг ИНС)
    double m_wDotVar;
    gnc::CfBodyVector m_wDotV;      // 20 ms smoothed derivative for the lever noise variance               // V90.36: дисперсия неучтённого углового ускорения, (рад/с^2)^2
};

}   // namespace gnc

#endif // GNCLIB_LIB_ESKF_H

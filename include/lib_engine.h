//============================================================================
// lib_engine.h -- силовая установка из нескольких сопел (1.7.0): описание
// сопел, закон тяги, направления тяги и их производные, суммарные сила,
// момент и расход.
//
// Сопла одинаковые по тяге и закону дросселирования, управляются
// индивидуально: дроссель (относительный секундный расход) и до двух
// сервоприводов отклонения. Суммарная сила -- сумма векторов тяг сопел
// без интерференции. Линия тяги проходит через шарнир сопла, поэтому
// момент определяется положением шарнира и направлением тяги.
//
// Связанная система -- CfBodyVector (bu -- ось аппарата к носу, be, bn).
// "Тяга" -- проекция суммарной силы на bu, крен -- момент вокруг bu.
//
// Закон тяги сопла: F = kf * g0 * Ivac * mdot - s * pa * Sa, F >= 0;
// mdot = u * mdotMax, u -- относительная команда дросселя в
// [throttleMin, throttleMax]; kf -- множитель разброса (только модель
// объекта, распределитель его не знает); s = 1 -- потеря на противодавление
// не зависит от дросселя (ENGINE_LOSS_CONSTANT), s = u -- пропорциональна
// дросселю (ENGINE_LOSS_THROTTLE: F = u (Fvac - (Fvac - Fsl) pa / p0) при
// Sa = (Fvac - Fsl) / p0 -- линейная по давлению интерполяция между
// паспортными тягами одного двигателя в пустоте и у земли).
//
// Направление тяги (сила на аппарат):
//   неподвижное сопло:   d = a;
//   в одной плоскости:   d = cos(d1) a + sin(d1) e1;
//   произвольно (кардан, d1 -- внешняя рамка, d2 -- внутренняя):
//                        d = cos(d2) (cos(d1) a + sin(d1) e1) + sin(d2) e2;
// a -- направление при нулевых углах, e1, e2 -- направления положительных
// отклонений (единичные, взаимно перпендикулярные и перпендикулярные a).
// Угол между d и a: cos = cos(d1) cos(d2) -- по нему предел конусом.
//
// C++17, без STL, без кучи, без исключений.
//============================================================================
#ifndef GNCLIB_LIB_ENGINE_H
#define GNCLIB_LIB_ENGINE_H

#include "lib_frames.h"

namespace gnc
{

enum
{
    ENGINE_MAX_NOZZLES = 16, // предел числа сопел
    ENGINE_DELAY_MAX = 128   // предел запаздывания, шагов модели
};

// Вид сопла по отклонению
enum NozzleKind
{
    NOZZLE_FIXED = 0, // неподвижное
    NOZZLE_PLANE = 1, // в одной плоскости (один сервопривод)
    NOZZLE_FREE = 2   // в произвольном направлении (два сервопривода)
};

// Закон потери тяги на противодавление
enum EngineLossLaw
{
    ENGINE_LOSS_CONSTANT = 0, // pa * Sa, не зависит от дросселя
    ENGINE_LOSS_THROTTLE = 1  // u * pa * Sa, пропорциональна дросселю
};

// Стандартное ускорение свободного падения (определение Iуд), м/с^2
inline constexpr double ENGINE_G0 = 9.80665;

// Описание сопла
struct NozzleConfig
{
    int kind;              // NozzleKind
    CfBodyVector hinge;    // шарнир, м
    CfBodyVector axis;     // направление тяги при нулевых углах
    CfBodyVector e1;       // положительное отклонение 1
    CfBodyVector e2;       // положительное отклонение 2 (NOZZLE_FREE)
    double limit1;         // предел |d1|, рад
    double limit2;         // предел |d2|, рад
    double cone;           // предел полного отклонения (NOZZLE_FREE), рад
    double rate;           // предел скорости сервопривода, рад/с (<= 0 -- нет)
    double delay;          // запаздывание сервоприводов, с
    double tau;            // постоянная времени сервопривода, с (0 -- нет)
    double wn;             // собственная частота привода 2-го порядка, рад/с
                           // (0 -- нет; при wn > 0 tau не действует)
    double zeta;           // относительное демпфирование привода 2-го порядка
    double thrustFactor;   // модель объекта: множитель тяги (1 -- номинал)
    CfBodyVector misalign; // модель объекта: перекос оси тяги, малый поворот, рад
};

// Описание силовой установки
struct EngineConfig
{
    int count;                               // число сопел
    NozzleConfig nozzle[ENGINE_MAX_NOZZLES]; // сопла
    double ivac;                             // удельный импульс в пустоте, с
    double mdotMax;                          // расход сопла при u = 1, кг/с
    double exitArea;                         // площадь среза сопла, м^2
    double throttleMin;                      // нижний предел дросселя
    double throttleMax;                      // верхний предел дросселя
    double throttleRate;                     // предел скорости дросселя, 1/с (<= 0 -- нет)
    double throttleDelay;                    // запаздывание дросселя, с
    double throttleTau;                      // постоянная времени дросселя, с (0 -- нет)
    int lossLaw;                             // EngineLossLaw
    double startRate;                        // модель объекта: нарастание дросселя от нуля
                                             // после включения сопла, 1/с (0 -- сразу)
    double shutdownRate;                     // модель объекта: спад дросселя до нуля после
                                             // выключения сопла, 1/с (0 -- сразу)
};

// Команды (или фактические положения) органов: по соплу -- дроссель и
// углы сервоприводов (неиспользуемые углы -- нули)
struct EngineCommand
{
    double throttle[ENGINE_MAX_NOZZLES]; // относительный расход
    double angle1[ENGINE_MAX_NOZZLES];   // угол d1, рад
    double angle2[ENGINE_MAX_NOZZLES];   // угол d2, рад
};

// Суммарные силовые характеристики: вектор тяги, его составляющие и
// направление, момент относительно заданной точки, расход, тяги сопел
struct EngineForces
{
    CfBodyVector force;                // суммарный вектор тяги, Н
    double norm;                       // модуль вектора тяги, Н
    double axial;                      // осевая составляющая (по bu), Н
    double lateral;                    // боковая сила (модуль составляющей be, bn), Н
    double tilt;                       // угол вектора тяги от оси bu, рад
    double azimuth;                    // азимут боковой силы от be к bn, рад
    CfBodyVector point;                // точка приведения момента, м
    CfBodyMoment moment;               // момент относительно point, Н*м
    double mdot;                       // суммарный секундный расход, кг/с
    double thrust[ENGINE_MAX_NOZZLES]; // тяга сопел, Н
};

// Сопло по умолчанию: неподвижное, тяга вдоль bu, e1 = be, e2 = bn,
// пределы нулевые, без запаздывания, без привода 2-го порядка, без разброса.
//   n -- заполняемое описание
void NozzleInit(NozzleConfig* n);

// Установка по умолчанию: сопел нет, закон потерь ENGINE_LOSS_CONSTANT,
// пределы и темпы нулевые (заполняет потребитель).
//   c -- заполняемое описание
void EngineInit(EngineConfig* c);

// Проверка описания установки: число сопел, единичность и
// перпендикулярность направлений, пределы, запаздывания. Диагностика --
// stderr, ASCII. Возврат: true -- описание годно.
//   c -- описание
bool EngineConfigCheck(const EngineConfig& c);

// Тяга сопла по закону F = kf * g0 * Ivac * u * mdotMax - s * pa * Sa,
// F >= 0 (s -- по EngineConfig::lossLaw).
//   c -- описание; u -- дроссель; pa -- давление среды, Па;
//   factor -- множитель разброса kf (1 -- номинал)
double NozzleThrust(const EngineConfig& c, double u, double pa, double factor);

// Дроссель, при котором номинальная тяга сопла равна f (обратный закон
// тяги без ограничения диапазоном).
//   c -- описание; f -- тяга сопла, Н; pa -- давление среды, Па
double NozzleThrottleFor(const EngineConfig& c, double f, double pa);

// Производная тяги по дросселю (0, если тяга упирается в ноль).
//   c, u, pa, factor -- как у NozzleThrust
double NozzleThrustSlope(const EngineConfig& c, double u, double pa, double factor);

// Номинальное направление тяги сопла (без перекоса).
//   n -- сопло; d1, d2 -- углы сервоприводов, рад
CfBodyVector NozzleDirection(const NozzleConfig& n, double d1, double d2);

// Производные номинального направления по углам (у неподвижного --
// нули, у сопла в плоскости dd2 -- ноль).
//   n -- сопло; d1, d2 -- углы, рад; dd1, dd2 -- результаты
void NozzleDirectionDerivs(const NozzleConfig& n, double d1, double d2,
                           CfBodyVector* dd1, CfBodyVector* dd2);

// Направление с перекосом сопла: d + misalign x d, нормированное.
//   n -- сопло; d -- номинальное направление
CfBodyVector NozzleMisaligned(const NozzleConfig& n, const CfBodyVector& d);

// Признак: углы внутри конуса, cos(d1) cos(d2) >= cos(cone) - 1e-15
// (допуск округления: углы на границе считаются внутри).
//   d1, d2 -- углы, рад; cone -- предел, рад
bool ConeInside(double d1, double d2, double cone);

// Пропорциональное уменьшение углов до границы конуса (если снаружи).
// Возврат: true -- углы уменьшены.
//   d1, d2 -- углы, рад (вход и выход); cone -- предел, рад
bool ConeScale(double* d1, double* d2, double cone);

// Ограничение углов сопла пределами сервоприводов и конусом.
//   n -- сопло; d1, d2 -- углы, рад (вход и выход)
void NozzleClampAngles(const NozzleConfig& n, double* d1, double* d2);

// Суммарные сила (вектор, модуль, осевая и боковая составляющие,
// направление), момент относительно точки и расход по положениям органов.
// Выключенные сопла (on[i] == 0) тяги и расхода не дают.
//   c -- описание; x -- положения органов; on -- признаки работы сопел
//   (0 -- все работают); pa -- давление среды, Па; point -- точка
//   приведения момента; dispersion -- учитывать разброс и перекос (модель
//   объекта); out -- результат
void EngineForcesOf(const EngineConfig& c, const EngineCommand& x,
                    const unsigned char* on, double pa, const CfBodyVector& point,
                    bool dispersion, EngineForces* out);

} // namespace gnc

#endif // GNCLIB_LIB_ENGINE_H

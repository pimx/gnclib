//============================================================================
// lib_engine_alloc.h -- распределение требований к силовой установке по
// органам (1.7.0): дросселям и сервоприводам сопел (lib_engine.h).
//
// Два режима (AllocRequest::mode):
//   ALLOC_THRUST   -- заданы тяга T (проекция суммарной силы на bu) и
//                     момент M относительно точки; неизвестные -- дроссели
//                     работающих сопел и углы;
//   ALLOC_THROTTLE -- заданы один дроссель для всех работающих сопел и
//                     момент M; неизвестные -- только углы, тяга -- следствие.
// Неизвестных больше, чем уравнений; выбор -- взвешенными наименьшими
// квадратами с границами и ограничениями пар углов (lib_lsc.h):
//   min  gamma (wT^2 eT^2 + wR^2 eR^2 + wP^2 (eY^2 + eZ^2))
//      + wU^2 sum ((u - ud) / su)^2 + wA^2 sum (d / sd)^2 + wL^2 (Fbe^2 + Fbn^2) / Tn^2,
// eT = (Fbu - T) / Tn (в режиме ALLOC_THROTTLE строки нет), eR, eY, eZ --
// невязки моментов по bu, be, bn, делённые на Mn; Tn -- наибольшая тяга
// работающих сопел в пустоте, Mn = Tn (|point| + наибольший вынос шарнира);
// ud -- равномерный дроссель, дающий T без отклонений; su -- диапазон
// дросселя, sd -- предел угла. Достижимое требование выполняется с
// невязкой ~ 1 / gamma; недостижимое -- со строгим приоритетом групп в
// порядке убывания весов (по умолчанию тангаж и рыскание wP, затем крен wR,
// затем тяга wT): группы решаются по очереди -- старшие с найденными
// достижимыми значениями, текущая с весом 1e-3 от своего (не возмущает
// старшие), младшие исключены; достигнутое значение текущей группы
// становится её целью; в конце -- общее решение по достижимым целям.
//
// Границы: пределы положения (углы -- не больше программного предела
// angleMax, если он задан) и окно скорости: предыдущая команда +- rate * dt
// (команды выполнимы приводами без насыщения скорости; в режиме
// ALLOC_THROTTLE заданный дроссель тоже ведётся с этим темпом). Если
// предыдущая команда вне пределов положения (сузился программный предел)
// и окно до них не достаёт, предел важнее: команда переходит на предел за
// один шаг (у сопла с двумя приводами -- если окно не пересекает конус). Конус сопла
// с двумя приводами (cos d1 cos d2 >= cos c) -- выпуклая область углов;
// она заменяется касательными (линейные ограничения пары углов, внешние
// приближения конуса): восемь через 45 град и касательная в точке
// линеаризации, если сопло стоит на границе; решение, вышедшее из конуса,
// добавляет касательную в своей радиальной проекции на границу, задача
// решается повторно (отсекающие плоскости, до 24 раз за итерацию). Задача
// каждой итерации -- выпуклая: сопло уходит с границы или скользит по ней
// по знаку множителя; остаток выхода за конус (< 1e-9 c) снимается
// радиальным масштабированием.
// Нелинейность -- итерациями (до iterations, выход по сходимости): вдали от
// решения производные направлений без составляющей вдоль оси сопла, вблизи
// (невязка < 1e-3) -- точные; две последние итерации -- восстановление
// выполнимости (близость к текущей точке с весом 100).
// Распределитель знает номинальную модель (без разброса и перекоса) и
// запаздывания не компенсирует. Результат всегда содержит команды,
// суммарный вектор тяги (модуль, осевая и боковая составляющие,
// направление), момент относительно точки и расход.
// C++17, без STL, без кучи (объект ~60 КБ -- объявлять static).
//============================================================================
#ifndef GNCLIB_LIB_ENGINE_ALLOC_H
#define GNCLIB_LIB_ENGINE_ALLOC_H

#include "lib_engine.h"
#include "lib_lsc.h"

namespace gnc
{

enum
{
    ALLOC_MAX_VARS = 3 * ENGINE_MAX_NOZZLES, // дроссель и до двух углов на сопло
    ALLOC_MAX_CUTS = 32,                     // касательных к конусу на сопло (линеаризация)
    ALLOC_SEED_CUTS = 8                      // из них начальных (через 45 град)
};

// Режим распределения
enum AllocMode
{
    ALLOC_THRUST = 0,  // заданы тяга и момент
    ALLOC_THROTTLE = 1 // заданы дроссель (один для всех) и момент
};

// Настройки распределителя
struct AllocConfig
{
    double wThrust;   // вес невязки тяги (младший приоритет)
    double wRoll;     // вес невязки момента крена
    double wPitch;    // вес невязок моментов тангажа и рыскания (старший)
    double wThrottle; // вес отклонения дросселей от равномерного
    double wAngle;    // вес отклонения углов от нуля
    double wLateral;  // вес боковой силы (0 -- не учитывается)
    double gamma;     // множитель приоритетной части
    int iterations;   // предел итераций по нелинейной модели
    int rateWindow;   // 1 -- окно скорости по предыдущей команде
    double dt;        // шаг управления (для окна скорости), с
    double feasTol;   // допуск нормированной невязки "требование выполнено"
};

// Умолчания: wThrust 1, wRoll 30, wPitch 1000, wThrottle 1, wAngle 1,
// wLateral 0, gamma 1e8, iterations 8, rateWindow 1, dt 0.01, feasTol 1e-6.
//   a -- заполняемые настройки
void AllocConfigDefaults(AllocConfig* a);

// Требование к установке (обнулённое -- режим ALLOC_THRUST, без
// программного предела)
struct AllocRequest
{
    int mode;                             // AllocMode
    double thrust;                        // ALLOC_THRUST: тяга (проекция на bu), Н
    double throttle;                      // ALLOC_THROTTLE: дроссель работающих сопел
    CfBodyMoment moment;                  // момент относительно point, Н*м
    CfBodyVector point;                   // точка приведения момента, м
    double pa;                            // давление среды, Па
    double angleMax;                      // программный предел углов и конуса, рад (0 -- нет)
    unsigned char on[ENGINE_MAX_NOZZLES]; // 1 -- сопло работает
};

// Результат распределения (по номинальной модели)
struct AllocResult
{
    int mode;              // AllocMode требования
    EngineCommand command; // команды органов (управление)
    EngineForces forces;   // вектор тяги, его составляющие и направление, момент
                           // относительно точки требования, расход, тяги сопел
    double throttle;       // средний дроссель работающих сопел
    int throttleLimited;   // ALLOC_THROTTLE: 1 -- заданный дроссель ограничен
                           // диапазоном или темпом
    double errThrust;      // нормированная невязка тяги (ALLOC_THROTTLE -- 0)
    double errRoll;        // нормированная невязка крена
    double errPitch;       // наибольшая из нормированных невязок тангажа, рыскания
    double scaleThrust;    // Tn, Н
    double scaleMoment;    // Mn, Н*м
    int feasible;          // 1 -- все невязки в допуске
    int saturated;         // активных ограничений: границы и касательные конуса
    int iterations;        // итераций активного набора, всего
    int outer;             // итераций по нелинейной модели
};

class EngineAllocator
{
public:
    EngineAllocator();

    // Инициализация. Описания не копируются (время жизни -- у
    // вызывающего). Возврат: false -- описание негодно.
    //   cfg -- описание установки; ac -- настройки
    bool init(const EngineConfig* cfg, const AllocConfig* ac);

    // Предыдущая команда (начало окна скорости), ограничивается пределами.
    //   x -- положения органов
    void reset(const EngineCommand& x);

    // Распределение. Возврат: false -- не инициализирован или решатель не
    // сошёлся (команда -- последняя допустимая точка).
    //   req -- требование; cmd -- команды органов; res -- результат (может быть 0)
    bool allocate(const AllocRequest& req, EngineCommand* cmd, AllocResult* res);

    // Последняя выданная команда
    const EngineCommand& last() const
    {
        return m_prev;
    }

    // Число переменных задачи
    int variables() const
    {
        return m_n;
    }

private:
    void buildMap();
    void setBounds(const AllocRequest& req, int* limited);
    void desired(const AllocRequest& req, int nOn, double* xd) const;
    void linearize(const AllocRequest& req, const double* x, bool exact);
    double residual(const AllocRequest& req, double tn, double mn) const;
    int buildLs(const AllocRequest& req, const double* x, const double* xd, double tn, double mn,
                double mu);
    bool solveLs(int rows, double* x, int* iters);
    int iterate(const AllocRequest& rq, const double* xd, double tn, double mn, double* x, int* iters,
                bool* ok);
    int prioritize(const AllocRequest& req, const double* xd, double tn, double mn, double* x,
                   int* iters, bool* ok);
    void setupCones();
    void coneStart(double* x) const;
    bool addCut(int i, double b1, double b2);
    void resetCuts(const double* x);
    bool solveCone(int rows, const double* x, double* xn, int* iters);
    void toPhysical(double* xn) const;
    void toVars(const EngineCommand& c, double* x) const;
    void fromVars(const double* x, EngineCommand* c) const;
    void evaluate(const AllocRequest& req, const EngineCommand& c, double tn, double mn,
                  AllocResult* res) const;

    const EngineConfig* m_cfg;                 // описание установки
    AllocConfig m_ac;                          // настройки
    int m_n;                                   // число переменных
    int m_vu[ENGINE_MAX_NOZZLES];              // переменная дросселя сопла
    int m_va1[ENGINE_MAX_NOZZLES];             // переменная угла 1 (-1 -- нет)
    int m_va2[ENGINE_MAX_NOZZLES];             // переменная угла 2 (-1 -- нет)
    double m_hingeMax;                         // наибольший вынос шарнира, м
    unsigned char m_thr[ALLOC_MAX_VARS];       // 1 -- переменная дросселя
    double m_sc[ALLOC_MAX_VARS];               // масштабы переменных
    double m_lb[ALLOC_MAX_VARS];               // границы вызова: пределы и окно скорости
    double m_ub[ALLOC_MAX_VARS];
    double m_cone[ENGINE_MAX_NOZZLES];         // действующий конус вызова, рад
    int m_pair[ENGINE_MAX_NOZZLES];            // пара решателя (-1 -- конус не действует)
    int m_cutN[ENGINE_MAX_NOZZLES];            // число касательных к конусу (первые -- начальные)
    double m_cutB[ENGINE_MAX_NOZZLES][ALLOC_MAX_CUTS][2]; // точки касания (d1, d2), рад
    double m_rowScale[6];                      // множители весов строк (0 -- строки нет)
    double m_regScale;                         // множитель весов регуляризации
    double m_g[6];                             // Fbu, Mbu, Mbe, Mbn, Fbe, Fbn
    double m_J[6][ALLOC_MAX_VARS];             // производные m_g
    double m_A[LSC_MAX_ROWS * LSC_MAX_COLS];   // матрица задачи (масштаб)
    double m_b[LSC_MAX_ROWS];                  // правая часть
    double m_z[ALLOC_MAX_VARS];                // переменные в масштабе
    EngineCommand m_prev;                      // предыдущая команда
    unsigned char m_on[ENGINE_MAX_NOZZLES];    // признаки работы текущего вызова
    Lsc m_ls;                                  // решатель
};

} // namespace gnc

#endif // GNCLIB_LIB_ENGINE_ALLOC_H

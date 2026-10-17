//============================================================================
// lib_eso.h -- линейные расширенные наблюдатели состояния (ESO) одного
// канала (1.6.0).
//
// Eso2 -- канал 1-го порядка dv/dt = u + f (скорость и возмущение).
//   Двойной полюс ошибки в -wo, дискретизация ТОЧНАЯ по полюсу:
//   z = exp(-wo dt),
//     e = v - z1;  z1 += dt (u + z2) + 2 (1 - z) e;  z2 += (1 - z)^2 / dt e;
//   собственные числа ошибки -- ровно z (двойное) при любом dt > 0.
//   Сброс -- с нулевым либо с заданным возмущением.
//
// Eso3 -- канал 2-го порядка theta'' = u + f (угол, скорость, возмущение).
//   Форма Люенбергера, тройной полюс в -w, один явный шаг Эйлера:
//     z1' = z2 + 3 w (y - z1),  z2' = z3 + u + 3 w^2 (y - z1),
//     z3' = w^3 (y - z1).
//   Подавать ДОСТИГНУТОЕ управление, а не запрошенное (иначе насыщенный
//   орган записывает недобор в возмущение -- предельный цикл), и шагать
//   на ФАКТИЧЕСКУЮ длину шага.
//   При w dt = 0.1 явный шаг устойчив с запасом.
// C++17, без STL, без кучи, без исключений.
//============================================================================
#ifndef GNCLIB_LIB_ESO_H
#define GNCLIB_LIB_ESO_H

#include <math.h>

namespace gnc
{

class Eso2
{
public:
    // Конструктор: полоса 1 рад/с, состояние нулевое.
    Eso2() : m_wo(1.0), m_z1(0.0), m_z2(0.0), m_uLast(0.0)
    {
    }

    // Полоса наблюдателя.
    //   wo -- рад/с
    void setup(double wo)
    {
        m_wo = wo;
    }

    // Сброс: оценка скорости v0, возмущение 0.
    //   v0 -- измеренная скорость
    void reset(double v0)
    {
        m_z1 = v0;
        m_z2 = 0.0;
    }

    // Сброс с предустановленным возмущением.
    //   v0 -- скорость; f0 -- начальное возмущение
    void resetWithDist(double v0, double f0)
    {
        m_z1 = v0;
        m_z2 = f0;
    }

    // Такт наблюдателя.
    //   vMeas -- измеренная скорость; uApplied -- достигнутое управление
    //   dt -- такт, с (> 0)
    void update(double vMeas, double uApplied, double dt)
    {
        double q = 1.0 - exp(-m_wo * dt);
        double e = vMeas - m_z1;
        m_uLast = uApplied;
        m_z1 += dt * (uApplied + m_z2) + 2.0 * q * e;
        m_z2 += q * q / dt * e;
    }

    // Оценка возмущения f.
    double dist() const
    {
        return m_z2;
    }

    // Оценка скорости z1.
    double est() const
    {
        return m_z1;
    }

    // Оценка полной производной dv/dt = u + f.
    double accel() const
    {
        return m_uLast + m_z2;
    }

private:
    double m_wo;    // полоса, рад/с
    double m_z1;    // оценка скорости
    double m_z2;    // оценка возмущения
    double m_uLast; // последнее достигнутое управление
};

class Eso3
{
public:
    // Наблюдатель без полюса; Setup() делает его рабочим.
    Eso3()
    {
    }

    // Настройка и сброс.
    //   wObs -- тройной полюс, рад/с
    // Возвращает false при неположительном полюсе.
    bool Setup(double wObs)
    {
        m_ready = wObs > 0.0;
        m_w = m_ready ? wObs : 0.0;
        m_z1 = 0.0;
        m_z2 = 0.0;
        m_z3 = 0.0;
        return m_ready;
    }

    // Сброс оценок.
    //   angle -- измеренный угол, рад
    void Reset(double angle)
    {
        m_z1 = angle;
        m_z2 = 0.0;
        m_z3 = 0.0;
    }

    // Шаг.
    //   angle -- ИЗМЕРЕННЫЙ угол, рад
    //   uReal -- ДОСТИГНУТОЕ управляющее ускорение, рад/с^2
    //   dt -- ФАКТИЧЕСКАЯ длина шага, с
    void Update(double angle, double uReal, double dt)
    {
        if (!m_ready || dt <= 0.0)
        {
            return;
        }
        double e = angle - m_z1;
        double w = m_w;
        double d1 = m_z2 + 3.0 * w * e;
        double d2 = m_z3 + uReal + 3.0 * w * w * e;
        double d3 = w * w * w * e;
        m_z1 = m_z1 + dt * d1;
        m_z2 = m_z2 + dt * d2;
        m_z3 = m_z3 + dt * d3;
    }

    // Оценка возмущающего ускорения, рад/с^2.
    double Disturbance() const
    {
        return m_z3;
    }

    // Оценки угла и скорости, рад и рад/с (для проверок и диагностики).
    double Angle() const
    {
        return m_z1;
    }

    double Rate() const
    {
        return m_z2;
    }

private:
    double m_w = 0.0;
    double m_z1 = 0.0;
    double m_z2 = 0.0;
    double m_z3 = 0.0;
    bool m_ready = false;
};

} // namespace gnc

#endif // GNCLIB_LIB_ESO_H

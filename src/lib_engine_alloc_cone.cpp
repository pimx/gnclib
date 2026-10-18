//============================================================================
// lib_engine_alloc_cone.cpp -- распределитель силовой установки (1.7.0):
// конус сопла с двумя приводами (см. lib_engine_alloc.h).
//
// Область конуса cos d1 cos d2 >= cos c выпукла (подмножество уровня
// выпуклой функции -ln cos d1 - ln cos d2), поэтому её касательные --
// внешние приближения: касательная в точке b границы
//   tan b1 (d1 - b1) + tan b2 (d2 - b2) <= 0
// не отсекает допустимых углов. Задача наименьших квадратов решается с
// касательными как с линейными ограничениями пары углов (lib_lsc.h) --
// подзадача линеаризованной модели остаётся выпуклой. Начальный набор
// вызова -- восемь касательных через 45 град; решение, вышедшее из конуса,
// добавляет касательную в своём радиальном проекте на границу (отсекающая
// плоскость Келли) и задача решается снова, с тёплым стартом: прежнее
// решение, сопло перенесено на новую касательную по отрезку от точки
// линеаризации. Сопло на границе в точке линеаризации получает касательную
// в этой точке (продолжение по границе; выход внутрь -- по знаку
// множителя). Остаток выхода за конус (меньше 1e-9 угла) снимается
// радиальным масштабированием.
//============================================================================
#include "lib_engine_alloc.h"
#include "lib_util.h"

#include <math.h>
#include <string.h>

namespace gnc
{

// Полный угол отклонения сопла с двумя приводами.
//   d1, d2 -- углы, рад
static double ConeAngle(double d1, double d2)
{
    double c = cos(d1) * cos(d2);
    c = c > 1.0 ? 1.0 : c;
    return acos(c);
}

void EngineAllocator::setupCones()
{
    const EngineConfig& c = *m_cfg;
    for (int i = 0; i < c.count; i++)
    {
        m_pair[i] = -1;
        m_cutN[i] = 0;
        int v1 = m_va1[i];
        int v2 = m_va2[i];
        if (v2 < 0 || m_on[i] == 0)
        {
            continue;
        }
        double l1 = fmax(fabs(m_lb[v1]), fabs(m_ub[v1]));
        double l2 = fmax(fabs(m_lb[v2]), fabs(m_ub[v2]));
        if (ConeInside(l1, l2, m_cone[i]))
        {
            continue;
        }
        m_pair[i] = v1;
        for (int k = 0; k < ALLOC_SEED_CUTS; k++)
        {
            double p = k * (2.0 * PI / ALLOC_SEED_CUTS);
            double b1 = 2.0 * m_cone[i] * cos(p);
            double b2 = 2.0 * m_cone[i] * sin(p);
            ConeScale(&b1, &b2, m_cone[i]);
            m_cutB[i][k][0] = b1;
            m_cutB[i][k][1] = b2;
        }
    }
}

// Набор касательных точки линеаризации: начальные (через 45 град) и
// касательная в точке x сопла, стоящего на границе. Касательные прежних
// линеаризаций не переносятся: густой многоугольник вдали от решения
// заставляет решатель обходить его вершины по одной.
void EngineAllocator::resetCuts(const double* x)
{
    const EngineConfig& c = *m_cfg;
    double zl[ALLOC_MAX_VARS];
    double zu[ALLOC_MAX_VARS];
    for (int v = 0; v < m_n; v++)
    {
        zl[v] = m_lb[v] / m_sc[v];
        zu[v] = m_ub[v] / m_sc[v];
    }
    m_ls.setup(m_n, zl, zu);
    for (int i = 0; i < c.count; i++)
    {
        if (m_pair[i] < 0)
        {
            continue;
        }
        m_ls.pair(m_va1[i], m_va2[i]);
        m_cutN[i] = 0;
        for (int k = 0; k < ALLOC_SEED_CUTS; k++)
        {
            addCut(i, m_cutB[i][k][0], m_cutB[i][k][1]);
        }
        double d1 = x[m_va1[i]];
        double d2 = x[m_va2[i]];
        if (ConeAngle(d1, d2) >= m_cone[i] * (1.0 - 1e-13))
        {
            addCut(i, d1, d2);
        }
    }
}

// Начальная точка вне конуса (предыдущая команда при суженном программном
// пределе) переносится на границу конуса по отрезку к точке окна,
// ближайшей к нулю (она в конусе: функция -ln cos d1 - ln cos d2
// разделима и растёт с |d1|, |d2|, её минимум на окне -- в этой точке;
// setBounds снимает окно, если и она вне конуса).
void EngineAllocator::coneStart(double* x) const
{
    for (int i = 0; i < m_cfg->count; i++)
    {
        int v1 = m_va1[i];
        int v2 = m_va2[i];
        if (m_pair[i] < 0 || ConeInside(x[v1], x[v2], m_cone[i]))
        {
            continue;
        }
        double z1 = Clamp(0.0, m_lb[v1], m_ub[v1]);
        double z2 = Clamp(0.0, m_lb[v2], m_ub[v2]);
        double lo = 0.0;
        double hi = 1.0;
        for (int k = 0; k < 60; k++)
        {
            double t = 0.5 * (lo + hi);
            if (ConeInside(z1 + t * (x[v1] - z1), z2 + t * (x[v2] - z2), m_cone[i]))
            {
                lo = t;
            }
            else
            {
                hi = t;
            }
        }
        x[v1] = z1 + lo * (x[v1] - z1);
        x[v2] = z2 + lo * (x[v2] - z2);
    }
}

// Касательная к конусу сопла i в точке границы b. Точки ближе 1e-6 конуса
// к имеющимся не добавляются: между касательными в точках на расстоянии t
// граница выходит за них не больше чем на t^2 / (8 c) -- меньше допуска
// 1e-9 c, а почти параллельные касательные дают плохо обусловленные вершины.
bool EngineAllocator::addCut(int i, double b1, double b2)
{
    if (m_pair[i] < 0 || m_cutN[i] >= ALLOC_MAX_CUTS)
    {
        return false;
    }
    double tol = 1e-6 * m_cone[i];
    for (int k = 0; k < m_cutN[i]; k++)
    {
        if (fabs(m_cutB[i][k][0] - b1) <= tol && fabs(m_cutB[i][k][1] - b2) <= tol)
        {
            return false;
        }
    }
    double t1 = tan(b1);
    double t2 = tan(b2);
    int v1 = m_va1[i];
    int v2 = m_va2[i];
    if (!m_ls.cut(m_pair[i], t1 * m_sc[v1], t2 * m_sc[v2], t1 * b1 + t2 * b2))
    {
        return false;
    }
    m_cutB[i][m_cutN[i]][0] = b1;
    m_cutB[i][m_cutN[i]][1] = b2;
    m_cutN[i]++;
    return true;
}

bool EngineAllocator::solveCone(int rows, const double* x, double* xn, int* iters)
{
    bool good = true;
    for (int round = 0; round < 24; round++)
    {
        good = solveLs(rows, xn, iters) && good;
        bool any = false;
        for (int i = 0; i < m_cfg->count; i++)
        {
            if (m_pair[i] < 0)
            {
                continue;
            }
            int v1 = m_va1[i];
            int v2 = m_va2[i];
            if (ConeAngle(xn[v1], xn[v2]) <= m_cone[i] * (1.0 + 1e-9))
            {
                continue;
            }
            double b1 = xn[v1];
            double b2 = xn[v2];
            ConeScale(&b1, &b2, m_cone[i]);
            if (!addCut(i, b1, b2))
            {
                continue;
            }
            double t1 = tan(b1);
            double t2 = tan(b2);
            double h = t1 * b1 + t2 * b2;
            double sx = t1 * x[v1] + t2 * x[v2];
            double sn = t1 * xn[v1] + t2 * xn[v2];
            double beta = sn > sx ? (h - sx) / (sn - sx) : 0.0;
            beta = Clamp(beta, 0.0, 1.0);
            xn[v1] = x[v1] + beta * (xn[v1] - x[v1]);
            xn[v2] = x[v2] + beta * (xn[v2] - x[v2]);
            any = true;
        }
        if (!any)
        {
            break;
        }
    }
    return good;
}

void EngineAllocator::toPhysical(double* xn) const
{
    for (int i = 0; i < m_cfg->count; i++)
    {
        if (m_pair[i] >= 0)
        {
            ConeScale(&xn[m_va1[i]], &xn[m_va2[i]], m_cone[i]);
        }
    }
}

} // namespace gnc

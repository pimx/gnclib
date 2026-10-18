//============================================================================
// lib_engine_alloc_prio.cpp -- распределитель силовой установки (1.7.0):
// строгий приоритет групп при недостижимом требовании (см.
// lib_engine_alloc.h).
//
// Группы: тангаж и рыскание (строки Mbe, Mbn), крен (Mbu), тяга (Fbu,
// только в режиме ALLOC_THRUST); порядок -- по убыванию весов wPitch, wRoll,
// wThrust (при равенстве -- тангаж, крен, тяга). Группы решаются по очереди:
// уже пройденные -- со своими весами и достижимыми целями, текущая -- с
// весом, уменьшенным до 1e3 / (sqrt(gamma) w) от своего (не меньше 1e-3;
// без уменьшения, если старших нет), младшие исключены, регуляризация
// ослаблена в 100 раз (иначе она тянет к желаемой точке и текущая группа не
// доходит до своего наибольшего значения). Достигнутое значение
// текущей группы становится её целью. Затем общее решение по достижимым
// целям: старшие группы выполняются с невязкой ~ 1 / gamma, недостижимое
// остаётся на младших.
//============================================================================
#include "lib_engine_alloc.h"
#include "lib_util.h"

#include <math.h>
#include <string.h>

namespace gnc
{

// Строки групп (-1 -- нет): 0 -- тангаж и рыскание, 1 -- крен, 2 -- тяга
static const int GROUP_ROWS[3][2] = { { 2, 3 }, { 1, -1 }, { 0, -1 } };

// Порядок групп по убыванию весов (устойчиво).
//   a -- настройки; order -- результат
static void GroupOrder(const AllocConfig& a, int order[3])
{
    double w[3] = { a.wPitch, a.wRoll, a.wThrust };
    order[0] = 0;
    order[1] = 1;
    order[2] = 2;
    for (int i = 1; i < 3; i++)
    {
        int g = order[i];
        int j = i - 1;
        while (j >= 0 && w[order[j]] < w[g])
        {
            order[j + 1] = order[j];
            j--;
        }
        order[j + 1] = g;
    }
}

// Цель строки требования.
//   q -- требование; r -- строка (0 -- тяга, 1..3 -- моменты)
static double Target(const AllocRequest& q, int r)
{
    return r == 0 ? q.thrust : (r == 1 ? q.moment.bu : (r == 2 ? q.moment.be : q.moment.bn));
}

// Новая цель строки требования.
//   q -- требование; r -- строка; v -- значение
static void SetTarget(AllocRequest* q, int r, double v)
{
    if (r == 0)
    {
        q->thrust = v;
    }
    else if (r == 1)
    {
        q->moment.bu = v;
    }
    else if (r == 2)
    {
        q->moment.be = v;
    }
    else
    {
        q->moment.bn = v;
    }
}

// Вес группы.
//   a -- настройки; g -- группа
static double GroupWeight(const AllocConfig& a, int g)
{
    return g == 0 ? a.wPitch : (g == 1 ? a.wRoll : a.wThrust);
}

int EngineAllocator::prioritize(const AllocRequest& req, const double* xd, double tn, double mn,
                                double* x, int* iters, bool* ok)
{
    AllocRequest rq = req;
    int order[3];
    GroupOrder(m_ac, order);
    bool done[3] = { false, false, false };
    int higher = 0;
    int k = 0;
    for (int li = 0; li < 3; li++)
    {
        int g = order[li];
        if (g == 2 && req.mode != ALLOC_THRUST)
        {
            continue;
        }
        linearize(rq, x, true);
        double nrm = g == 2 ? tn : mn;
        bool met = true;
        for (int j = 0; j < 2; j++)
        {
            int r = GROUP_ROWS[g][j];
            met = met && (r < 0 || fabs(m_g[r] - Target(rq, r)) / nrm < m_ac.feasTol);
        }
        if (!met)
        {
            double eps = 1e3 / (sqrt(m_ac.gamma) * GroupWeight(m_ac, g));
            eps = higher == 0 ? 1.0 : Clamp(eps, 1e-3, 1.0);
            for (int r = 0; r < 4; r++)
            {
                int h = r == 0 ? 2 : (r == 1 ? 1 : 0);
                m_rowScale[r] = h == g ? eps : (done[h] ? 1.0 : 0.0);
            }
            m_regScale = 1e-2;
            k += iterate(rq, xd, tn, mn, x, iters, ok);
            m_regScale = 1.0;
            linearize(rq, x, true);
            for (int j = 0; j < 2; j++)
            {
                int r = GROUP_ROWS[g][j];
                if (r >= 0)
                {
                    SetTarget(&rq, r, m_g[r]);
                }
            }
        }
        done[g] = true;
        higher++;
    }
    for (int r = 0; r < 6; r++)
    {
        m_rowScale[r] = 1.0;
    }
    k += iterate(rq, xd, tn, mn, x, iters, ok);
    return k;
}

} // namespace gnc

//============================================================================
// lib_engine_alloc_ls.cpp -- распределитель силовой установки (1.7.0):
// линеаризация номинальной модели, построение и решение задачи
// наименьших квадратов с границами, удержание углов в конусе, оценка
// результата (см. lib_engine_alloc.h).
//============================================================================
#include "lib_engine_alloc.h"
#include "lib_frame_ops.h"
#include "lib_util.h"

#include <math.h>
#include <string.h>

namespace gnc
{

// Столбец производных: сила и момент по одной переменной.
//   J -- матрица производных; v -- переменная; f -- производная силы;
//   m -- производная момента
static void PutColumn(double J[6][ALLOC_MAX_VARS], int v, const CfBodyVector& f,
                      const CfBodyMoment& m)
{
    J[0][v] = f.bu;
    J[1][v] = m.bu;
    J[2][v] = m.be;
    J[3][v] = m.bn;
    J[4][v] = f.be;
    J[5][v] = f.bn;
}

// exact = false: производные направления по углам берутся без
// составляющей вдоль оси сопла (она ~ sin d и обращается в ноль при d = 0):
// при большой невязке тяги (недостижимое требование) точные шаги
// Гаусса--Ньютона по углам колеблются между пределами. Значения сил --
// точные, поэтому неподвижная точка итераций та же. Вблизи решения
// (exact = true) -- точные производные и квадратичная сходимость.
void EngineAllocator::linearize(const AllocRequest& req, const double* x, bool exact)
{
    const EngineConfig& c = *m_cfg;
    memset(m_g, 0, sizeof(m_g));
    memset(m_J, 0, sizeof(m_J));
    for (int i = 0; i < c.count; i++)
    {
        if (m_on[i] == 0)
        {
            continue;
        }
        const NozzleConfig& n = c.nozzle[i];
        double u = x[m_vu[i]];
        double d1 = m_va1[i] >= 0 ? x[m_va1[i]] : 0.0;
        double d2 = m_va2[i] >= 0 ? x[m_va2[i]] : 0.0;
        double f = NozzleThrust(c, u, req.pa, 1.0);
        double df = NozzleThrustSlope(c, u, req.pa, 1.0);
        CfBodyVector d = NozzleDirection(n, d1, d2);
        CfBodyVector dd1;
        CfBodyVector dd2;
        NozzleDirectionDerivs(n, d1, d2, &dd1, &dd2);
        if (!exact)
        {
            dd1 = sub(dd1, scale(n.axis, dot(dd1, n.axis)));
            dd2 = sub(dd2, scale(n.axis, dot(dd2, n.axis)));
        }
        CfBodyVector rho = sub(n.hinge, req.point);
        CfBodyVector fv = scale(d, f);
        CfBodyMoment mv = Torque(rho, fv);
        m_g[0] += fv.bu;
        m_g[1] += mv.bu;
        m_g[2] += mv.be;
        m_g[3] += mv.bn;
        m_g[4] += fv.be;
        m_g[5] += fv.bn;
        CfBodyVector fu = scale(d, df);
        PutColumn(m_J, m_vu[i], fu, Torque(rho, fu));
        if (m_va1[i] >= 0)
        {
            CfBodyVector f1 = scale(dd1, f);
            PutColumn(m_J, m_va1[i], f1, Torque(rho, f1));
        }
        if (m_va2[i] >= 0)
        {
            CfBodyVector f2 = scale(dd2, f);
            PutColumn(m_J, m_va2[i], f2, Torque(rho, f2));
        }
    }
}

double EngineAllocator::residual(const AllocRequest& req, double tn, double mn) const
{
    bool thr = req.mode == ALLOC_THRUST && m_rowScale[0] > 0.0;
    double e = thr ? fabs(m_g[0] - req.thrust) / tn : 0.0;
    double r = m_rowScale[1] > 0.0 ? fabs(m_g[1] - req.moment.bu) / mn : 0.0;
    double y = m_rowScale[2] > 0.0 ? fabs(m_g[2] - req.moment.be) / mn : 0.0;
    double z = m_rowScale[3] > 0.0 ? fabs(m_g[3] - req.moment.bn) / mn : 0.0;
    e = r > e ? r : e;
    e = y > e ? y : e;
    return z > e ? z : e;
}

// mu > 0 -- к строкам регуляризации добавляется mu^2 (z - zk)^2, zk --
// текущая точка: w^2 (z - zd)^2 + mu^2 (z - zk)^2 = (w^2 + mu^2) (z - z*)^2 +
// const, z* = (w^2 zd + mu^2 zk) / (w^2 + mu^2) -- одна строка на переменную.
int EngineAllocator::buildLs(const AllocRequest& req, const double* x, const double* xd,
                             double tn, double mn, double mu)
{
    int n = m_n;
    double want[6] = { req.thrust, req.moment.bu, req.moment.be, req.moment.bn, 0.0, 0.0 };
    double sg = sqrt(m_ac.gamma);
    double cr[6] = { sg * m_ac.wThrust / tn, sg * m_ac.wRoll / mn, sg * m_ac.wPitch / mn,
                     sg * m_ac.wPitch / mn, m_ac.wLateral / tn, m_ac.wLateral / tn };
    int nr = m_ac.wLateral > 0.0 ? 6 : 4;
    int row = 0;
    for (int r = req.mode == ALLOC_THROTTLE ? 1 : 0; r < nr; r++)
    {
        if (m_rowScale[r] <= 0.0)
        {
            continue;
        }
        cr[r] *= m_rowScale[r];
        double lin = want[r] - m_g[r];
        for (int v = 0; v < n; v++)
        {
            lin += m_J[r][v] * x[v];
            m_A[row * n + v] = cr[r] * m_J[r][v] * m_sc[v];
        }
        m_b[row] = cr[r] * lin;
        row++;
    }
    for (int v = 0; v < n; v++)
    {
        double w = m_regScale * (m_thr[v] ? m_ac.wThrottle : m_ac.wAngle);
        double w2 = w * w + mu * mu;
        double target = (w * w * xd[v] + mu * mu * x[v]) / w2;
        for (int k = 0; k < n; k++)
        {
            m_A[row * n + k] = 0.0;
        }
        m_A[row * n + v] = sqrt(w2);
        m_b[row] = sqrt(w2) * target / m_sc[v];
        row++;
    }
    return row;
}

bool EngineAllocator::solveLs(int rows, double* x, int* iters)
{
    for (int v = 0; v < m_n; v++)
    {
        m_z[v] = x[v] / m_sc[v];
    }
    int maxIter = 6 * m_n + 100;
    int it = m_ls.solve(m_A, m_b, rows, m_z, maxIter);
    if (it == -3)
    {
        *iters += 1;
        return false;
    }
    for (int v = 0; v < m_n; v++)
    {
        x[v] = Clamp(m_z[v] * m_sc[v], m_lb[v], m_ub[v]);
    }
    *iters += it > 0 ? it : maxIter;
    return it > 0;
}

void EngineAllocator::evaluate(const AllocRequest& req, const EngineCommand& cmd, double tn,
                               double mn, AllocResult* res) const
{
    res->mode = req.mode;
    res->command = cmd;
    EngineForcesOf(*m_cfg, cmd, m_on, req.pa, req.point, false, &res->forces);
    const EngineForces& f = res->forces;
    double su = 0.0;
    int n = 0;
    for (int i = 0; i < m_cfg->count; i++)
    {
        su += m_on[i] ? cmd.throttle[i] : 0.0;
        n += m_on[i] ? 1 : 0;
    }
    res->throttle = n > 0 ? su / (double)n : 0.0;
    res->throttleLimited = 0;
    res->errThrust = req.mode == ALLOC_THROTTLE ? 0.0 : fabs(f.axial - req.thrust) / tn;
    res->errRoll = fabs(f.moment.bu - req.moment.bu) / mn;
    double ey = fabs(f.moment.be - req.moment.be) / mn;
    double ez = fabs(f.moment.bn - req.moment.bn) / mn;
    res->errPitch = ey > ez ? ey : ez;
    res->scaleThrust = tn;
    res->scaleMoment = mn;
    double tol = m_ac.feasTol;
    res->feasible = (res->errThrust < tol && res->errRoll < tol && res->errPitch < tol) ? 1 : 0;
}

} // namespace gnc

//============================================================================
// lib_engine_alloc.cpp -- распределение тяги и момента по органам силовой
// установки (1.7.0), см. lib_engine_alloc.h. Построение и решение
// задачи наименьших квадратов -- lib_engine_alloc_ls.cpp.
//============================================================================
#include "lib_engine_alloc.h"
#include "lib_frame_ops.h"
#include "lib_util.h"

#include <math.h>
#include <string.h>

namespace gnc
{

void AllocConfigDefaults(AllocConfig* a)
{
    a->wThrust = 1.0;
    a->wRoll = 30.0;
    a->wPitch = 1000.0;
    a->wThrottle = 1.0;
    a->wAngle = 1.0;
    a->wLateral = 0.0;
    a->gamma = 1e8;
    a->iterations = 8;
    a->rateWindow = 1;
    a->dt = 0.01;
    a->feasTol = 1e-6;
}

EngineAllocator::EngineAllocator()
{
    m_cfg = 0;
    AllocConfigDefaults(&m_ac);
    m_n = 0;
    m_hingeMax = 0.0;
    memset(m_vu, 0, sizeof(m_vu));
    memset(m_va1, 0, sizeof(m_va1));
    memset(m_va2, 0, sizeof(m_va2));
    memset(m_thr, 0, sizeof(m_thr));
    memset(m_g, 0, sizeof(m_g));
    for (int r = 0; r < 6; r++)
    {
        m_rowScale[r] = 1.0;
    }
    m_regScale = 1.0;
    memset(&m_prev, 0, sizeof(m_prev));
    memset(m_on, 0, sizeof(m_on));
    memset(m_cone, 0, sizeof(m_cone));
    memset(m_pair, 0, sizeof(m_pair));
    memset(m_cutN, 0, sizeof(m_cutN));
}

bool EngineAllocator::init(const EngineConfig* cfg, const AllocConfig* ac)
{
    m_cfg = 0;
    if (cfg == 0 || ac == 0 || !EngineConfigCheck(*cfg))
    {
        return false;
    }
    m_cfg = cfg;
    m_ac = *ac;
    m_ac.wThrottle = m_ac.wThrottle < 1e-9 ? 1e-9 : m_ac.wThrottle;
    m_ac.wAngle = m_ac.wAngle < 1e-9 ? 1e-9 : m_ac.wAngle;
    m_ac.iterations = m_ac.iterations < 1 ? 1 : m_ac.iterations;
    buildMap();
    EngineCommand x;
    memset(&x, 0, sizeof(x));
    for (int i = 0; i < cfg->count; i++)
    {
        x.throttle[i] = cfg->throttleMax;
    }
    reset(x);
    return true;
}

void EngineAllocator::buildMap()
{
    const EngineConfig& c = *m_cfg;
    m_n = 0;
    m_hingeMax = 0.0;
    for (int i = 0; i < c.count; i++)
    {
        const NozzleConfig& n = c.nozzle[i];
        m_vu[i] = m_n;
        m_thr[m_n] = 1;
        m_sc[m_n] = c.throttleMax - c.throttleMin;
        m_n++;
        m_va1[i] = -1;
        m_va2[i] = -1;
        if (n.kind == NOZZLE_PLANE || n.kind == NOZZLE_FREE)
        {
            m_va1[i] = m_n;
            m_thr[m_n] = 0;
            m_sc[m_n] = n.limit1;
            m_n++;
        }
        if (n.kind == NOZZLE_FREE)
        {
            m_va2[i] = m_n;
            m_thr[m_n] = 0;
            m_sc[m_n] = n.limit2;
            m_n++;
        }
        double h = norm(n.hinge);
        m_hingeMax = h > m_hingeMax ? h : m_hingeMax;
    }
}

void EngineAllocator::reset(const EngineCommand& x)
{
    if (m_cfg == 0)
    {
        return;
    }
    m_prev = x;
    for (int i = 0; i < m_cfg->count; i++)
    {
        m_prev.throttle[i] = Clamp(x.throttle[i], m_cfg->throttleMin, m_cfg->throttleMax);
        NozzleClampAngles(m_cfg->nozzle[i], &m_prev.angle1[i], &m_prev.angle2[i]);
    }
}

// Границы переменной: положение и окно скорости вокруг предыдущей команды.
//   p -- предыдущая команда; lo, hi -- пределы положения; rate -- предел
//   скорости (<= 0 -- нет); dt -- шаг; win -- окно включено; l, u -- границы
static void Window(double p, double lo, double hi, double rate, double dt, bool win,
                   double* l, double* u)
{
    double a = lo;
    double b = hi;
    if (win && rate > 0.0)
    {
        a = p - rate * dt > lo ? p - rate * dt : lo;
        b = p + rate * dt < hi ? p + rate * dt : hi;
    }
    if (a > b)
    {
        a = Clamp(p, lo, hi);
        b = a;
    }
    *l = a;
    *u = b;
}

// Границы угла сопла (выключенное сопло -- угол закреплён).
//   v -- переменная (-1 -- нет); p -- предыдущая команда; lim -- предел;
//   rate, dt, win -- окно; on -- сопло работает; lb, ub -- массивы
static void AngleBounds(int v, double p, double lim, double rate, double dt, bool win,
                        bool on, double* lb, double* ub)
{
    if (v < 0)
    {
        return;
    }
    if (on)
    {
        Window(p, -lim, lim, rate, dt, win, &lb[v], &ub[v]);
    }
    else
    {
        lb[v] = p;
        ub[v] = p;
    }
}

// Действующие пределы углов и конус сопла с учётом программного предела.
//   n -- сопло; angleMax -- программный предел (0 -- нет); lim1, lim2, cone -- результат
static void Limits(const NozzleConfig& n, double angleMax, double* lim1, double* lim2,
                   double* cone)
{
    *lim1 = n.limit1;
    *lim2 = n.limit2;
    *cone = n.cone;
    if (angleMax > 0.0)
    {
        *lim1 = *lim1 < angleMax ? *lim1 : angleMax;
        *lim2 = *lim2 < angleMax ? *lim2 : angleMax;
        *cone = *cone < angleMax ? *cone : angleMax;
    }
}

void EngineAllocator::setBounds(const AllocRequest& req, int* limited)
{
    const EngineConfig& c = *m_cfg;
    bool win = m_ac.rateWindow != 0 && m_ac.dt > 0.0;
    bool fixedThr = req.mode == ALLOC_THROTTLE;
    double ut = Clamp(req.throttle, c.throttleMin, c.throttleMax);
    *limited = (fixedThr && ut != req.throttle) ? 1 : 0;
    for (int i = 0; i < c.count; i++)
    {
        const NozzleConfig& n = c.nozzle[i];
        bool on = m_on[i] != 0;
        int v = m_vu[i];
        double p = m_prev.throttle[i];
        m_lb[v] = p;
        m_ub[v] = p;
        if (on)
        {
            Window(p, c.throttleMin, c.throttleMax, c.throttleRate, m_ac.dt, win, &m_lb[v], &m_ub[v]);
        }
        if (on && fixedThr)
        {
            double u = Clamp(ut, m_lb[v], m_ub[v]);
            *limited = u != ut ? 1 : *limited;
            m_lb[v] = u;
            m_ub[v] = u;
        }
        double lim1 = 0.0;
        double lim2 = 0.0;
        Limits(n, req.angleMax, &lim1, &lim2, &m_cone[i]);
        AngleBounds(m_va1[i], m_prev.angle1[i], lim1, n.rate, m_ac.dt, win, on, m_lb, m_ub);
        AngleBounds(m_va2[i], m_prev.angle2[i], lim2, n.rate, m_ac.dt, win, on, m_lb, m_ub);
        int v1 = m_va1[i];
        int v2 = m_va2[i];
        if (on && v2 >= 0)
        {
            // Окно скорости мимо конуса (сузился программный предел): точка
            // окна, ближайшая к нулю, вне конуса -- окно этого сопла
            // снимается, предел положения важнее скорости (как у одного
            // привода в Window).
            double z1 = Clamp(0.0, m_lb[v1], m_ub[v1]);
            double z2 = Clamp(0.0, m_lb[v2], m_ub[v2]);
            if (!ConeInside(z1, z2, m_cone[i]))
            {
                m_lb[v1] = -lim1;
                m_ub[v1] = lim1;
                m_lb[v2] = -lim2;
                m_ub[v2] = lim2;
            }
        }
    }
}

void EngineAllocator::desired(const AllocRequest& req, int nOn, double* xd) const
{
    const EngineConfig& c = *m_cfg;
    double ud = NozzleThrottleFor(c, req.thrust / (double)nOn, req.pa);
    ud = Clamp(ud, c.throttleMin, c.throttleMax);
    for (int v = 0; v < m_n; v++)
    {
        xd[v] = 0.0;
    }
    for (int i = 0; i < c.count; i++)
    {
        bool on = m_on[i] != 0;
        int v = m_vu[i];
        xd[v] = on ? ud : m_prev.throttle[i];
        xd[v] = req.mode == ALLOC_THROTTLE ? m_lb[v] : xd[v];
        if (m_va1[i] >= 0)
        {
            xd[m_va1[i]] = on ? 0.0 : m_prev.angle1[i];
        }
        if (m_va2[i] >= 0)
        {
            xd[m_va2[i]] = on ? 0.0 : m_prev.angle2[i];
        }
    }
}

void EngineAllocator::toVars(const EngineCommand& c, double* x) const
{
    for (int i = 0; i < m_cfg->count; i++)
    {
        x[m_vu[i]] = c.throttle[i];
        if (m_va1[i] >= 0)
        {
            x[m_va1[i]] = c.angle1[i];
        }
        if (m_va2[i] >= 0)
        {
            x[m_va2[i]] = c.angle2[i];
        }
    }
}

void EngineAllocator::fromVars(const double* x, EngineCommand* c) const
{
    memset(c, 0, sizeof(*c));
    for (int i = 0; i < m_cfg->count; i++)
    {
        c->throttle[i] = x[m_vu[i]];
        c->angle1[i] = m_va1[i] >= 0 ? x[m_va1[i]] : 0.0;
        c->angle2[i] = m_va2[i] >= 0 ? x[m_va2[i]] : 0.0;
    }
}

bool EngineAllocator::allocate(const AllocRequest& req, EngineCommand* cmd, AllocResult* res)
{
    if (m_cfg == 0)
    {
        return false;
    }
    const EngineConfig& c = *m_cfg;
    int nOn = 0;
    for (int i = 0; i < c.count; i++)
    {
        m_on[i] = req.on[i] != 0 ? 1 : 0;
        nOn += m_on[i];
    }
    double tn = (double)nOn * ENGINE_G0 * c.ivac * c.mdotMax * c.throttleMax;
    double lc = norm(req.point) + m_hingeMax;
    double mn = tn * (lc > 1e-3 ? lc : 1e-3);
    if (nOn == 0)
    {
        *cmd = m_prev;
        if (res != 0)
        {
            memset(res, 0, sizeof(*res));
            evaluate(req, *cmd, 1.0, 1.0, res);
        }
        return true;
    }
    int limited = 0;
    setBounds(req, &limited);
    setupCones();
    double x[ALLOC_MAX_VARS];
    double xd[ALLOC_MAX_VARS];
    toVars(m_prev, x);
    for (int v = 0; v < m_n; v++)
    {
        x[v] = Clamp(x[v], m_lb[v], m_ub[v]);
    }
    coneStart(x);
    desired(req, nOn, xd);
    int iters = 0;
    bool ok = true;
    int k = iterate(req, xd, tn, mn, x, &iters, &ok);
    linearize(req, x, true);
    if (residual(req, tn, mn) >= m_ac.feasTol)
    {
        k += prioritize(req, xd, tn, mn, x, &iters, &ok);
    }
    fromVars(x, cmd);
    m_prev = *cmd;
    if (res != 0)
    {
        evaluate(req, *cmd, tn, mn, res);
        res->throttleLimited = limited;
        res->saturated = m_ls.activeCount();
        res->iterations = iters;
        res->outer = k;
    }
    return ok;
}

int EngineAllocator::iterate(const AllocRequest& rq, const double* xd, double tn, double mn, double* x,
                             int* iters, bool* ok)
{
    int k = 0;
    bool exact = false;
    while (k < m_ac.iterations)
    {
        linearize(rq, x, exact);
        if (!exact && residual(rq, tn, mn) < 1e-3)
        {
            exact = true;
            linearize(rq, x, exact);
        }
        double mu = (k >= m_ac.iterations - 2 && k >= 1) ? 100.0 : 0.0;
        double xn[ALLOC_MAX_VARS];
        memcpy(xn, x, sizeof(double) * (size_t)m_n);
        int rows = buildLs(rq, x, xd, tn, mn, mu);
        resetCuts(x);
        *ok = solveCone(rows, x, xn, iters) && *ok;
        toPhysical(xn);
        double dx = 0.0;
        for (int v = 0; v < m_n; v++)
        {
            double e = fabs(xn[v] - x[v]) / m_sc[v];
            dx = e > dx ? e : dx;
            x[v] = xn[v];
        }
        k++;
        if (dx < 1e-12)
        {
            break;
        }
    }
    return k;
}

} // namespace gnc

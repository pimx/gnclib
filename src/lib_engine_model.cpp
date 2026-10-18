//============================================================================
// lib_engine_model.cpp -- модель объекта силовой установки (1.7.0), см.
// lib_engine_model.h.
//============================================================================
#include "lib_engine_model.h"
#include "lib_util.h"

#include <math.h>
#include <string.h>

namespace gnc
{

// Число шагов запаздывания (округление), -1 -- вне предела.
//   delay -- запаздывание, с; dt -- шаг, с
static int DelaySteps(double delay, double dt)
{
    double n = floor(delay / dt + 0.5);
    if (n < 0.0 || n >= (double)ENGINE_DELAY_MAX)
    {
        return -1;
    }
    return (int)n;
}

// Доля звена 1-го порядка за шаг (1 -- звена нет).
//   tau -- постоянная времени, с; dt -- шаг, с
static double LagShare(double tau, double dt)
{
    if (tau <= 0.0)
    {
        return 1.0;
    }
    return 1.0 - exp(-dt / tau);
}

EngineModel::EngineModel()
{
    m_cfg = 0;
    m_dt = 0.0;
    m_nThr = 0;
    m_aThr = 1.0;
    m_head = 0;
    memset(m_nSrv, 0, sizeof(m_nSrv));
    memset(m_aSrv, 0, sizeof(m_aSrv));
    memset(m_ring, 0, sizeof(m_ring));
    memset(&m_act, 0, sizeof(m_act));
    memset(m_v1, 0, sizeof(m_v1));
    memset(m_v2, 0, sizeof(m_v2));
    memset(m_level, 0, sizeof(m_level));
    memset(m_starting, 0, sizeof(m_starting));
    memset(m_on, 0, sizeof(m_on));
}

bool EngineModel::init(const EngineConfig* cfg, double dt)
{
    m_cfg = 0;
    if (cfg == 0 || dt <= 0.0 || !EngineConfigCheck(*cfg))
    {
        return false;
    }
    m_dt = dt;
    m_nThr = DelaySteps(cfg->throttleDelay, dt);
    m_aThr = LagShare(cfg->throttleTau, dt);
    if (m_nThr < 0)
    {
        return false;
    }
    for (int i = 0; i < cfg->count; i++)
    {
        m_nSrv[i] = DelaySteps(cfg->nozzle[i].delay, dt);
        m_aSrv[i] = LagShare(cfg->nozzle[i].tau, dt);
        if (m_nSrv[i] < 0)
        {
            return false;
        }
    }
    m_cfg = cfg;
    reset(cfg->throttleMax);
    return true;
}

void EngineModel::reset(double throttle)
{
    if (m_cfg == 0)
    {
        return;
    }
    EngineCommand x;
    memset(&x, 0, sizeof(x));
    for (int i = 0; i < m_cfg->count; i++)
    {
        x.throttle[i] = throttle;
    }
    reset(x);
}

void EngineModel::reset(const EngineCommand& x)
{
    if (m_cfg == 0)
    {
        return;
    }
    m_act = x;
    for (int i = 0; i < m_cfg->count; i++)
    {
        m_act.throttle[i] = Clamp(x.throttle[i], m_cfg->throttleMin, m_cfg->throttleMax);
        NozzleClampAngles(m_cfg->nozzle[i], &m_act.angle1[i], &m_act.angle2[i]);
        m_on[i] = 1;
        m_v1[i] = 0.0;
        m_v2[i] = 0.0;
        m_level[i] = 0.0;
        m_starting[i] = 0;
    }
    for (int k = 0; k < ENGINE_DELAY_MAX; k++)
    {
        m_ring[k] = m_act;
    }
    m_head = 0;
}

void EngineModel::setOn(int i, bool on, bool immediate)
{
    if (m_cfg == 0 || i < 0 || i >= m_cfg->count || (m_on[i] != 0) == on)
    {
        return;
    }
    if (on)
    {
        m_on[i] = 1;
        if (m_cfg->startRate > 0.0 && !immediate)
        {
            m_act.throttle[i] = m_level[i];
            m_starting[i] = 1;
        }
        m_level[i] = 0.0;
        return;
    }
    m_on[i] = 0;
    m_starting[i] = 0;
    m_level[i] = (m_cfg->shutdownRate > 0.0 && !immediate) ? m_act.throttle[i] : 0.0;
}

const EngineCommand& EngineModel::delayed(int n) const
{
    int k = (m_head - n + ENGINE_DELAY_MAX) % ENGINE_DELAY_MAX;
    return m_ring[k];
}

double EngineModel::follow(double x, double target, double rate, double alpha) const
{
    double d = target - x;
    if (alpha < 1.0)
    {
        d = d * alpha;
    }
    if (rate > 0.0)
    {
        double r = rate * m_dt;
        if (d > r)
        {
            return x + r;
        }
        if (d < -r)
        {
            return x - r;
        }
    }
    return alpha < 1.0 ? x + d : target;
}

void EngineModel::servo2(const NozzleConfig& n, double* x, double* v, double target,
                         double limit) const
{
    int ns = (int)ceil(n.wn * m_dt / 0.05);
    ns = ns < 1 ? 1 : ns;
    double h = m_dt / (double)ns;
    double w2 = n.wn * n.wn;
    double c = 2.0 * n.zeta * n.wn;
    for (int k = 0; k < ns; k++)
    {
        *v += h * (w2 * (target - *x) - c * *v);
        if (n.rate > 0.0)
        {
            *v = Limit(*v, n.rate);
        }
        *x += h * *v;
        if (*x > limit || *x < -limit)
        {
            *x = Limit(*x, limit);
            *v = (*v > 0.0) == (*x > 0.0) ? 0.0 : *v;
        }
    }
}

void EngineModel::servoStep(int i, const EngineCommand& cs)
{
    const NozzleConfig& n = m_cfg->nozzle[i];
    double t1 = cs.angle1[i];
    double t2 = cs.angle2[i];
    NozzleClampAngles(n, &t1, &t2);
    if (n.wn <= 0.0)
    {
        m_act.angle1[i] = follow(m_act.angle1[i], t1, n.rate, m_aSrv[i]);
        if (n.kind == NOZZLE_FREE)
        {
            m_act.angle2[i] = follow(m_act.angle2[i], t2, n.rate, m_aSrv[i]);
            ConeScale(&m_act.angle1[i], &m_act.angle2[i], n.cone);
        }
        return;
    }
    servo2(n, &m_act.angle1[i], &m_v1[i], t1, n.limit1);
    if (n.kind != NOZZLE_FREE)
    {
        return;
    }
    servo2(n, &m_act.angle2[i], &m_v2[i], t2, n.limit2);
    double* a1 = &m_act.angle1[i];
    double* a2 = &m_act.angle2[i];
    if (ConeScale(a1, a2, n.cone))
    {
        double r = sqrt(*a1 * *a1 + *a2 * *a2);
        double vr = r > 0.0 ? (m_v1[i] * *a1 + m_v2[i] * *a2) / r : 0.0;
        if (vr > 0.0)
        {
            m_v1[i] -= vr * *a1 / r;
            m_v2[i] -= vr * *a2 / r;
        }
    }
}

void EngineModel::throttleStep(int i, const EngineCommand& ct)
{
    const EngineConfig& c = *m_cfg;
    double u = Clamp(ct.throttle[i], c.throttleMin, c.throttleMax);
    if (m_on[i] == 0)
    {
        m_level[i] = m_level[i] - c.shutdownRate * m_dt;
        m_level[i] = m_level[i] > 0.0 ? m_level[i] : 0.0;
    }
    if (m_starting[i] != 0)
    {
        double x = m_act.throttle[i] + c.startRate * m_dt;
        m_starting[i] = x < u ? 1 : 0;
        m_act.throttle[i] = x < u ? x : u;
        return;
    }
    m_act.throttle[i] = follow(m_act.throttle[i], u, c.throttleRate, m_aThr);
}

void EngineModel::step(const EngineCommand& cmd)
{
    if (m_cfg == 0)
    {
        return;
    }
    m_head = (m_head + 1) % ENGINE_DELAY_MAX;
    m_ring[m_head] = cmd;
    const EngineCommand& ct = delayed(m_nThr);
    for (int i = 0; i < m_cfg->count; i++)
    {
        throttleStep(i, ct);
        if (m_cfg->nozzle[i].kind != NOZZLE_FIXED)
        {
            servoStep(i, delayed(m_nSrv[i]));
        }
    }
}

void EngineModel::forces(double pa, const CfBodyVector& point, EngineForces* out) const
{
    if (m_cfg == 0)
    {
        memset(out, 0, sizeof(*out));
        return;
    }
    EngineCommand x = m_act;
    unsigned char on[ENGINE_MAX_NOZZLES];
    for (int i = 0; i < m_cfg->count; i++)
    {
        on[i] = (m_on[i] != 0 || m_level[i] > 0.0) ? 1 : 0;
        x.throttle[i] = throttleOut(i);
    }
    EngineForcesOf(*m_cfg, x, on, pa, point, true, out);
}

} // namespace gnc

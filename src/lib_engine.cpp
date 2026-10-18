//============================================================================
// lib_engine.cpp -- силовая установка из нескольких сопел (1.7.0):
// закон тяги, направления тяги, суммарные сила и момент (см. lib_engine.h).
//============================================================================
#include "lib_engine.h"
#include "lib_frame_ops.h"
#include "lib_util.h"

#include <math.h>
#include <stdio.h>

namespace gnc
{

void NozzleInit(NozzleConfig* n)
{
    n->kind = NOZZLE_FIXED;
    n->hinge = mkBody(0.0, 0.0, 0.0);
    n->axis = mkBody(1.0, 0.0, 0.0);
    n->e1 = mkBody(0.0, 1.0, 0.0);
    n->e2 = mkBody(0.0, 0.0, 1.0);
    n->limit1 = 0.0;
    n->limit2 = 0.0;
    n->cone = 0.0;
    n->rate = 0.0;
    n->delay = 0.0;
    n->tau = 0.0;
    n->wn = 0.0;
    n->zeta = 0.0;
    n->thrustFactor = 1.0;
    n->misalign = mkBody(0.0, 0.0, 0.0);
}

void EngineInit(EngineConfig* c)
{
    c->count = 0;
    for (int i = 0; i < ENGINE_MAX_NOZZLES; i++)
    {
        NozzleInit(&c->nozzle[i]);
    }
    c->ivac = 0.0;
    c->mdotMax = 0.0;
    c->exitArea = 0.0;
    c->throttleMin = 0.0;
    c->throttleMax = 0.0;
    c->throttleRate = 0.0;
    c->throttleDelay = 0.0;
    c->throttleTau = 0.0;
    c->lossLaw = ENGINE_LOSS_CONSTANT;
    c->startRate = 0.0;
    c->shutdownRate = 0.0;
}

// Проверка единичности вектора с диагностикой.
//   v -- вектор; what -- имя для диагностики; i -- номер сопла
static bool CheckUnit(const CfBodyVector& v, const char* what, int i)
{
    double q = norm(v);
    if (fabs(q - 1.0) > 1e-9)
    {
        fprintf(stderr, "engine: nozzle %d: %s is not unit (|v| = %.12g)\n", i, what, q);
        return false;
    }
    return true;
}

// Проверка перпендикулярности двух векторов с диагностикой.
//   a, b -- векторы; what -- имя пары; i -- номер сопла
static bool CheckOrtho(const CfBodyVector& a, const CfBodyVector& b, const char* what, int i)
{
    double c = dot(a, b);
    if (fabs(c) > 1e-9)
    {
        fprintf(stderr, "engine: nozzle %d: %s not orthogonal (dot = %.12g)\n", i, what, c);
        return false;
    }
    return true;
}

// Проверка одного сопла.
//   n -- сопло; i -- номер для диагностики
static bool CheckNozzle(const NozzleConfig& n, int i)
{
    bool ok = CheckUnit(n.axis, "axis", i);
    if (n.kind == NOZZLE_FIXED)
    {
        return ok;
    }
    if (n.kind != NOZZLE_PLANE && n.kind != NOZZLE_FREE)
    {
        fprintf(stderr, "engine: nozzle %d: unknown kind %d\n", i, n.kind);
        return false;
    }
    ok = CheckUnit(n.e1, "e1", i) && ok;
    ok = CheckOrtho(n.axis, n.e1, "axis, e1", i) && ok;
    if (n.limit1 <= 0.0 || n.delay < 0.0 || n.tau < 0.0 || n.wn < 0.0)
    {
        fprintf(stderr, "engine: nozzle %d: limit1 must be > 0, delay, tau and wn >= 0\n", i);
        ok = false;
    }
    if (n.wn > 0.0 && n.zeta <= 0.0)
    {
        fprintf(stderr, "engine: nozzle %d: zeta must be > 0 with wn > 0\n", i);
        ok = false;
    }
    if (n.kind == NOZZLE_FREE)
    {
        ok = CheckUnit(n.e2, "e2", i) && ok;
        ok = CheckOrtho(n.axis, n.e2, "axis, e2", i) && ok;
        ok = CheckOrtho(n.e1, n.e2, "e1, e2", i) && ok;
        if (n.limit2 <= 0.0 || n.cone <= 0.0)
        {
            fprintf(stderr, "engine: nozzle %d: limit2 and cone must be > 0\n", i);
            ok = false;
        }
    }
    return ok;
}

bool EngineConfigCheck(const EngineConfig& c)
{
    bool ok = true;
    if (c.count < 1 || c.count > ENGINE_MAX_NOZZLES)
    {
        fprintf(stderr, "engine: nozzle count %d out of 1..%d\n", c.count, (int)ENGINE_MAX_NOZZLES);
        return false;
    }
    for (int i = 0; i < c.count; i++)
    {
        ok = CheckNozzle(c.nozzle[i], i) && ok;
    }
    if (c.ivac <= 0.0 || c.mdotMax <= 0.0 || c.exitArea < 0.0)
    {
        fprintf(stderr, "engine: ivac and mdotMax must be > 0, exitArea >= 0\n");
        ok = false;
    }
    if (c.throttleMin < 0.0 || c.throttleMax <= c.throttleMin)
    {
        fprintf(stderr, "engine: throttle range must satisfy 0 <= min < max\n");
        ok = false;
    }
    if (c.throttleDelay < 0.0 || c.throttleTau < 0.0)
    {
        fprintf(stderr, "engine: throttle delay and tau must be >= 0\n");
        ok = false;
    }
    if (c.lossLaw != ENGINE_LOSS_CONSTANT && c.lossLaw != ENGINE_LOSS_THROTTLE)
    {
        fprintf(stderr, "engine: unknown loss law %d\n", c.lossLaw);
        ok = false;
    }
    if (c.startRate < 0.0 || c.shutdownRate < 0.0)
    {
        fprintf(stderr, "engine: start and shutdown rates must be >= 0\n");
        ok = false;
    }
    return ok;
}

double NozzleThrust(const EngineConfig& c, double u, double pa, double factor)
{
    double loss = c.lossLaw == ENGINE_LOSS_THROTTLE ? u * pa * c.exitArea : pa * c.exitArea;
    double f = factor * ENGINE_G0 * c.ivac * u * c.mdotMax - loss;
    return f > 0.0 ? f : 0.0;
}

double NozzleThrottleFor(const EngineConfig& c, double f, double pa)
{
    double k = ENGINE_G0 * c.ivac * c.mdotMax;
    if (c.lossLaw == ENGINE_LOSS_THROTTLE)
    {
        double kn = k - pa * c.exitArea;
        return kn > 0.0 ? f / kn : 0.0;
    }
    return (f + pa * c.exitArea) / k;
}

double NozzleThrustSlope(const EngineConfig& c, double u, double pa, double factor)
{
    if (NozzleThrust(c, u, pa, factor) <= 0.0)
    {
        return 0.0;
    }
    double k = factor * ENGINE_G0 * c.ivac * c.mdotMax;
    return c.lossLaw == ENGINE_LOSS_THROTTLE ? k - pa * c.exitArea : k;
}

CfBodyVector NozzleDirection(const NozzleConfig& n, double d1, double d2)
{
    if (n.kind == NOZZLE_FIXED)
    {
        return n.axis;
    }
    CfBodyVector p = add(scale(n.axis, cos(d1)), scale(n.e1, sin(d1)));
    if (n.kind == NOZZLE_PLANE)
    {
        return p;
    }
    return add(scale(p, cos(d2)), scale(n.e2, sin(d2)));
}

void NozzleDirectionDerivs(const NozzleConfig& n, double d1, double d2,
                           CfBodyVector* dd1, CfBodyVector* dd2)
{
    *dd1 = mkBody(0.0, 0.0, 0.0);
    *dd2 = mkBody(0.0, 0.0, 0.0);
    if (n.kind == NOZZLE_FIXED)
    {
        return;
    }
    CfBodyVector p = add(scale(n.axis, cos(d1)), scale(n.e1, sin(d1)));
    CfBodyVector dp = add(scale(n.axis, -sin(d1)), scale(n.e1, cos(d1)));
    if (n.kind == NOZZLE_PLANE)
    {
        *dd1 = dp;
        return;
    }
    *dd1 = scale(dp, cos(d2));
    *dd2 = add(scale(p, -sin(d2)), scale(n.e2, cos(d2)));
}

CfBodyVector NozzleMisaligned(const NozzleConfig& n, const CfBodyVector& d)
{
    if (norm(n.misalign) == 0.0)
    {
        return d;
    }
    return unit(add(d, cross(n.misalign, d)));
}

bool ConeInside(double d1, double d2, double cone)
{
    return cos(d1) * cos(d2) >= cos(cone) - 1e-15;
}

bool ConeScale(double* d1, double* d2, double cone)
{
    if (ConeInside(*d1, *d2, cone))
    {
        return false;
    }
    double lo = 0.0;
    double hi = 1.0;
    for (int k = 0; k < 60; k++)
    {
        double s = 0.5 * (lo + hi);
        if (ConeInside(s * *d1, s * *d2, cone))
        {
            lo = s;
        }
        else
        {
            hi = s;
        }
    }
    *d1 = lo * *d1;
    *d2 = lo * *d2;
    return true;
}

void NozzleClampAngles(const NozzleConfig& n, double* d1, double* d2)
{
    if (n.kind == NOZZLE_FIXED)
    {
        *d1 = 0.0;
        *d2 = 0.0;
        return;
    }
    *d1 = Limit(*d1, n.limit1);
    if (n.kind == NOZZLE_PLANE)
    {
        *d2 = 0.0;
        return;
    }
    *d2 = Limit(*d2, n.limit2);
    ConeScale(d1, d2, n.cone);
}

// Производные величины вектора тяги: модуль, составляющие, направление.
//   out -- результат (force задана)
static void ForceSummary(EngineForces* out)
{
    const CfBodyVector& f = out->force;
    out->norm = norm(f);
    out->axial = f.bu;
    out->lateral = sqrt(f.be * f.be + f.bn * f.bn);
    out->tilt = atan2(out->lateral, f.bu);
    out->azimuth = out->lateral > 0.0 ? atan2(f.bn, f.be) : 0.0;
}

void EngineForcesOf(const EngineConfig& c, const EngineCommand& x,
                    const unsigned char* on, double pa, const CfBodyVector& point,
                    bool dispersion, EngineForces* out)
{
    out->force = mkBody(0.0, 0.0, 0.0);
    out->point = point;
    out->moment = mkMom(0.0, 0.0, 0.0);
    out->mdot = 0.0;
    for (int i = 0; i < c.count; i++)
    {
        out->thrust[i] = 0.0;
        if (on != 0 && on[i] == 0)
        {
            continue;
        }
        const NozzleConfig& n = c.nozzle[i];
        double kf = dispersion ? n.thrustFactor : 1.0;
        double f = NozzleThrust(c, x.throttle[i], pa, kf);
        CfBodyVector d = NozzleDirection(n, x.angle1[i], x.angle2[i]);
        if (dispersion)
        {
            d = NozzleMisaligned(n, d);
        }
        CfBodyVector fv = scale(d, f);
        out->thrust[i] = f;
        out->force = add(out->force, fv);
        out->moment = Add(out->moment, Torque(sub(n.hinge, point), fv));
        out->mdot += x.throttle[i] * c.mdotMax;
    }
    ForceSummary(out);
}

} // namespace gnc

// test_engine_modes.cpp -- gnclib 1.7.0: thrust law with the back-pressure
// loss proportional to the throttle, force summary (vector, tilt, azimuth,
// lateral force), second-order servo, start and shutdown of a nozzle,
// allocation with a given throttle and with a program angle limit.
// Neutral configuration: a central two-servo nozzle and six outer nozzles
// with two servos each.
// Build: g++ -std=c++17 -O2 -Iinclude tests/test_engine_modes.cpp
//        build/libgnc.a -o build/test_engine_modes
#include "lib_engine.h"
#include "lib_engine_alloc.h"
#include "lib_engine_model.h"
#include "lib_frame_ops.h"
#include "lib_util.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

using namespace gnc;

static int g_bad = 0;
static EngineConfig g_cfg;
static unsigned long long g_seed = 777;

// Check of a condition with a print.
static void Check(const char* name, bool ok)
{
    printf("%-64s %s\n", name, ok ? "OK" : "FAIL");
    g_bad += ok ? 0 : 1;
}

// Uniform random number in [a, b] (reproducible LCG).
static double Uni(double a, double b)
{
    g_seed = g_seed * 6364136223846793005ULL + 1442695040888963407ULL;
    return a + (b - a) * ((double)(g_seed >> 11) / 9007199254740992.0);
}

// Neutral configuration: centre + 6 outer two-servo nozzles, loss law
// proportional to the throttle, second-order servos.
static void MakeConfig(EngineConfig* c)
{
    EngineInit(c);
    c->count = 7;
    c->ivac = 340.0;
    c->mdotMax = 150.0;
    c->exitArea = 4.0;
    c->throttleMin = 0.4;
    c->throttleMax = 1.0;
    c->throttleRate = 1.0;
    c->lossLaw = ENGINE_LOSS_THROTTLE;
    for (int i = 0; i < c->count; i++)
    {
        NozzleConfig& n = c->nozzle[i];
        n.kind = NOZZLE_FREE;
        n.limit1 = Deg(7.0);
        n.limit2 = Deg(7.0);
        n.cone = Deg(7.0);
        n.rate = Deg(25.0);
        n.wn = 30.0;
        n.zeta = 0.7;
        double phi = Deg(60.0 * (i - 1));
        n.hinge = i == 0 ? mkBody(0.0, 0.0, 0.0) : mkBody(0.0, 1.2 * cos(phi), 1.2 * sin(phi));
    }
}

static void TestThrustLaw()
{
    double k = ENGINE_G0 * 340.0 * 150.0;
    double f = NozzleThrust(g_cfg, 0.6, 50000.0, 1.0);
    Check("loss proportional to throttle: F = u (k - pa Sa)", fabs(f - 0.6 * (k - 2e5)) < 1e-6);
    double u = NozzleThrottleFor(g_cfg, f, 50000.0);
    Check("inverse thrust law", fabs(u - 0.6) < 1e-14);
    double sl = NozzleThrustSlope(g_cfg, 0.6, 50000.0, 1.0);
    Check("thrust slope k - pa Sa", fabs(sl - (k - 2e5)) < 1e-6);
}

static void TestSummary()
{
    EngineCommand x;
    memset(&x, 0, sizeof(x));
    for (int i = 0; i < g_cfg.count; i++)
    {
        x.throttle[i] = 1.0;
        x.angle1[i] = Deg(3.0);
    }
    EngineForces f;
    EngineForcesOf(g_cfg, x, 0, 0.0, mkBody(10.0, 0.0, 0.0), false, &f);
    double F = 7.0 * NozzleThrust(g_cfg, 1.0, 0.0, 1.0);
    bool ok = fabs(f.norm - F) < 1e-6 && fabs(f.tilt - Deg(3.0)) < 1e-12;
    ok = ok && fabs(f.azimuth) < 1e-12 && fabs(f.lateral - F * sin(Deg(3.0))) < 1e-6;
    ok = ok && fabs(f.moment.bn + 10.0 * f.lateral) < 1e-5 && f.point.bu == 10.0;
    Check("summary: vector, tilt, azimuth, lateral force, moment", ok);
}

static void TestServo2()
{
    static EngineModel m;
    Check("model with second-order servo: init", m.init(&g_cfg, 0.01));
    m.reset(0.7);
    EngineCommand c = m.actual();
    c.angle1[1] = Deg(2.0);
    double peak = 0.0;
    double vmax = 0.0;
    double prev = 0.0;
    for (int k = 0; k < 100; k++)
    {
        m.step(c);
        double a = m.actual().angle1[1];
        peak = fmax(peak, a);
        vmax = fmax(vmax, fabs(a - prev) / 0.01);
        prev = a;
    }
    double over = peak / Deg(2.0) - 1.0;
    Check("servo2: settles on the command", fabs(prev - Deg(2.0)) < 1e-6);
    Check("servo2: overshoot of zeta 0.7 (4.6 % +- 1.5 %)", over > 0.031 && over < 0.061);
    c.angle1[1] = Deg(7.0);
    c.angle2[1] = Deg(7.0);
    for (int k = 0; k < 100; k++)
    {
        m.step(c);
        double a = m.actual().angle1[1];
        vmax = fmax(vmax, fabs(a - prev) / 0.01);
        prev = a;
    }
    const EngineCommand& a = m.actual();
    Check("servo2: rate within the limit", vmax <= Deg(25.0) * (1.0 + 1e-12));
    Check("servo2: held in the cone", ConeInside(a.angle1[1], a.angle2[1], Deg(7.0) + 1e-12));
}

static void TestStartShutdown()
{
    static EngineModel m;
    EngineConfig c = g_cfg;
    c.startRate = 0.5;
    c.shutdownRate = 1.0;
    m.init(&c, 0.01);
    m.reset(0.8);
    EngineCommand cmd = m.actual();
    m.setOn(2, false);
    EngineForces f0;
    EngineForces f1;
    m.step(cmd);
    m.forces(0.0, mkBody(0.0, 0.0, 0.0), &f0);
    for (int k = 0; k < 100; k++)
    {
        m.step(cmd);
    }
    m.forces(0.0, mkBody(0.0, 0.0, 0.0), &f1);
    bool down = fabs(m.throttleOut(2) - 0.0) < 1e-12 && f1.thrust[2] == 0.0;
    Check("shutdown: burns down at 1/s (0.79 after 0.01 s, 0 after 1.01 s)",
          fabs(f0.thrust[2] - NozzleThrust(c, 0.79, 0.0, 1.0)) < 1e-6 && down);
    m.setOn(2, true);
    m.step(cmd);
    double u1 = m.throttleOut(2);
    for (int k = 0; k < 200; k++)
    {
        m.step(cmd);
    }
    Check("start: ramps from zero at 0.5/s to the command",
          fabs(u1 - 0.005) < 1e-12 && fabs(m.throttleOut(2) - 0.8) < 1e-12);
}

// Allocation in the throttle mode of requests made by admissible commands
// with equal throttles: always feasible.
static void TestThrottleMode()
{
    static EngineAllocator al;
    AllocConfig ac;
    AllocConfigDefaults(&ac);
    ac.rateWindow = 0;
    al.init(&g_cfg, &ac);
    int bad = 0;
    double worst = 0.0;
    for (int s = 0; s < 2000; s++)
    {
        EngineCommand x;
        memset(&x, 0, sizeof(x));
        double u = Uni(0.4, 1.0);
        for (int i = 0; i < g_cfg.count; i++)
        {
            x.throttle[i] = u;
            x.angle1[i] = Uni(-Deg(7.0), Deg(7.0));
            x.angle2[i] = Uni(-Deg(7.0), Deg(7.0));
            NozzleClampAngles(g_cfg.nozzle[i], &x.angle1[i], &x.angle2[i]);
        }
        AllocRequest q;
        memset(&q, 0, sizeof(q));
        memset(q.on, 1, sizeof(q.on));
        q.mode = ALLOC_THROTTLE;
        q.throttle = u;
        q.point = mkBody(9.0, 0.0, 0.0);
        q.pa = Uni(0.0, 101325.0);
        EngineForces f;
        EngineForcesOf(g_cfg, x, 0, q.pa, q.point, false, &f);
        q.moment = f.moment;
        EngineCommand c;
        AllocResult r;
        bool ok = al.allocate(q, &c, &r);
        bool equal = true;
        for (int i = 0; i < g_cfg.count; i++)
        {
            equal = equal && c.throttle[i] == u;
        }
        worst = fmax(worst, fmax(r.errRoll, r.errPitch));
        bad += (ok && r.feasible && equal && r.throttleLimited == 0) ? 0 : 1;
    }
    printf("  worst normalized residual %.3e\n", worst);
    Check("throttle mode: 2000 feasible requests, equal throttles", bad == 0);
    AllocRequest q;
    memset(&q, 0, sizeof(q));
    memset(q.on, 1, sizeof(q.on));
    q.mode = ALLOC_THROTTLE;
    q.throttle = 1.3;
    q.point = mkBody(9.0, 0.0, 0.0);
    q.moment = mkMom(5e7, 2e5, 0.0);
    EngineCommand c;
    AllocResult r;
    al.allocate(q, &c, &r);
    printf("  roll beyond capacity: pitch residual %.2e, roll residual %.2e\n", r.errPitch, r.errRoll);
    Check("throttle mode: 1.3 limited to 1.0, roll short, pitch kept (strict, 1e-9)",
          r.throttleLimited == 1 && c.throttle[0] == 1.0 && r.errPitch < 1e-9 && r.errRoll > 1e-3);
}

static void TestAngleMax()
{
    static EngineAllocator al;
    AllocConfig ac;
    AllocConfigDefaults(&ac);
    ac.rateWindow = 0;
    al.init(&g_cfg, &ac);
    AllocRequest q;
    memset(&q, 0, sizeof(q));
    memset(q.on, 1, sizeof(q.on));
    q.mode = ALLOC_THROTTLE;
    q.throttle = 0.8;
    q.point = mkBody(9.0, 0.0, 0.0);
    q.angleMax = Deg(3.0);
    double F = 7.0 * NozzleThrust(g_cfg, 0.8, 0.0, 1.0);
    q.moment = mkMom(0.0, 9.0 * F * sin(Deg(2.5)), 0.0);
    EngineCommand c;
    AllocResult r;
    al.allocate(q, &c, &r);
    Check("program limit 3 deg: 2.5 deg equivalent met", r.feasible == 1);
    q.moment = mkMom(0.0, 9.0 * F * sin(Deg(4.0)), 0.0);
    al.allocate(q, &c, &r);
    bool inside = true;
    for (int i = 0; i < g_cfg.count; i++)
    {
        inside = inside && ConeInside(c.angle1[i], c.angle2[i], Deg(3.0) + 1e-12);
    }
    Check("program limit 3 deg: 4 deg request -- angles in the 3 deg cone, short",
          inside && r.feasible == 0 && fabs(r.forces.tilt - Deg(3.0)) < Deg(0.01));
    // The limit narrows under the rate window: the previous command at the
    // 7 deg cone, the window (25 deg/s * 0.01 s) misses the 3 deg cone.
    static EngineAllocator aw;
    ac.rateWindow = 1;
    aw.init(&g_cfg, &ac);
    EngineCommand p = c;
    for (int i = 0; i < g_cfg.count; i++)
    {
        p.angle1[i] = Deg(5.0);
        p.angle2[i] = Deg(5.0);
        NozzleClampAngles(g_cfg.nozzle[i], &p.angle1[i], &p.angle2[i]);
        p.throttle[i] = 0.8;
    }
    aw.reset(p);
    q.moment = mkMom(0.0, 9.0 * F * sin(Deg(1.0)), 0.0);
    bool good = aw.allocate(q, &c, &r);
    inside = true;
    for (int i = 0; i < g_cfg.count; i++)
    {
        inside = inside && ConeInside(c.angle1[i], c.angle2[i], Deg(3.0) + 1e-12);
    }
    Check("program limit narrowed under the rate window: in the 3 deg cone at once", good && inside);
}

int main()
{
    printf("gnclib test_engine_modes\n");
    MakeConfig(&g_cfg);
    Check("config check", EngineConfigCheck(g_cfg));
    TestThrustLaw();
    TestSummary();
    TestServo2();
    TestStartShutdown();
    TestThrottleMode();
    TestAngleMax();
    printf("%s (%d failed)\n", g_bad == 0 ? "TEST_ENGINE_MODES PASS" : "TEST_ENGINE_MODES FAIL", g_bad);
    return g_bad == 0 ? 0 : 1;
}

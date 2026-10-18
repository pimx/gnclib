// test_engine.cpp -- gnclib 1.7.0: multi-nozzle engine (thrust law, nozzle
// kinematics, plant model with delays and rate limits), bounded least
// squares (with and without constraints on pairs of variables), thrust and
// moment allocation. Neutral configuration: a central
// two-servo nozzle and four outer nozzles deflecting in tangential planes.
// Build: g++ -std=c++17 -O2 -Iinclude tests/test_engine.cpp build/libgnc.a
//        -o build/test_engine
#include "lib_bvls.h"
#include "lib_lsc.h"
#include "lib_engine.h"
#include "lib_engine_alloc.h"
#include "lib_engine_model.h"
#include "lib_frame_ops.h"
#include "lib_random.h"
#include "lib_util.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

using namespace gnc;

static int g_bad = 0;

// Check of a condition with a print.
static void Check(const char* name, bool ok)
{
    printf("%-64s %s\n", name, ok ? "OK" : "FAIL");
    g_bad += ok ? 0 : 1;
}

// Uniform random number in [a, b] (simple LCG, reproducible).
static unsigned long long g_seed = 12345;
static double Uni(double a, double b)
{
    g_seed = g_seed * 6364136223846793005ULL + 1442695040888963407ULL;
    double r = (double)(g_seed >> 11) / 9007199254740992.0;
    return a + (b - a) * r;
}

static EngineConfig g_cfg;
static void TestAllocator();

// Neutral configuration: centre + 4 outer tangential nozzles.
static void MakeConfig(EngineConfig* c)
{
    memset(c, 0, sizeof(*c));
    c->count = 5;
    c->ivac = 300.0;
    c->mdotMax = 100.0;
    c->exitArea = 0.5;
    c->throttleMin = 0.4;
    c->throttleMax = 1.0;
    c->throttleRate = 1.0;
    for (int i = 0; i < c->count; i++)
    {
        NozzleConfig& n = c->nozzle[i];
        NozzleInit(&n);
        n.limit1 = Deg(6.0);
        n.rate = Deg(10.0);
        if (i == 0)
        {
            n.kind = NOZZLE_FREE;
            n.limit2 = Deg(6.0);
            n.cone = Deg(6.0);
            continue;
        }
        double phi = Deg(90.0 * (i - 1));
        n.kind = NOZZLE_PLANE;
        n.hinge = mkBody(0.0, cos(phi), sin(phi));
        n.e1 = mkBody(0.0, -sin(phi), cos(phi));
    }
}

static void TestThrustAndKinematics()
{
    double k = ENGINE_G0 * 300.0 * 100.0;
    Check("thrust law, vacuum full throttle", NozzleThrust(g_cfg, 1.0, 0.0, 1.0) == k);
    Check("thrust law, back pressure", fabs(NozzleThrust(g_cfg, 0.5, 1e5, 1.0) - (0.5 * k - 5e4)) < 1e-9);
    Check("thrust law, floor at zero", NozzleThrust(g_cfg, 0.1, 1e6, 1.0) == 0.0);
    double eu = 0.0;
    double ec = 0.0;
    double ed = 0.0;
    const NozzleConfig& n = g_cfg.nozzle[0];
    for (int s = 0; s < 2000; s++)
    {
        double d1 = Uni(-0.2, 0.2);
        double d2 = Uni(-0.2, 0.2);
        CfBodyVector d = NozzleDirection(n, d1, d2);
        eu = fmax(eu, fabs(norm(d) - 1.0));
        ec = fmax(ec, fabs(d.bu - cos(d1) * cos(d2)));
        CfBodyVector a1;
        CfBodyVector a2;
        NozzleDirectionDerivs(n, d1, d2, &a1, &a2);
        double h = 1e-6;
        CfBodyVector f1 = scale(sub(NozzleDirection(n, d1 + h, d2), NozzleDirection(n, d1 - h, d2)), 0.5 / h);
        CfBodyVector f2 = scale(sub(NozzleDirection(n, d1, d2 + h), NozzleDirection(n, d1, d2 - h)), 0.5 / h);
        ed = fmax(ed, fmax(norm(sub(f1, a1)), norm(sub(f2, a2))));
    }
    Check("direction is unit (2000 angle pairs)", eu < 1e-15);
    Check("cone relation cos = cos(d1) cos(d2)", ec < 1e-15);
    Check("direction derivatives vs central differences", ed < 1e-8);
    double d1 = Deg(5.0);
    double d2 = Deg(5.0);
    bool clipped = ConeScale(&d1, &d2, Deg(6.0));
    Check("cone scale: on the boundary, direction kept",
          clipped && fabs(cos(d1) * cos(d2) - cos(Deg(6.0))) < 1e-15 && fabs(d1 - d2) < 1e-15);
}

static void TestForces()
{
    EngineCommand x;
    memset(&x, 0, sizeof(x));
    for (int i = 0; i < 5; i++)
    {
        x.throttle[i] = 1.0;
    }
    CfBodyVector p = mkBody(3.0, 0.0, 0.0);
    EngineForces f;
    EngineForcesOf(g_cfg, x, 0, 0.0, p, false, &f);
    double F = NozzleThrust(g_cfg, 1.0, 0.0, 1.0);
    Check("forces: zero deflection", fabs(f.force.bu - 5.0 * F) < 1e-6 && norm(momVec(f.moment)) < 1e-6);
    Check("forces: mass flow", fabs(f.mdot - 500.0) < 1e-12);
    double d = Deg(4.0);
    for (int i = 1; i < 5; i++)
    {
        x.angle1[i] = d;
    }
    EngineForcesOf(g_cfg, x, 0, 0.0, p, false, &f);
    bool roll = fabs(f.moment.bu - 4.0 * F * sin(d)) < 1e-6 && fabs(f.force.be) < 1e-6 && fabs(f.force.bn) < 1e-6;
    Check("forces: tangential deflection gives pure roll 4 R F sin d", roll);
    memset(&x.angle1, 0, sizeof(x.angle1));
    x.angle1[0] = d;
    EngineForcesOf(g_cfg, x, 0, 0.0, p, false, &f);
    bool pitch = fabs(f.force.be - F * sin(d)) < 1e-6 && fabs(f.moment.bn + 3.0 * F * sin(d)) < 1e-6;
    Check("forces: centre deflection, moment arm 3 m", pitch);
    unsigned char on[ENGINE_MAX_NOZZLES];
    memset(on, 1, sizeof(on));
    on[2] = 0;
    memset(&x.angle1, 0, sizeof(x.angle1));
    EngineForcesOf(g_cfg, x, on, 0.0, p, false, &f);
    Check("forces: failed nozzle gives no thrust", f.thrust[2] == 0.0 && fabs(f.moment.be + F) < 1e-6);
}

static void TestModel()
{
    static EngineModel m;
    EngineConfig c = g_cfg;
    c.throttleDelay = 0.05;
    c.nozzle[0].delay = 0.03;
    Check("model: init", m.init(&c, 0.01));
    m.reset(0.5);
    EngineCommand cmd = m.actual();
    cmd.throttle[1] = 1.0;
    cmd.angle1[0] = Deg(5.0);
    cmd.angle2[0] = Deg(5.0);
    bool delayOk = true;
    bool rateOk = true;
    for (int k = 1; k <= 80; k++)
    {
        m.step(cmd);
        double u = m.actual().throttle[1];
        double want = k <= 5 ? 0.5 : fmin(1.0, 0.5 + 0.01 * (k - 5));
        delayOk = delayOk && fabs(u - want) < 1e-12;
        double a = m.actual().angle1[0];
        rateOk = rateOk && (k <= 3 ? a == 0.0 : (a > 0.0 && a <= Deg(10.0) * 0.01 * (k - 3) + 1e-15));
    }
    Check("model: throttle delay 5 steps then rate 1/s", delayOk);
    Check("model: servo delay 3 steps then rate limit", rateOk);
    const EngineCommand& a = m.actual();
    Check("model: actual angles held in the cone", ConeInside(a.angle1[0], a.angle2[0], Deg(6.0) + 1e-12));
}

// KKT conditions of a bounded least-squares solution.
static bool Kkt(const double* A, const double* b, int m, int n, const double* lb,
                const double* ub, const double* x)
{
    for (int j = 0; j < n; j++)
    {
        double g = 0.0;
        for (int i = 0; i < m; i++)
        {
            double r = -b[i];
            for (int k = 0; k < n; k++)
            {
                r += A[i * n + k] * x[k];
            }
            g += A[i * n + j] * r;
        }
        bool atL = x[j] <= lb[j];
        bool atU = x[j] >= ub[j];
        bool ok = atL ? g >= -1e-9 : (atU ? g <= 1e-9 : fabs(g) < 1e-9);
        if (!ok || x[j] < lb[j] || x[j] > ub[j])
        {
            return false;
        }
    }
    return true;
}

static void TestBvls()
{
    static Bvls ls;
    bool ok = true;
    for (int s = 0; s < 300; s++)
    {
        int n = 3 + s % 10;
        int m = n + 2 + s % 5;
        double A[BVLS_MAX_ROWS * BVLS_MAX_COLS];
        double b[BVLS_MAX_ROWS];
        double lb[BVLS_MAX_COLS];
        double ub[BVLS_MAX_COLS];
        double x[BVLS_MAX_COLS];
        for (int i = 0; i < m * n; i++)
        {
            A[i] = Uni(-1.0, 1.0);
        }
        for (int i = 0; i < m; i++)
        {
            b[i] = Uni(-3.0, 3.0);
        }
        for (int j = 0; j < n; j++)
        {
            lb[j] = Uni(-1.0, 0.0);
            ub[j] = lb[j] + Uni(0.1, 1.0);
            x[j] = 0.5 * (lb[j] + ub[j]);
        }
        int it = ls.solve(A, b, m, n, lb, ub, x, 200);
        ok = ok && it > 0 && Kkt(A, b, m, n, lb, ub, x);
    }
    Check("bvls: KKT conditions on 300 random problems", ok);
}


// Optimality of a pair: -g in the cone of the normals of the active
// constraints (2-D: one normal or a pair of them with non-negative
// multipliers). g -- gradient of the pair; an, slack -- normals (unit) and
// slacks of its constraints; na -- count
static bool PairKkt(const double g[2], const double an[][2], const double* slack, int na)
{
    double gn = sqrt(g[0] * g[0] + g[1] * g[1]);
    double tol = 1e-8 * (1.0 + gn);
    if (gn < tol)
    {
        return true;
    }
    for (int a = 0; a < na; a++)
    {
        if (slack[a] > 1e-8)
        {
            continue;
        }
        double l = -(g[0] * an[a][0] + g[1] * an[a][1]);
        double r0 = -g[0] - l * an[a][0];
        double r1 = -g[1] - l * an[a][1];
        if (l >= -tol && sqrt(r0 * r0 + r1 * r1) < tol)
        {
            return true;
        }
        for (int c = a + 1; c < na; c++)
        {
            if (slack[c] > 1e-8)
            {
                continue;
            }
            double det = an[a][0] * an[c][1] - an[c][0] * an[a][1];
            if (fabs(det) < 1e-6)
            {
                continue;
            }
            double la = (-g[0] * an[c][1] + g[1] * an[c][0]) / det;
            double lc = (-an[a][0] * g[1] + an[a][1] * g[0]) / det;
            if (la >= -tol / fabs(det) && lc >= -tol / fabs(det))
            {
                return true;
            }
        }
    }
    return false;
}

// Random problems with bounds and linear constraints on pairs (0, 1),
// (2, 3), (4, 5): KKT conditions; one pair of each third problem has a
// 64-gon (tangents of a circle).
static void TestLsc()
{
    static Lsc ls;
    bool ok = true;
    int worstIt = 0;
    for (int s = 0; s < 400; s++)
    {
        int n = 4 + s % 9;
        int m = n + 2 + s % 5;
        double A[LSC_MAX_ROWS * LSC_MAX_COLS];
        double b[LSC_MAX_ROWS];
        double lb[LSC_MAX_COLS];
        double ub[LSC_MAX_COLS];
        double x[LSC_MAX_COLS];
        for (int i = 0; i < m * n; i++)
        {
            A[i] = Uni(-1.0, 1.0);
        }
        for (int i = 0; i < m; i++)
        {
            b[i] = Uni(-3.0, 3.0);
        }
        for (int j = 0; j < n; j++)
        {
            lb[j] = Uni(-1.0, -0.2);
            ub[j] = Uni(0.2, 1.0);
            x[j] = 0.0;
        }
        ls.setup(n, lb, ub);
        int np = n >= 8 ? 3 : 2;
        static double cn[3][80][2];
        static double ch[3][80];
        int nc[3] = { 0, 0, 0 };
        for (int p = 0; p < np; p++)
        {
            int id = ls.pair(2 * p, 2 * p + 1);
            bool gon = p == 0 && s % 3 == 0;
            int k = gon ? 64 : 2 + s % 9;
            double r = Uni(0.15, 0.8);
            for (int c = 0; c < k; c++)
            {
                double th = gon ? c * (2.0 * PI / 64.0) : Uni(0.0, 2.0 * PI);
                double h = gon ? r : Uni(0.15, 0.8);
                cn[p][c][0] = cos(th);
                cn[p][c][1] = sin(th);
                ch[p][c] = h;
                ls.cut(id, cos(th), sin(th), h);
            }
            nc[p] = k;
        }
        int it = ls.solve(A, b, m, x, 400);
        worstIt = it > worstIt ? it : worstIt;
        ok = ok && it > 0;
        double r[LSC_MAX_ROWS];
        for (int i = 0; i < m; i++)
        {
            r[i] = -b[i];
            for (int k = 0; k < n; k++)
            {
                r[i] += A[i * n + k] * x[k];
            }
        }
        double g[LSC_MAX_COLS];
        for (int j = 0; j < n; j++)
        {
            g[j] = 0.0;
            for (int i = 0; i < m; i++)
            {
                g[j] += A[i * n + j] * r[i];
            }
            ok = ok && x[j] >= lb[j] - 1e-12 && x[j] <= ub[j] + 1e-12;
        }
        for (int j = 2 * np; j < n; j++)
        {
            bool atL = x[j] <= lb[j];
            bool atU = x[j] >= ub[j];
            ok = ok && (atL ? g[j] >= -1e-9 : (atU ? g[j] <= 1e-9 : fabs(g[j]) < 1e-9));
        }
        for (int p = 0; p < np; p++)
        {
            double an[84][2];
            double sl[84];
            int na = 0;
            for (int v = 0; v < 2; v++)
            {
                int j = 2 * p + v;
                an[na][0] = v == 0 ? -1.0 : 0.0;
                an[na][1] = v == 0 ? 0.0 : -1.0;
                sl[na++] = x[j] - lb[j];
                an[na][0] = v == 0 ? 1.0 : 0.0;
                an[na][1] = v == 0 ? 0.0 : 1.0;
                sl[na++] = ub[j] - x[j];
            }
            for (int c = 0; c < nc[p]; c++)
            {
                an[na][0] = cn[p][c][0];
                an[na][1] = cn[p][c][1];
                sl[na] = ch[p][c] - cn[p][c][0] * x[2 * p] - cn[p][c][1] * x[2 * p + 1];
                ok = ok && sl[na] >= -1e-12;
                na++;
            }
            double gp[2] = { g[2 * p], g[2 * p + 1] };
            ok = ok && PairKkt(gp, an, sl, na);
        }
    }
    printf("  lsc: largest iteration count %d\n", worstIt);
    Check("lsc: KKT conditions on 400 random problems with pair constraints", ok);
}

int main()
{
    printf("gnclib test_engine\n");
    MakeConfig(&g_cfg);
    Check("config check", EngineConfigCheck(g_cfg));
    TestThrustAndKinematics();
    TestForces();
    TestModel();
    TestBvls();
    TestLsc();
    TestAllocator();
    printf("%s (%d failed)\n", g_bad == 0 ? "TEST_ENGINE PASS" : "TEST_ENGINE FAIL", g_bad);
    return g_bad == 0 ? 0 : 1;
}

// Random admissible command (within limits and the cone).
static void RandomCommand(EngineCommand* x, const unsigned char* on)
{
    memset(x, 0, sizeof(*x));
    for (int i = 0; i < g_cfg.count; i++)
    {
        const NozzleConfig& n = g_cfg.nozzle[i];
        x->throttle[i] = Uni(g_cfg.throttleMin, g_cfg.throttleMax);
        x->angle1[i] = n.kind != NOZZLE_FIXED ? Uni(-n.limit1, n.limit1) : 0.0;
        x->angle2[i] = n.kind == NOZZLE_FREE ? Uni(-n.limit2, n.limit2) : 0.0;
        NozzleClampAngles(n, &x->angle1[i], &x->angle2[i]);
    }
    (void)on;
}

// Allocation of requests produced by admissible commands: always feasible.
static int AllocFeasible(EngineAllocator& al, const unsigned char* on, double pa, int count,
                         double* worst)
{
    int bad = 0;
    *worst = 0.0;
    for (int s = 0; s < count; s++)
    {
        EngineCommand x;
        RandomCommand(&x, on);
        AllocRequest q;
        memset(&q, 0, sizeof(q));
        memcpy(q.on, on, sizeof(q.on));
        q.point = mkBody(3.0, 0.0, 0.0);
        q.pa = pa;
        EngineForces f;
        EngineForcesOf(g_cfg, x, on, pa, q.point, false, &f);
        q.thrust = f.force.bu;
        q.moment = f.moment;
        EngineCommand cmd;
        AllocResult r;
        bool ok = al.allocate(q, &cmd, &r);
        double e = fmax(r.errThrust, fmax(r.errRoll, r.errPitch));
        *worst = fmax(*worst, e);
        bool inCone = ConeInside(cmd.angle1[0], cmd.angle2[0], Deg(6.0) + 1e-12);
        bad += (ok && r.feasible && inCone) ? 0 : 1;
    }
    return bad;
}

static void TestAllocator()
{
    static EngineAllocator al;
    AllocConfig ac;
    AllocConfigDefaults(&ac);
    ac.rateWindow = 0;
    Check("alloc: init", al.init(&g_cfg, &ac));
    unsigned char on[ENGINE_MAX_NOZZLES];
    memset(on, 1, sizeof(on));
    double worst = 0.0;
    int bad = AllocFeasible(al, on, 0.0, 2000, &worst);
    printf("  worst normalized residual %.3e\n", worst);
    Check("alloc: 2000 feasible requests, vacuum", bad == 0);
    bad = AllocFeasible(al, on, 101325.0, 1000, &worst);
    Check("alloc: 1000 feasible requests, sea level", bad == 0);
    on[3] = 0;
    bad = AllocFeasible(al, on, 0.0, 1000, &worst);
    Check("alloc: 1000 feasible requests, nozzle 3 failed", bad == 0);
    on[3] = 1;
    AllocRequest q;
    memset(&q, 0, sizeof(q));
    memset(q.on, 1, sizeof(q.on));
    q.point = mkBody(3.0, 0.0, 0.0);
    q.thrust = 2.0 * 5.0 * NozzleThrust(g_cfg, 1.0, 0.0, 1.0);
    q.moment = mkMom(1e4, 2e4, 0.0);
    EngineCommand cmd;
    AllocResult r;
    al.allocate(q, &cmd, &r);
    Check("alloc: thrust over maximum -- moments kept, thrust short",
          r.errPitch < 1e-6 && r.errRoll < 1e-6 && r.errThrust > 0.4);
    q.thrust = 3.0 * NozzleThrust(g_cfg, 1.0, 0.0, 1.0);
    q.moment = mkMom(1e7, 2e4, 0.0);
    al.allocate(q, &cmd, &r);
    Check("alloc: roll over capacity -- pitch kept, roll short", r.errPitch < 1e-6 && r.errRoll > 1e-3);
    ac.rateWindow = 1;
    ac.dt = 0.01;
    al.init(&g_cfg, &ac);
    static EngineModel m;
    m.init(&g_cfg, 0.01);
    m.reset(0.7);
    al.reset(m.actual());
    bool rateOk = true;
    bool follow = true;
    for (int k = 0; k < 500; k++)
    {
        q.thrust = (3.5 + 0.5 * sin(0.02 * k)) * NozzleThrust(g_cfg, 1.0, 0.0, 1.0);
        q.moment = mkMom(2e4 * sin(0.03 * k), 1.5e5 * sin(0.05 * k), 1e5 * cos(0.04 * k));
        EngineCommand prev = al.last();
        al.allocate(q, &cmd, &r);
        for (int i = 0; i < 5; i++)
        {
            rateOk = rateOk && fabs(cmd.throttle[i] - prev.throttle[i]) <= 0.01 + 1e-15;
            rateOk = rateOk && fabs(cmd.angle1[i] - prev.angle1[i]) <= Deg(10.0) * 0.01 + 1e-15;
            rateOk = rateOk && fabs(cmd.angle2[i] - prev.angle2[i]) <= Deg(10.0) * 0.01 + 1e-15;
        }
        m.step(cmd);
        for (int i = 0; i < 5; i++)
        {
            const EngineCommand& a = m.actual();
            follow = follow && a.throttle[i] == cmd.throttle[i];
            follow = follow && a.angle1[i] == cmd.angle1[i] && a.angle2[i] == cmd.angle2[i];
        }
    }
    Check("alloc: rate window -- command steps within rate * dt", rateOk);
    Check("alloc: plant without delay reproduces commands bit-exactly", follow);
}

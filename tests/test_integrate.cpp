// test_integrate.cpp -- gnclib 1.1.0: IntegrateRK4 (stage forces, aux
// states) and the rotating Earth (lib_earth.h); 1.7.1: attitude order (f);
// 1.11.0: Earth rate in the rotational dynamics (g).
// Build: g++ -std=c++17 -O2 -Iinclude tests/test_integrate.cpp -o build/test_integrate
#include "lib.h"
#include "lib_integrate.h"
#include "lib_earth.h"
#include <stdio.h>
#include <math.h>

using namespace gnc;

static int g_bad = 0;
static void Check(const char* name, double val, double lim)
{
    int ok = (fabs(val) <= lim) ? 1 : 0;
    printf("%-54s %12.3e  (limit %.1e)  %s\n", name, val, lim, ok ? "OK" : "FAIL");
    g_bad += ok ? 0 : 1;
}

static INR Inr(double m, double a, double b, double c)
{
    INR r { };
    r.m = m;
    r.j_bubu = a;
    r.j_bebe = b;
    r.j_bnbn = c;
    return r;
}

// (a) constant force and moment
struct CtxConst { CfLocalVector f; CfBodyMoment m; INR inr; };
static void FnConst(void* c, const VehState&, const AuxState&, DynDeriv* d)
{
    CtxConst* k = (CtxConst*)c;
    d->f = k->f;
    d->mom = k->m;
    d->inr = k->inr;
}
// (b) torque-free body
static void FnFree(void* c, const VehState&, const AuxState&, DynDeriv* d)
{
    CtxConst* k = (CtxConst*)c;
    CfLocalVector z { };
    CfBodyMoment zm { };
    d->f = z;
    d->mom = zm;
    d->inr = k->inr;
}
// (c) aux harmonic oscillator x'' = -w^2 x
static void FnOsc(void*, const VehState&, const AuxState& a, DynDeriv* d)
{
    const double w = 5.0;
    d->inr = Inr(1.0, 1.0, 1.0, 1.0);
    d->xdot[0] = a.x[1];
    d->xdot[1] = -w * w * a.x[0];
}
// (d) state-dependent force: spring on the position (order check)
static void FnSpring(void*, const VehState& s, const AuxState&, DynDeriv* d)
{
    d->inr = Inr(1.0, 1.0, 1.0, 1.0);
    d->f.e = -4.0 * s.PosLop().e;   // w = 2 rad/s
}
// (e) free fall in the rotating Earth
struct CtxEarth { const GEO* geo; };
static void FnEarth(void* c, const VehState& s, const AuxState&, DynDeriv* d)
{
    CtxEarth* k = (CtxEarth*)c;
    d->inr = Inr(1.0, 1.0, 1.0, 1.0);
    CfLocalVector g = EarthGravityLoc(*k->geo, s.PosLop());
    CfLocalVector cor = CoriolisAccLoc(EarthRateLoc(*k->geo), s.VelLoc());
    d->f.e = g.e + cor.e;
    d->f.n = g.n + cor.n;
    d->f.u = g.u + cor.u;
}

static VehState Init(double pe, double vn, const CfBodyOmega& w)
{
    VehState s;
    CfLocalPoint p { };
    p.e = pe;
    CfLocalVector v { };
    v.n = vn;
    CfLocalVector a { };
    s.SetLinear(p, v, a);
    CfRotation att = s.Att();
    CfBodyOmegaAcc e { };
    s.SetAngular(att, w, e);
    return s;
}

// (f) torque-free body: angle between two attitudes (vector part of
// conj(a) * b, accurate for small angles), rad.
static double AngleBetween(const CfRotation& a, const CfRotation& b)
{
    double cx = a.Qe() * b.Qu() - a.Qu() * b.Qe();
    double cy = a.Qu() * b.Qn() - a.Qn() * b.Qu();
    double cz = a.Qn() * b.Qe() - a.Qe() * b.Qn();
    double vx = a.Qw() * b.Qn() - b.Qw() * a.Qn() - cx;
    double vy = a.Qw() * b.Qe() - b.Qw() * a.Qe() - cy;
    double vz = a.Qw() * b.Qu() - b.Qw() * a.Qu() - cz;
    return 2.0 * sqrt(vx * vx + vy * vy + vz * vz);
}

// (f) angular momentum in the LOCAL axes (depends on the attitude).
static CfLocalVector HLocal(const VehState& s, const INR& j)
{
    const CfBodyOmega& w = s.OmgBdw();
    CfBodyVector hb { };
    hb.bn = j.j_bnbn * w.bn;
    hb.be = j.j_bebe * w.be;
    hb.bu = j.j_bubu * w.bu;
    return s.Att().B2L(hb);
}

// (f) torque-free run of T seconds with the step dt from the rate w0.
static VehState RunFree(CtxConst* k, const CfBodyOmega& w0, double dt, double T)
{
    VehState s = Init(0.0, 0.0, w0);
    int n = (int)(T / dt + 0.5);
    for (int i = 0; i < n; ++i)
    {
        s = IntegrateRK4(s, 0, dt, FnFree, k);
    }
    return s;
}

// (g) 1.11.0: torque-free body on the rotating Earth (DynDeriv::earthRate)
struct CtxRot { INR inr; CfLocalOmega wE; };
static void FnRot(void* c, const VehState&, const AuxState&, DynDeriv* d)
{
    CtxRot* k = (CtxRot*)c;
    d->inr = k->inr;
    d->earthRate = k->wE;
}

// (g) local vector v turned about the unit axis k by ang (Rodrigues).
static CfLocalVector Turn(const CfLocalVector& v, const CfLocalVector& k, double ang)
{
    double c = cos(ang);
    double s = sin(ang);
    double kv = k.e * v.e + k.n * v.n + k.u * v.u;
    CfLocalVector x { };
    x.e = k.n * v.u - k.u * v.n;
    x.n = k.u * v.e - k.e * v.u;
    x.u = k.e * v.n - k.n * v.e;
    CfLocalVector r { };
    r.e = v.e * c + x.e * s + k.e * kv * (1.0 - c);
    r.n = v.n * c + x.n * s + k.n * kv * (1.0 - c);
    r.u = v.u * c + x.u * s + k.u * kv * (1.0 - c);
    return r;
}

// (g) nose (bu) of a torque-free body after T s; the body rate is relative
// to the Earth, w0 -- its initial value.
static CfLocalVector NoseAfter(CtxRot* k, const CfBodyOmega& w0, double dt, double T)
{
    VehState s = Init(0.0, 0.0, w0);
    int n = (int)(T / dt + 0.5);
    for (int i = 0; i < n; ++i)
    {
        s = IntegrateRK4(s, 0, dt, FnRot, k);
    }
    CfBodyVector bu { };
    bu.bu = 1.0;
    return s.Att().B2L(bu);
}

int main()
{
    // (a) constant F, M: IntegrateRK4 == Integrate (translation exact, rotation RK4)
    {
        CtxConst k { };
        k.f.e = 3.0; k.f.n = -1.0; k.f.u = 2.0;
        k.m.bu = 0.2; k.m.be = -0.1; k.m.bn = 0.05;
        k.inr = Inr(10.0, 1.0, 2.0, 3.0);
        CfBodyOmega w0 { };
        w0.bu = 0.3; w0.be = 0.1; w0.bn = -0.2;
        VehState a = Init(0.0, 1.0, w0);
        VehState b = a;
        for (int i = 0; i < 2000; ++i)
        {
            a = IntegrateRK4(a, 0, 0.001, FnConst, &k);
            b = Integrate(b, 0.001, k.f, k.m, k.inr);
        }
        double dp = fabs(a.PosLop().e - b.PosLop().e) + fabs(a.PosLop().n - b.PosLop().n)
                    + fabs(a.PosLop().u - b.PosLop().u);
        double dw = fabs(a.OmgBdw().bu - b.OmgBdw().bu) + fabs(a.OmgBdw().be - b.OmgBdw().be)
                    + fabs(a.OmgBdw().bn - b.OmgBdw().bn);
        Check("(a) const F,M: |p_RK4 - p_Integrate| after 2 s, m", dp, 1e-9);
        Check("(a) const F,M: |w_RK4 - w_Integrate| after 2 s, rad/s", dw, 1e-12);
    }
    // (b) torque-free asymmetric body: energy and |H| conserved
    {
        CtxConst k { };
        k.inr = Inr(1.0, 1.0, 2.0, 3.0);
        CfBodyOmega w0 { };
        w0.bu = 0.01; w0.be = 2.0; w0.bn = 0.01;   // near the intermediate axis
        VehState s = Init(0.0, 0.0, w0);
        double E0 = 0.5 * (1.0 * w0.bu * w0.bu + 2.0 * w0.be * w0.be + 3.0 * w0.bn * w0.bn);
        double H0 = sqrt(1.0 * w0.bu * w0.bu * 1.0 + 4.0 * w0.be * w0.be + 9.0 * w0.bn * w0.bn);
        for (int i = 0; i < 10000; ++i)
        {
            s = IntegrateRK4(s, 0, 0.001, FnFree, &k);
        }
        CfBodyOmega w = s.OmgBdw();
        double E = 0.5 * (1.0 * w.bu * w.bu + 2.0 * w.be * w.be + 3.0 * w.bn * w.bn);
        double H = sqrt(w.bu * w.bu + 4.0 * w.be * w.be + 9.0 * w.bn * w.bn);
        Check("(b) torque-free 10 s: relative energy drift", (E - E0) / E0, 1e-9);
        Check("(b) torque-free 10 s: relative |H| drift", (H - H0) / H0, 1e-9);
    }
    // (c) auxiliary states: x'' = -25 x, x(0) = 1
    {
        AuxState aux { };
        aux.n = 2;
        aux.x[0] = 1.0;
        aux.x[1] = 0.0;
        CfBodyOmega w0 { };
        VehState s = Init(0.0, 0.0, w0);
        for (int i = 0; i < 10000; ++i)
        {
            s = IntegrateRK4(s, &aux, 0.001, FnOsc, 0);
        }
        Check("(c) aux oscillator 10 s: |x - cos(5t)|", aux.x[0] - cos(50.0), 1e-9);
    }
    // (d) state-dependent force: 4th order
    {
        double err[2];
        double dts[2] = { 0.02, 0.01 };
        for (int j = 0; j < 2; ++j)
        {
            CfBodyOmega w0 { };
            VehState s = Init(1.0, 0.0, w0);
            int n = (int)(5.0 / dts[j] + 0.5);
            for (int i = 0; i < n; ++i)
            {
                s = IntegrateRK4(s, 0, dts[j], FnSpring, 0);
            }
            err[j] = fabs(s.PosLop().e - cos(10.0));
        }
        double order = log(err[0] / err[1]) / log(2.0);
        Check("(d) spring on position: convergence order - 4", order - 4.0, 0.15);
    }
    // (e) free fall from 100 m at 45 N: east deflection (1/3) w g t^3 cos(lat)
    {
        GEO geo;
        geo.SetOrigin(mkGdt(45.0 * PI / 180.0, 30.0 * PI / 180.0, 100.0));
        CtxEarth k { &geo };
        CfBodyOmega w0 { };
        VehState s = Init(0.0, 0.0, w0);
        CfLocalPoint p { };
        p.u = 100.0;
        CfLocalVector z { };
        s.SetLinear(p, z, z);
        double t = 0.0;
        while (s.PosLop().u > 0.0)
        {
            s = IntegrateRK4(s, 0, 0.001, FnEarth, &k);
            t += 0.001;
        }
        double g = GravityNormalWgs84(45.0 * PI / 180.0, 100.0);
        double tf = sqrt(2.0 * 100.0 / g);
        double east = OMEGA_E * g * tf * tf * tf * cos(45.0 * PI / 180.0) / 3.0;
        printf("    free fall %.3f s: east %.3f mm (analytic %.3f mm)\n", t, s.PosLop().e * 1e3,
               east * 1e3);
        Check("(e) Coriolis free fall: relative east deflection error", s.PosLop().e / east - 1.0,
              0.02);
        CfLocalVector gl = EarthGravityLoc(geo, mkLop(0.0, 0.0, 0.0));
        Check("(e) |g| at the site vs GravityNormalWgs84, m/s^2",
              sqrt(gl.e * gl.e + gl.n * gl.n + gl.u * gl.u) - g, 1e-3);
    }
    // (f) torque-free asymmetric body tumbling at about 2.3 rad/s: the
    // attitude converges with the 4th order (RKMK4, 1.7.1; up to 1.7.0 --
    // 2nd), the angular momentum in the LOCAL axes is conserved
    {
        CtxConst k { };
        k.inr = Inr(1.0, 1.0, 2.0, 3.0);
        CfBodyOmega w0 { };
        w0.bu = 0.7;
        w0.be = 2.0;
        w0.bn = -1.0;
        VehState ref = RunFree(&k, w0, 1.0e-4, 10.0);
        VehState s1 = RunFree(&k, w0, 0.02, 10.0);
        VehState s2 = RunFree(&k, w0, 0.01, 10.0);
        double e1 = AngleBetween(s1.Att(), ref.Att());
        double e2 = AngleBetween(s2.Att(), ref.Att());
        printf("    attitude error after 10 s: dt 20 ms %.3e rad, dt 10 ms %.3e rad\n", e1, e2);
        Check("(f) torque-free: attitude convergence order - 4", log(e1 / e2) / log(2.0) - 4.0, 0.3);
        VehState s0 = Init(0.0, 0.0, w0);
        VehState s3 = RunFree(&k, w0, 0.002, 20.0);
        CfLocalVector h0 = HLocal(s0, k.inr);
        CfLocalVector h3 = HLocal(s3, k.inr);
        Check("(f) dt 2 ms, 20 s: rel. drift of H, local axes",
              Norm(Sub(h3, h0)) / Norm(h0), 1.0e-11);
    }
    // (g) 1.11.0: torque-free bodies on the rotating Earth (Plesetsk). The
    // body rate is relative to the Earth; a body with zero ABSOLUTE rate
    // (relative rate -om) keeps its direction in the inertial space, so in
    // the local axes its nose turns at -om about the Earth axis; a fast
    // spinner keeps its axis in the inertial space (gyroscopic terms of the
    // absolute rate); a body at rest relative to the Earth co-rotates
    {
        GEO geo;
        CfGeodetic site { };
        site.lat = 62.928008 * PI / 180.0;
        site.lon = 40.456137 * PI / 180.0;
        geo.SetOrigin(site);
        CtxRot k { };
        k.inr = Inr(1.0, 2.0, 5.0, 5.0);
        k.wE = EarthRateOmg(geo);
        CtxRot k0 = k;
        CfLocalOmega z { };
        k0.wE = z;
        CfLocalVector ax = EarthRateLoc(geo);
        double om = Norm(ax);
        ax.e = ax.e / om;
        ax.n = ax.n / om;
        ax.u = ax.u / om;
        CfLocalVector up { };
        up.u = 1.0;
        // zero absolute rate (identity attitude: bu = u, be = e, bn = n), 158 s
        CfBodyOmega wi { };
        wi.bu = -k.wE.u;
        wi.be = -k.wE.e;
        wi.bn = -k.wE.n;
        CfLocalVector want = Turn(up, ax, -om * 158.0);
        Check("(g) zero absolute rate, 158 s: nose vs inertial, rad",
              Norm(Sub(NoseAfter(&k, wi, 0.01, 158.0), want)), 1.0e-9);
        // spin 10 rad/s (absolute) about the symmetry axis bu, 60 s
        CfBodyOmega ws = wi;
        ws.bu = 10.0 - k.wE.u;
        CfLocalVector want60 = Turn(up, ax, -om * 60.0);
        double eNew = Norm(Sub(NoseAfter(&k, ws, 0.001, 60.0), want60));
        double eOld = Norm(Sub(NoseAfter(&k0, ws, 0.001, 60.0), want60));
        printf("    spinner, 60 s: nose off inertial %.3e rad (earthRate 0: %.3e rad)\n",
               eNew, eOld);
        Check("(g) spinning 10 rad/s, 60 s: nose vs inertial, rad", eNew, 1.0e-8);
        // at rest relative to the Earth, spherical body: co-rotates, stays put
        // in the local axes (an unsymmetric one nutates: rotation with the
        // Earth about a non-principal axis is not steady)
        CtxRot ks = k;
        ks.inr = Inr(1.0, 5.0, 5.0, 5.0);
        CfBodyOmega w0 { };
        Check("(g) spherical body at rest on the Earth, 158 s: nose vs up, rad",
              Norm(Sub(NoseAfter(&ks, w0, 0.01, 158.0), up)), 1.0e-12);
    }
    printf(g_bad == 0 ? "TEST_INTEGRATE PASS\n" : "TEST_INTEGRATE FAIL\n");
    return g_bad;
}

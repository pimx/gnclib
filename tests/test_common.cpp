// test_common.cpp -- gnclib 1.6.0: limits and angles, tables, observers,
// notch filter, gusts, wind profile, orbit helpers. Where a reference
// implementation is written out below, the library function must give
// its numbers to the bit.
// Build: g++ -std=c++17 -O2 -Iinclude tests/test_common.cpp build/libgnc.a
//        -o build/test_common
#include "lib.h"
#include "lib_eso.h"
#include "lib_gust.h"
#include "lib_notch.h"
#include "lib_orbit.h"
#include "lib_random.h"
#include "lib_table.h"
#include "lib_util.h"
#include "lib_wind_profile.h"
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

// Two doubles equal to the bit.
static bool Same(double a, double b)
{
    return memcmp(&a, &b, sizeof(double)) == 0;
}

// ---- references ---------------------------------------------------------
static double RefLimit(double v, double limit)
{
    return v > limit ? limit : (v < -limit ? -limit : v);
}

static double RefClamp(double v, double lo, double hi)
{
    if (v < lo)
    {
        return lo;
    }
    if (v > hi)
    {
        return hi;
    }
    return v;
}

static constexpr double RefDeg(double d)
{
    return d * (PI / 180.0);
}

static const int RN = 4;
static const int RM = 6;

static double RefBilinear(const double tab[RN][RM], int iD, double fD, int iM, double fM)
{
    double a = tab[iD][iM] + fM * (tab[iD][iM + 1] - tab[iD][iM]);
    double b = tab[iD + 1][iM] + fM * (tab[iD + 1][iM + 1] - tab[iD + 1][iM]);
    return a + fD * (b - a);
}

static double RefPhase(int seed, int component, int tone)
{
    unsigned int x = (unsigned int)(seed * 2654435761u) ^ (unsigned int)(component * 40503u)
                     ^ (unsigned int)(tone * 2246822519u);
    x = x ^ (x >> 13);
    x = x * 1274126177u;
    x = x ^ (x >> 16);
    return 2.0 * 3.14159265358979323846 * ((double)(x % 100000u) / 100000.0);
}

static double RefComponent(double sigma, double tau, int seed, double time, int component)
{
    (void)sigma;
    double sum = 0.0;
    double norm = 0.0;
    for (int k = 0; k < 8; ++k)
    {
        double period = tau * (1.0 + 19.0 * (double)k / (double)(8 - 1));
        double amp = sqrt(period / tau);
        double w = 2.0 * 3.14159265358979323846 / period;
        sum = sum + amp * sin(w * time + RefPhase(seed, component, k));
        norm = norm + 0.5 * amp * amp;
    }
    return norm > 1.0e-12 ? sum / sqrt(norm) : 0.0;
}

static double RefAzimuth(double lat, double target, double orbitAlt, bool south)
{
    const double PI_ = 3.14159265358979323846;
    const double MU = ELL_PZ9011.mu;
    const double RE = ELL_PZ9011.a;
    double sinA = cos(target) / cos(lat);
    double a0 = 0.5 * PI_;
    if (sinA < 1.0)
    {
        double root = asin(sinA);
        double ai = south ? PI_ - root : root;
        double v = sqrt(MU / (RE + orbitAlt));
        double vPad = ELL_PZ9011.om * RE * cos(lat);
        a0 = atan2(v * sin(ai) - vPad, v * cos(ai));
    }
    return a0;
}

struct RefEso2
{
    double wo = 1.0;
    double z1 = 0.0;
    double z2 = 0.0;
    void update(double vMeas, double uApplied, double dt)
    {
        double q = 1.0 - exp(-wo * dt);
        double e = vMeas - z1;
        z1 += dt * (uApplied + z2) + 2.0 * q * e;
        z2 += q * q / dt * e;
    }
};

// ---- groups -------------------------------------------------------------
static void TestUtil(RNG* rng)
{
    bool ok = true;
    for (int i = 0; i < 200000; ++i)
    {
        double v = (0.0 + 10.0 * rng->Gauss());
        double lim = fabs((0.0 + 5.0 * rng->Gauss()));
        double lo = (0.0 + 5.0 * rng->Gauss());
        double hi = lo + fabs((0.0 + 5.0 * rng->Gauss()));
        ok = ok && Same(Limit(v, lim), RefLimit(v, lim));
        ok = ok && Same(Clamp(v, lo, hi), RefClamp(v, lo, hi));
        ok = ok && Same(Deg(v * 30.0), RefDeg(v * 30.0));
    }
    Check("Limit, Clamp, Deg == references (200 000 values)", ok);
    Check("ToDeg(Deg(90)) == 90 within 1e-13", fabs(ToDeg(Deg(90.0)) - 90.0) < 1.0e-13);
}

static void TestTable(RNG* rng)
{
    double rows[RN] = { 0.0, 2.0, 5.0, 9.0 };
    double cols[RM] = { 0.3, 0.8, 1.1, 1.5, 3.0, 6.0 };
    double tab[RN][RM];
    for (int r = 0; r < RN; ++r)
    {
        for (int c = 0; c < RM; ++c)
        {
            tab[r][c] = (1.0 + 1.0 * rng->Gauss());
        }
    }
    bool ok = true;
    for (int i = 0; i < 100000; ++i)
    {
        double x = (4.0 + 6.0 * rng->Gauss());
        double y = (2.0 + 3.0 * rng->Gauss());
        double fR = 0.0;
        double fC = 0.0;
        int iR = TableLocate(rows, RN, x, &fR);
        int iC = TableLocate(cols, RM, y, &fC);
        ok = ok && iR >= 0 && iR <= RN - 2 && iC >= 0 && iC <= RM - 2;
        double got = TableBilinear(&tab[0][0], RM, iR, fR, iC, fC);
        ok = ok && Same(got, RefBilinear(tab, iR, fR, iC, fC));
    }
    Check("TableBilinear == reference (100 000 points, edges held)", ok);
    double f = 0.0;
    int i = TableLocate(rows, RN, 3.5, &f);
    Check("TableLocate(3.5) -> row 1, fraction 0.5", i == 1 && Same(f, 0.5));
}

static void TestEso(RNG* rng)
{
    Eso2 a;
    RefEso2 b;
    a.setup(12.0);
    b.wo = 12.0;
    a.reset(0.3);
    b.z1 = 0.3;
    bool ok = true;
    for (int i = 0; i < 20000; ++i)
    {
        double v = sin(0.01 * i) + (0.0 + 0.01 * rng->Gauss());
        double u = cos(0.013 * i);
        a.update(v, u, 0.01);
        b.update(v, u, 0.01);
        ok = ok && Same(a.est(), b.z1) && Same(a.dist(), b.z2) && Same(a.accel(), u + b.z2);
    }
    Check("Eso2 == reference (20 000 steps)", ok);
    // Eso3 on a plant theta'' = u + f with f = 0.7: the estimate converges.
    Eso3 e;
    Check("Eso3 refuses a pole <= 0", !e.Setup(0.0));
    e.Setup(10.0);
    e.Reset(0.0);
    double th = 0.0;
    double w = 0.0;
    for (int i = 0; i < 3000; ++i)
    {
        double u = -2.0 * th - 1.5 * w;
        w = w + 0.001 * (u + 0.7);
        th = th + 0.001 * w;
        e.Update(th, u, 0.001);
    }
    Check("Eso3 disturbance -> 0.7 within 1e-3", fabs(e.Disturbance() - 0.7) < 1.0e-3);
}

static void TestNotch()
{
    Notch n;
    NotchInit(&n);
    double peak = 0.0;
    for (int i = 0; i < 40000; ++i)
    {
        double y = NotchStep(&n, sin(4.0 * 0.0025 * i), 4.0, 0.02, 0.3, 0.0025);
        if (i > 20000)
        {
            peak = fabs(y) > peak ? fabs(y) : peak;
        }
    }
    Check("Notch: a tone at w0 attenuated below zz / zp (0.067)", peak < 0.07);
    Check("Notch: pass-through at w0 dt >= 3", NotchStep(&n, 1.25, 4.0, 0.02, 0.3, 1.0) == 1.25);
}

static void TestGust()
{
    GustTones g;
    g.sigma = 3.0;
    g.tau = 5.0;
    g.seed = 4321;
    g.vertical = 0.5;
    bool ok = true;
    double s2 = 0.0;
    int n = 0;
    for (double t = 0.0; t < 2000.0; t += 0.37)
    {
        CfLocalVector v = GustTonesAt(g, t);
        ok = ok && Same(v.n, 3.0 * RefComponent(3.0, 5.0, g.seed, t, 0));
        ok = ok && Same(v.e, 3.0 * RefComponent(3.0, 5.0, g.seed, t, 1));
        ok = ok && Same(v.u, 3.0 * 0.5 * RefComponent(3.0, 5.0, g.seed, t, 2));
        s2 = s2 + v.n * v.n;
        n = n + 1;
    }
    Check("GustTonesAt == reference (5 406 moments)", ok);
    double rms = sqrt(s2 / n);
    Check("GustTonesAt: rms of n near sigma (2.5 .. 3.5 m/s)", rms > 2.5 && rms < 3.5);
    GustTones off;
    CfLocalVector z = GustTonesAt(off, 10.0);
    Check("GustTonesAt: sigma 0 -- no gust", z.n == 0.0 && z.e == 0.0 && z.u == 0.0);
}

static void TestWindProfile()
{
    const char* path = "build/test_wind_profile.txt";
    FILE* f = fopen(path, "w");
    if (f == 0)
    {
        Check("wind profile: temporary file", false);
        return;
    }
    fprintf(f, "; profile\naltitude,speed,from\n0, 3, 90\n2000 8 120\n\n10000,25,200\n");
    fclose(f);
    WindProfile p;
    char err[256];
    bool ok = ReadWindProfile(path, &p, err, (int)sizeof(err));
    Check("ReadWindProfile: comments, header, commas and spaces", ok && p.count == 3);
    Check("ReadWindProfile: azimuth in rad by Deg", Same(p.row[1].dir, RefDeg(120.0)));
    double sp = 0.0;
    double fr = 0.0;
    WindProfileAt(p, 1000.0, &sp, &fr);
    Check("WindProfileAt: linear between rows", Same(sp, 3.0 + 0.5 * (8.0 - 3.0)));
    WindProfileAt(p, 50000.0, &sp, &fr);
    Check("WindProfileAt: held above the last row", sp == 25.0 && Same(fr, RefDeg(200.0)));
    f = fopen(path, "w");
    fprintf(f, "0 3 90\n2000 8 120\n1500 9 130\n");
    fclose(f);
    ok = ReadWindProfile(path, &p, err, (int)sizeof(err));
    Check("ReadWindProfile: refuses falling altitudes", !ok && strstr(err, "do not rise") != 0);
    remove(path);
}

static void TestOrbit(RNG* rng)
{
    bool ok = true;
    for (int i = 0; i < 20000; ++i)
    {
        double lat = (0.0 + 0.8 * rng->Gauss());
        double inc = fabs((1.0 + 0.5 * rng->Gauss()));
        double alt = 200000.0 + fabs((0.0 + 1.0e5 * rng->Gauss()));
        bool south = (i % 2) == 1;
        double az = LaunchAzimuthStart(lat, inc, alt, south, ELL_PZ9011);
        ok = ok && Same(az, RefAzimuth(lat, inc, alt, south));
    }
    Check("LaunchAzimuthStart == reference (20 000)", ok);
    GEO geo;
    geo.SetEllipsoid(ELL_PZ9011);
    CfGeodetic org { };
    org.lat = Deg(45.0);
    org.lon = Deg(30.0);
    org.alt = 100.0;
    geo.SetOrigin(org);
    CfGeodetic g = GeoAtLevel(geo, 0.0, 1000.0, 100.0);
    Check("GeoAtLevel: height within 1e-9 m", fabs(g.alt - 100.0) < 1.0e-9);
    CfLocalPoint p { };
    CfLocalVector still { };
    InertialState is = InertialOf(geo, p, still);
    double vPad = ELL_PZ9011.om * is.radius * cos(asin(geo.L2E(p).z / is.radius));
    Check("InertialOf: a pad at rest moves with the Earth (1e-9 rel)",
          fabs(is.speed - vPad) < 1.0e-9 * vPad);
    double geocLat = asin(geo.L2E(p).z / is.radius);
    Check("InertialOf: pad at rest -- inclination = geocentric latitude",
          fabs(is.inclination - geocLat) < 1.0e-12);
    Check("CircularSpeed", Same(CircularSpeed(ELL_PZ9011.mu, 7.0e6), sqrt(ELL_PZ9011.mu / 7.0e6)));
}

int main()
{
    RNG rng(20261017ULL);
    printf("gnclib %s: utilities, tables, observers, gusts, wind, orbit\n", GNCLIB_VERSION);
    TestUtil(&rng);
    TestTable(&rng);
    TestEso(&rng);
    TestNotch();
    TestGust();
    TestWindProfile();
    TestOrbit(&rng);
    printf("TEST_COMMON: %s (%d failed)\n", g_bad == 0 ? "PASS" : "FAIL", g_bad);
    return g_bad == 0 ? 0 : 1;
}

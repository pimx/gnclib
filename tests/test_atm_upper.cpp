// test_atm_upper.cpp -- gnclib 1.10.0: EnvAtmUpper, GOST R 25645.166-2004.
// (a) rho_n, K'0 ... K'4, K''4 against tables 4-11 of the standard;
// (b) formula (1) against an independent implementation (400 points);
// (c) Sun and sidereal time against astropy;
// (d) the diurnal maximum at 12 h + phi1 local solar time;
// (e) Smooth: Gost outside the windows to the bit, no steps at 120 km,
//     at the group bounds and across the F0 levels;
// (f) F81, kpp, Kp(Ap), indices from series;
// (g) above 1500 km: continuous, falling; the height found = given.
// Build: g++ -std=c++17 -O2 -Iinclude tests/test_atm_upper.cpp
//        -o build/test_atm_upper
#include "gnclib_version.h"
#include "lib_atm_upper.h"
#include "test_atm_upper_data.h"
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

// A point at a geocentric latitude and longitude, radius from the height
// above a sphere of 6371 km (formula (1) uses the direction only; the
// height is passed apart).
static CfEcefPoint PointAt(double lat, double lon, double hM)
{
    double r = 6371000.0 + hM;
    CfEcefPoint p { };
    p.x = r * cos(lat) * cos(lon);
    p.y = r * cos(lat) * sin(lon);
    p.z = r * sin(lat);
    return p;
}

// Half a unit of the third significant digit of v (table 4 is printed so).
static double HalfUnit3(double v)
{
    double e = floor(log10(fabs(v)));
    return 0.5 * pow(10.0, e - 2.0);
}

// (a) Tables 4-11.
static void TestTables()
{
    double wR = 0.0;
    double w[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    int i = 0;
    for (i = 0; i < 70; i = i + 1)
    {
        double h = 120.0 + 20.0 * (double)i;
        int l = 0;
        for (l = 0; l < 7; l = l + 1)
        {
            double r = EnvAtmUpper::NightDensity(l, h) - TAB_RHON[i][l];
            double u = fabs(r) / HalfUnit3(TAB_RHON[i][l]);
            wR = u > wR ? u : wR;
            double m[5] = { EnvAtmUpper::KP0(l, h), EnvAtmUpper::KP1(l, h), EnvAtmUpper::KP2(l, h),
                            EnvAtmUpper::KP3(l, h), EnvAtmUpper::KP4(l, h) };
            const double* t[5] = { TAB_K0[i], TAB_K1[i], TAB_K2[i], TAB_K3[i], TAB_K4[i] };
            int k = 0;
            for (k = 0; k < 5; k = k + 1)
            {
                //  Table 7 prints 4.466 at 780 km, F0 = 125 (neighbours 2.442, 2.487).
                if (k == 2 && i == 33 && l == 2)
                {
                    continue;
                }
                double d = fabs(fabs(m[k]) - fabs(t[k][l]));
                w[k] = d > w[k] ? d : w[k];
            }
        }
    }
    double wpp = 0.0;
    int j = 0;
    for (j = 0; j < 22; j = j + 1)
    {
        int l = 0;
        for (l = 0; l < 7; l = l + 1)
        {
            double d1 = fabs(EnvAtmUpper::KPP4(l, TAB_K4PP_DAY_KP[j], false) - TAB_K4PP_DAY[j][l]);
            double d2 = fabs(EnvAtmUpper::KPP4(l, TAB_K4PP_3H_KP[j], true) - TAB_K4PP_3H[j][l]);
            wpp = d1 > wpp ? d1 : wpp;
            wpp = d2 > wpp ? d2 : wpp;
        }
    }
    printf("  rho_n: worst %.2f of half the 3rd digit; K'0 %.4f K'1 %.4f K'2 %.4f K'3 %.4f K'4 %.4f; K''4 %.4f\n",
           wR, w[0], w[1], w[2], w[3], w[4], wpp);
    Check("Table 4: night density within the printed digits", wR <= 1.02);
    Check("Tables 5, 7-9: K'0, K'2, K'3, K'4 within 0.0008", w[0] < 8.0e-4 && w[2] < 8.0e-4 && w[3] < 8.0e-4 && w[4] < 8.0e-4);
    Check("Table 6: K'1 within 0.002 (6-digit coefficients, F0 = 175)", w[1] < 2.0e-3);
    Check("Tables 10, 11: K''4 within 0.0005", wpp <= 5.0e-4 + 1.0e-12);
}

// (b) Formula (1) against the independent implementation.
static void TestReference()
{
    EnvAtmUpper atm;
    int n = (int)(sizeof(REF_CASES) / sizeof(REF_CASES[0]));
    double worst = 0.0;
    int i = 0;
    for (i = 0; i < n; i = i + 1)
    {
        const double* c = REF_CASES[i];
        CfEcefPoint p { };
        p.x = c[0];
        p.y = c[1];
        p.z = c[2];
        UpperAtmTime tm { };
        tm.ra = c[4];
        tm.dec = c[5];
        tm.sStar = c[6];
        tm.tUt = c[7];
        tm.day = c[8];
        UpperAtmIndices ix { };
        ix.f107 = c[9];
        ix.f81 = c[10];
        ix.kp = c[11];
        ix.kp3h = c[12] > 0.5;
        double rho = atm.Density(p, c[3], tm, ix);
        double d = fabs(rho / c[13] - 1.0);
        worst = d > worst ? d : worst;
    }
    printf("  %d points (0-1500 km, F81 65-280, both Kp kinds): worst relative %.2e\n", n, worst);
    Check("Formula (1) and appendix A: the independent code within 1e-12", worst < 1.0e-12);
}

// (c) Sun and sidereal time.
static void TestTime()
{
    int n = (int)(sizeof(REF_EPH) / sizeof(REF_EPH[0]));
    double wRa = 0.0;
    double wDec = 0.0;
    double wG = 0.0;
    int badDay = 0;
    int i = 0;
    for (i = 0; i < n; i = i + 1)
    {
        UpperAtmTime tm = UpperAtmTimeOf(REF_EPH[i][0]);
        double dra = fabs(remainder(tm.ra - REF_EPH[i][1], 2.0 * PI));
        double ddec = fabs(tm.dec - REF_EPH[i][2]);
        double dg = fabs(remainder(tm.sStar - REF_EPH[i][3], 2.0 * PI));
        wRa = dra > wRa ? dra : wRa;
        wDec = ddec > wDec ? ddec : wDec;
        wG = dg > wG ? dg : wG;
        badDay = badDay + ((int)floor(tm.day) == (int)REF_EPH[i][4] ? 0 : 1);
    }
    printf("  Sun against astropy (8 dates 2000-2028): ra %.4f deg, dec %.4f deg; S* %.2e rad\n",
           wRa * 180.0 / PI, wDec * 180.0 / PI, wG);
    //  S*: astropy goes by UT1 (UT1 - UTC up to 0.9 s = 6.6e-5 rad), here UTC is taken for UT.
    Check("Sun within 0.01 deg, sidereal time within 1e-4 rad (UT1-UTC)", wRa < 1.75e-4 && wDec < 1.75e-4 && wG < 1.0e-4);
    Check("Day of the year", badDay == 0);
    UpperAtmTime t0 = UpperAtmTimeOf(MjdOfDate(2026, 10, 21) + 0.25);
    Check("2026-10-21 06:00 UT: tUt 21600 s, day 294.25", fabs(t0.tUt - 21600.0) < 1.0e-6 && fabs(t0.day - 294.25) < 1.0e-9);
}

// (d) The diurnal maximum: 12 h + phi1 local solar time at the Sun's latitude.
static void TestDiurnal()
{
    EnvAtmUpper atm;
    UpperAtmTime tm = UpperAtmTimeOf(MjdOfDate(2026, 6, 21) + 0.5);
    UpperAtmIndices ix { };
    ix.f107 = 150.0;
    ix.f81 = 150.0;
    ix.kp = 2.667;
    ix.kp3h = false;
    double lonSun = tm.ra - tm.sStar - EnvAtmUpper::OMEGA_3 * tm.tUt;
    double best = 0.0;
    double bestLon = 0.0;
    double worstRho = 1.0e30;
    int k = 0;
    for (k = 0; k < 1440; k = k + 1)
    {
        double lon = lonSun + (double)k * PI / 720.0;
        double rho = atm.Density(PointAt(tm.dec, lon, 400000.0), 400000.0, tm, ix);
        if (rho > best)
        {
            best = rho;
            bestLon = (double)k * PI / 720.0;
        }
        worstRho = rho < worstRho ? rho : worstRho;
    }
    double lst = 12.0 + bestLon * 12.0 / PI;
    printf("  400 km, F81 = 150: maximum at %.2f h local solar time (12 + phi1: %.2f h), max/min %.2f\n",
           lst, 12.0 + 0.5585 * 12.0 / PI, best / worstRho);
    Check("Diurnal maximum at 12 h + phi1 (within 0.05 h)", fabs(lst - (12.0 + 0.5585 * 12.0 / PI)) < 0.05);
}

// Largest |ln rho(h + dh) - ln rho(h)| / dh over [h0, h1] (heights in m).
static double MaxLogSlope(const EnvAtmUpper& a, const CfEcefPoint& dir, double h0, double h1, double dh,
                          const UpperAtmTime& tm, const UpperAtmIndices& ix)
{
    double w = 0.0;
    double r0 = sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
    double h = h0;
    double f0 = log(a.Density(dir, h, tm, ix));
    while (h < h1)
    {
        double f1 = log(a.Density(dir, h + dh, tm, ix));
        double g = fabs(f1 - f0) / dh;
        w = g > w ? g : w;
        f0 = f1;
        h = h + dh;
    }
    (void)r0;
    return w;
}

// (e) Smooth mode.
static void TestSmooth()
{
    EnvAtmUpper gost;
    EnvAtmUpper smooth(AtmMode::Smooth);
    UpperAtmTime tm = UpperAtmTimeOf(MjdOfDate(2026, 4, 5) + 0.6);
    UpperAtmIndices ix { };
    ix.f107 = 160.0;
    ix.f81 = 150.0;
    ix.kp = 3.0;
    ix.kp3h = true;
    CfEcefPoint dir = PointAt(0.3, 1.0, 0.0);

    //  Outside every window -- Gost to the bit.
    long bad = 0;
    long n = 0;
    double h = 0.0;
    for (h = 130000.0; h < 1500000.0; h = h + 997.0)
    {
        double hk = h / 1000.0;
        bool inside = fabs(hk - 500.0) < 10.0 || fabs(hk - 600.0) < 10.0 || fabs(hk - 640.0) < 10.0
            || fabs(hk - 660.0) < 10.0 || fabs(hk - 700.0) < 10.0 || fabs(hk - 740.0) < 10.0
            || fabs(hk - 760.0) < 10.0 || fabs(hk - 780.0) < 10.0 || fabs(hk - 800.0) < 10.0
            || fabs(hk - 820.0) < 10.0 || fabs(hk - 860.0) < 10.0 || fabs(hk - 900.0) < 10.0
            || fabs(hk - 920.0) < 10.0 || fabs(hk - 980.0) < 10.0 || fabs(hk - 1000.0) < 10.0;
        if (inside)
        {
            continue;
        }
        bad = bad + (Same(gost.Density(dir, h, tm, ix), smooth.Density(dir, h, tm, ix)) ? 0 : 1);
        n = n + 1;
    }
    Check("Smooth: outside the windows Gost to the bit", bad == 0 && n > 1000);

    //  Steps in height: 120 km and every group bound of every level.
    double wG = 0.0;
    double wS = 0.0;
    int lev = 0;
    for (lev = 0; lev < 7; lev = lev + 1)
    {
        UpperAtmIndices iy = ix;
        iy.f81 = EnvAtmUpper::F0(lev);
        iy.f107 = iy.f81 + 10.0;
        double g1 = MaxLogSlope(gost, dir, 110000.0, 130000.0, 1.0, tm, iy);
        double s1 = MaxLogSlope(smooth, dir, 110000.0, 130000.0, 1.0, tm, iy);
        double g2 = MaxLogSlope(gost, dir, 480000.0, 1020000.0, 1.0, tm, iy);
        double s2 = MaxLogSlope(smooth, dir, 480000.0, 1020000.0, 1.0, tm, iy);
        wG = g1 > wG ? g1 : wG;
        wG = g2 > wG ? g2 : wG;
        wS = s1 > wS ? s1 : wS;
        wS = s2 > wS ? s2 : wS;
    }
    printf("  |d ln rho / dh| over 1 m, 110-130 and 480-1020 km, all levels: Gost %.2e, Smooth %.2e 1/m\n", wG, wS);
    Check("Gost: the steps are there (above 1e-3 1/m)", wG > 1.0e-3);
    //  The profile itself falls by up to 1.9e-4 1/m near 110 km (appendix A: k1 = -0.18553 1/km).
    Check("Smooth: no steps in height (below 3e-4 1/m)", wS < 3.0e-4);

    //  Steps in F81 across the level middles, 400 km.
    double wFG = 0.0;
    double wFS = 0.0;
    double f = 0.0;
    double pG = 0.0;
    double pS = 0.0;
    for (f = 70.0; f <= 260.0; f = f + 0.01)
    {
        UpperAtmIndices iy = ix;
        iy.f81 = f;
        iy.f107 = f;
        double rg = log(gost.Density(dir, 400000.0, tm, iy));
        double rs = log(smooth.Density(dir, 400000.0, tm, iy));
        if (f > 70.0)
        {
            wFG = fabs(rg - pG) > wFG ? fabs(rg - pG) : wFG;
            wFS = fabs(rs - pS) > wFS ? fabs(rs - pS) : wFS;
        }
        pG = rg;
        pS = rs;
    }
    printf("  F81 70-260 by 0.01 at 400 km: largest |d ln rho| Gost %.2e, Smooth %.2e\n", wFG, wFS);
    Check("Smooth: no steps across the F0 levels (below 2e-3 per 0.01)", wFS < 2.0e-3 && wFG > 2.0e-2);
}

// (f) Indices.
static void TestIndices()
{
    static double f[200];
    static double kd[200];
    static double k3[1600];
    static double kpp[1600];
    int i = 0;
    for (i = 0; i < 200; i = i + 1)
    {
        f[i] = 100.0 + 0.5 * (double)i;
        kd[i] = 1.0 + 0.01 * (double)i;
    }
    for (i = 0; i < 1600; i = i + 1)
    {
        k3[i] = (i / 8) % 2 == 0 ? 2.0 : 5.0;
    }
    KppOf(k3, 1600, kpp);
    //  F81 of a linear series: the value at the weighted mean lag.
    double sw = 0.0;
    double sj = 0.0;
    int j = 0;
    for (j = -80; j <= 0; j = j + 1)
    {
        double w = 1.0 + 0.5 * (double)j / 80.0;
        sw = sw + w;
        sj = sj + w * (double)j;
    }
    double expect = 100.0 + 0.5 * (120.0 + sj / sw);
    Check("F81 of a linear series: the weighted mean", fabs(F81Of(f, 120) - expect) < 1.0e-12);
    Check("kpp: first = kp; rise by 0.7 of the step, fall by 0.3",
          kpp[0] == 2.0 && fabs(kpp[8] - (2.0 + 0.7 * 3.0)) < 1.0e-12 && fabs(kpp[16] - (kpp[15] - 0.3 * (kpp[15] - 2.0))) < 1.0e-12);
    Check("Kp of Ap: table A.1 rows and between",
          KpOfAp(0.0) == 0.0 && fabs(KpOfAp(15.0) - 3.0) < 1.0e-12 && fabs(KpOfAp(400.0) - 9.0) < 1.0e-12
          && fabs(KpOfAp(10.5) - (7.0 + 0.5) / 3.0) < 1.0e-12);
    UpperAtmSeries s { };
    s.mjd0 = 61000.0;
    s.nDays = 200;
    s.f107 = f;
    s.kpDay = kd;
    s.kpp3h = nullptr;
    UpperAtmIndices ix { };
    double t = 61150.3;
    bool ok = UpperAtmIndicesAt(s, t, &ix);
    double fx = 100.0 + 0.5 * (t - EnvAtmUpper::DT_F - (61000.0 + 20.0 / 24.0));
    double kx = 1.0 + 0.01 * (t - EnvAtmUpper::DT_KP_DAY - 61000.5);
    Check("Indices at a moment: delayed and interpolated (daily Kp)",
          ok && fabs(ix.f107 - fx) < 1.0e-9 && fabs(ix.kp - kx) < 1.0e-9 && !ix.kp3h);
    Check("Indices: F81 needs 80 days before", !UpperAtmIndicesAt(s, 61050.0, &ix));
    s.kpp3h = kpp;
    ok = UpperAtmIndicesAt(s, 61100.0 + 0.3125 + 0.25, &ix);
    Check("Indices: 3-hour kpp at the middle of its interval, 0.25 d delay", ok && ix.kp3h && Same(ix.kp, kpp[8 * 100 + 2]));
}

// (g) Above 1500 km; the height found.
static void TestTop()
{
    EnvAtmUpper atm;
    UpperAtmTime tm = UpperAtmTimeOf(MjdOfDate(2026, 1, 10) + 0.1);
    UpperAtmIndices ix { };
    ix.f107 = 120.0;
    ix.f81 = 125.0;
    ix.kp = 1.0;
    ix.kp3h = false;
    CfEcefPoint dir = PointAt(-0.5, 2.0, 0.0);
    double r = sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
    CfEcefPoint pa = dir;
    double k = (r + 1500000.0) / r;
    pa.x = dir.x * k;
    pa.y = dir.y * k;
    pa.z = dir.z * k;
    CfEcefPoint pb = dir;
    double kb = (r + 1500000.001) / r;
    pb.x = dir.x * kb;
    pb.y = dir.y * kb;
    pb.z = dir.z * kb;
    double a = atm.Density(pa, 1500000.0, tm, ix);
    double b = atm.Density(pb, 1500000.001, tm, ix);
    CfEcefPoint pc = dir;
    double kc = (r + 2000000.0) / r;
    pc.x = dir.x * kc;
    pc.y = dir.y * kc;
    pc.z = dir.z * kc;
    double c = atm.Density(pc, 2000000.0, tm, ix);
    printf("  1500 km: %.4e, 1 mm above: rel %.1e; 2000 km: %.3e kg/m^3\n", a, b / a - 1.0, c);
    Check("Above 1500 km: continuous and falling", fabs(b / a - 1.0) < 1.0e-6 && c < a && c > 0.0);
    CfEcefPoint q = PointAt(0.7, -1.2, 420000.0);
    Check("Density(p): the height above PZ-90.11 found",
          Same(atm.Density(q, tm, ix), atm.Density(q, EcpToAlt(q, ELL_PZ9011), tm, ix)));
}

// Entry point.
int main()
{
    printf("TEST_ATM_UPPER (gnclib %s)\n", GNCLIB_VERSION);
    TestTables();
    TestReference();
    TestTime();
    TestDiurnal();
    TestSmooth();
    TestIndices();
    TestTop();
    printf("TEST_ATM_UPPER %s\n", g_bad == 0 ? "PASS" : "FAIL");
    return g_bad == 0 ? 0 : 1;
}

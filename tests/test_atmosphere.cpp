// test_atmosphere.cpp -- gnclib 1.8.0/1.9.0: EnvAtmosphere against GOST 4401-81.
// (a) table 4 of the standard (geometric height), every row available
//     from 80 to 250 km, plus sea level and two sound-speed rows;
// (b) the steps the standard itself makes at 120, 150, 200 km (printed);
// (c) pressure never rises with height, -2 ... 250 km every metre;
// (d) the temperature offset keeps the pressure and p / (rho T);
// (e) State and MolarMass against the separate calls, to the bit;
// (f) 1.9.0 Smooth mode: Gost outside the windows to the bit, no steps,
//     C2 density across the knots, off the standard by less than a step;
// (g) RK4 order of a fall with drag through the knots, both modes.
// Build: g++ -std=c++17 -O2 -Iinclude tests/test_atmosphere.cpp
//        -o build/test_atmosphere
#include "gnclib_version.h"
#include "lib_atmosphere.h"
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

// One row of GOST 4401-81 table 4: geometric height, m; temperature, K;
// pressure, Pa; density, kg/m^3; molar mass, kg/kmol (0 - not checked).
// side: 0 - the height itself; -1 - the limit from below (the table gives
// the knot the value of the interval below it).
struct GostRow
{
    double h;
    double t;
    double p;
    double rho;
    double m;
    int side;
};

static const GostRow ROWS[] =
{
    { 0.0,      288.150, 1.01325e5,  1.22500,    28.964420, 0 },
    { 80000.0,  198.639, 1.05247,    1.84580e-5, 28.964420, 0 },
    { 90000.0,  186.650, 1.83140e-1, 3.41817e-6, 28.964420, 0 },
    { 91000.0,  186.650, 1.53292e-1, 2.86108e-6, 28.964420, 0 },
    { 92000.0,  186.650, 1.28308e-1, 2.39477e-6, 28.964420, 0 },
    { 93000.0,  186.650, 1.07416e-1, 2.00484e-6, 28.964420, 0 },
    { 94000.0,  186.650, 8.99258e-2, 1.67840e-6, 28.964420, 0 },
    { 95000.0,  186.622, 7.52834e-2, 1.40510e-6, 28.960, 0 },
    { 96000.0,  188.240, 6.30558e-2, 1.16617e-6, 28.945, 0 },
    { 97000.0,  190.917, 5.29545e-2, 9.64450e-7, 28.910, 0 },
    { 98000.0,  192.912, 4.45890e-2, 7.99924e-7, 28.774, 0 },
    { 99000.0,  194.775, 3.76413e-2, 6.65312e-7, 28.623, 0 },
    { 100000.0, 196.605, 3.18606e-2, 5.54951e-7, 28.472, 0 },
    { 101000.0, 198.404, 2.70326e-2, 4.64111e-7, 28.321, 0 },
    { 102000.0, 200.172, 2.29898e-2, 3.89128e-7, 28.170, 0 },
    { 103000.0, 201.911, 1.95962e-2, 3.27067e-7, 28.019, 0 },
    { 104000.0, 203.618, 1.67431e-2, 2.75614e-7, 27.868, 0 },
    { 110000.0, 255.487, 7.35900e-3, 9.34035e-8, 26.961, 0 },
    { 114000.0, 288.427, 4.70615e-3, 5.17245e-8, 26.357, 0 },
    { 115000.0, 296.377, 4.24849e-3, 4.51814e-8, 26.206, 0 },
    { 116000.0, 304.207, 3.84826e-3, 3.96418e-8, 26.055, 0 },
    { 117000.0, 311.926, 3.49640e-3, 3.49222e-8, 25.904, 0 },
    { 118000.0, 319.535, 3.18584e-3, 3.08814e-8, 25.752, 0 },
    { 119000.0, 327.033, 2.91073e-3, 2.74061e-8, 25.601, 0 },
    { 120000.0, 334.417, 2.66618e-3, 2.44041e-8, 25.450, -1 },
    { 121000.0, 345.68,  2.50215e-3, 2.20814e-8, 25.36, 0 },
    { 122000.0, 356.94,  2.36291e-3, 2.01267e-8, 25.28, 0 },
    { 123000.0, 368.20,  2.22898e-3, 1.83438e-8, 25.19, 0 },
    { 124000.0, 379.46,  2.10075e-3, 1.67199e-8, 25.11, 0 },
    { 125000.0, 390.72,  1.97851e-3, 1.52430e-8, 25.03, 0 },
    { 126000.0, 401.97,  1.86245e-3, 1.39015e-8, 24.95, 0 },
    { 130000.0, 447.01,  1.46076e-3, 9.68064e-9, 24.63, 0 },
    { 140000.0, 559.60,  8.26375e-4, 4.24614e-9, 23.91, 0 },
    { 145000.0, 593.60,  6.08841e-4, 2.90864e-9, 23.58, 0 },
    { 146000.0, 600.40,  5.72689e-4, 2.69768e-9, 23.51, 0 },
    { 147000.0, 607.20,  5.38618e-4, 2.50210e-9, 23.45, 0 },
    { 148000.0, 614.00,  5.06620e-4, 2.32127e-9, 23.39, 0 },
    { 149000.0, 620.80,  4.76771e-4, 2.15494e-9, 23.33, 0 },
    { 150000.0, 627.60,  4.49233e-4, 2.00329e-9, 23.27, -1 },
    { 151000.0, 634.40,  4.31065e-4, 1.89681e-9, 23.21, 0 },
    { 152000.0, 641.20,  4.14031e-4, 1.79797e-9, 23.15, 0 },
    { 153000.0, 648.00,  3.97769e-4, 1.70494e-9, 23.09, 0 },
    { 154000.0, 654.80,  3.82257e-4, 1.61741e-9, 23.04, 0 },
    { 155000.0, 661.60,  3.67473e-4, 1.53510e-9, 22.98, 0 },
    { 160000.0, 695.60,  3.03620e-4, 1.19204e-9, 22.71, 0 },
    { 170000.0, 735.30,  2.08053e-4, 7.55732e-10, 22.21, 0 },
    { 175000.0, 755.15,  1.76153e-4, 6.16603e-10, 21.98, 0 },
    { 180000.0, 775.00,  1.50574e-4, 5.08478e-10, 21.76, 0 },
    { 190000.0, 814.70,  1.10866e-4, 3.49499e-10, 21.35, 0 },
    { 191000.0, 818.67,  1.07561e-4, 3.36823e-10, 21.31, 0 },
    { 192000.0, 822.64,  1.04386e-4, 3.24717e-10, 21.28, 0 },
    { 193000.0, 826.61,  1.01351e-4, 3.13199e-10, 21.24, 0 },
    { 194000.0, 830.58,  9.84678e-5, 3.02294e-10, 21.20, 0 },
    { 195000.0, 834.55,  9.57500e-5, 2.92032e-10, 21.16, 0 },
    { 196000.0, 838.52,  9.32141e-5, 2.82451e-10, 21.13, 0 },
    { 197000.0, 842.49,  9.08788e-5, 2.73595e-10, 21.09, 0 },
    { 198000.0, 846.46,  8.87655e-5, 2.65513e-10, 21.05, 0 },
    { 199000.0, 850.43,  8.68979e-5, 2.58262e-10, 21.01, 0 },
    { 200000.0, 854.40,  8.53026e-5, 2.51904e-10, 20.98, -1 },
    { 201000.0, 856.15,  8.23175e-5, 2.42171e-10, 20.94, 0 },
    { 202000.0, 857.90,  8.00431e-5, 2.34594e-10, 20.91, 0 },
    { 203000.0, 859.65,  7.78409e-5, 2.27283e-10, 20.87, 0 },
    { 204000.0, 861.40,  7.57087e-5, 2.20229e-10, 20.83, 0 },
    { 205000.0, 863.15,  7.36446e-5, 2.13424e-10, 20.80, 0 },
    { 206000.0, 864.90,  7.16466e-5, 2.06859e-10, 20.76, 0 },
    { 207000.0, 866.65,  6.97127e-5, 2.00526e-10, 20.73, 0 },
    { 208000.0, 868.40,  6.78409e-5, 1.94416e-10, 20.69, 0 },
    { 209000.0, 870.15,  6.60295e-5, 1.88523e-10, 20.66, 0 },
    { 210000.0, 871.90,  6.42764e-5, 1.82838e-10, 20.62, 0 },
    { 211000.0, 873.65,  6.25800e-5, 1.77353e-10, 20.59, 0 },
    { 212000.0, 875.40,  6.09383e-5, 1.72062e-10, 20.55, 0 },
    { 220000.0, 889.40,  4.95697e-5, 1.35889e-10, 20.27, 0 },
    { 230000.0, 906.90,  3.88854e-5, 1.02724e-10, 19.92, 0 },
    { 240000.0, 924.40,  3.09161e-5, 7.86493e-11, 19.55, 0 },
    { 250000.0, 941.90,  2.47564e-5, 6.05679e-11, 19.16, 0 }
};

// Relative miss.
static double Rel(double a, double b)
{
    return fabs(a / b - 1.0);
}

// Molar mass printed with 2, 3 or 6 decimals: half of its last digit;
// 94-97 km: 0.001 (the square-root formula of GOST section 4 misses the
// 3-decimal table by up to 0.0009 there; pressure and density do not use M
// below 120 km, kinetic temperature moves by at most 0.005 K).
static double MTol(double h, double m)
{
    if (h > 94000.0 && h < 97000.0)
    {
        return 0.001;
    }
    double k = m * 100.0;
    if (fabs(k - floor(k + 0.5)) < 1.0e-9)
    {
        return 0.0051;
    }
    k = m * 1000.0;
    if (fabs(k - floor(k + 0.5)) < 1.0e-9)
    {
        return 0.00051;
    }
    return 1.0e-6;
}

// (a) Table 4.
static void TestTable()
{
    int n = (int)(sizeof(ROWS) / sizeof(ROWS[0]));
    double wt = 0.0;
    double wp = 0.0;
    double wr = 0.0;
    double hT = 0.0;
    double hP = 0.0;
    double hR = 0.0;
    int badM = 0;
    int i = 0;
    for (i = 0; i < n; i = i + 1)
    {
        const GostRow& g = ROWS[i];
        double h = g.h;
        if (g.side < 0)
        {
            h = g.h - 1.0e-6;
        }
        AtmState s = ATM_STD.State(h);
        double dt = fabs(s.t - g.t);
        double dp = Rel(s.p, g.p);
        double dr = Rel(s.rho, g.rho);
        if (dt > wt)
        {
            wt = dt;
            hT = g.h;
        }
        if (dp > wp)
        {
            wp = dp;
            hP = g.h;
        }
        if (dr > wr)
        {
            wr = dr;
            hR = g.h;
        }
        if (g.m > 0.0 && fabs(s.m - g.m) > MTol(g.h, g.m))
        {
            printf("  M at %.0f m: %.6f, table %.6f\n", g.h, s.m, g.m);
            badM = badM + 1;
        }
    }
    printf("  %d rows: worst T %.4f K (%.0f m), p %.2e (%.0f m), rho %.2e (%.0f m)\n",
           n, wt, hT, wp, hP, wr, hR);
    Check("GOST table 4: temperature within 0.006 K", wt < 0.006);
    Check("GOST table 4: pressure within 1e-4", wp < 1.0e-4);
    Check("GOST table 4: density within 1e-4", wr < 1.0e-4);
    Check("GOST table 4: molar mass within the printed digits", badM == 0);

    //  Sound speed, GOST table 2 (geometric 26 and 27 km) and sea level.
    double a0 = ATM_STD.SoundSpeed(0.0);
    double a26 = ATM_STD.SoundSpeed(26000.0);
    double a27 = ATM_STD.SoundSpeed(27000.0);
    printf("  sound speed: %.4f %.4f %.4f m/s (table 340.294 299.056 299.722)\n", a0, a26, a27);
    Check("GOST: sound speed within 0.001 m/s",
          fabs(a0 - 340.294) < 1.0e-3 && fabs(a26 - 299.056) < 1.0e-3 && fabs(a27 - 299.722) < 1.0e-3);
}

// (b) The steps of the standard at its knots: printed, and bounded.
static void TestSteps()
{
    double knots[3] = { 120000.0, 150000.0, 200000.0 };
    double worst = 0.0;
    int k = 0;
    for (k = 0; k < 3; k = k + 1)
    {
        AtmState lo = ATM_STD.State(knots[k] - 1.0e-6);
        AtmState hi = ATM_STD.State(knots[k]);
        double sp = hi.p / lo.p - 1.0;
        double sr = hi.rho / lo.rho - 1.0;
        printf("  step at %.0f km: T %+.4f K, p %+.3e, rho %+.3e, M %+.2e\n",
               knots[k] / 1000.0, hi.t - lo.t, sp, sr, hi.m - lo.m);
        if (fabs(sp) > worst)
        {
            worst = fabs(sp);
        }
    }
    Check("Steps of the standard at 120/150/200 km below 0.8 %", worst < 0.008);
    double m94 = ATM_STD.MolarMass(94000.0);
    double m94u = ATM_STD.MolarMass(94000.0 + 1.0e-3);
    double m97 = ATM_STD.MolarMass(97000.0);
    double m97u = ATM_STD.MolarMass(97000.0 + 1.0e-3);
    double m975 = ATM_STD.MolarMass(97500.0);
    double m975u = ATM_STD.MolarMass(97500.0 + 1.0e-3);
    printf("  molar mass steps at 94, 97, 97.5 km: %.2e %.2e %.2e kg/kmol\n",
           m94u - m94, m97u - m97, m975u - m975);
    Check("Molar mass continuous at 94, 97, 97.5 km (1e-4)",
          fabs(m94u - m94) < 1.0e-4 && fabs(m97u - m97) < 1.0e-4 && fabs(m975u - m975) < 1.0e-4);
}

// (c) Pressure and density never rise with height.
static void TestMonotonic()
{
    long badP = 0;
    long badR = 0;
    double pPrev = 1.0e30;
    double rPrev = 1.0e30;
    long i = 0;
    for (i = -1999; i <= 250000; i = i + 1)
    {
        double h = (double)i;
        AtmState s = ATM_STD.State(h);
        badP = badP + (s.p <= pPrev ? 0 : 1);
        badR = badR + (s.rho <= rPrev ? 0 : 1);
        pPrev = s.p;
        rPrev = s.rho;
    }
    AtmState top = ATM_STD.State(250000.5);
    printf("  -2 ... 250 km every 1 m: pressure rises %ld times, density %ld\n", badP, badR);
    Check("Pressure and density never rise with height", badP == 0 && badR == 0);
    Check("Above 250 km: no air, temperature held", top.p == 0.0 && top.rho == 0.0 && top.t == 941.90);
    Check("Below -2 km geopotential: held", Same(ATM_STD.Pressure(-5000.0), ATM_STD.Pressure(-1999.3709472292348)));
}

// (d) Temperature offset: pressure kept, p M / (rho T) = R* everywhere.
static void TestOffset()
{
    EnvAtmosphere hot(15.0);
    double worst = 0.0;
    long badP = 0;
    long i = 0;
    for (i = 0; i <= 250; i = i + 1)
    {
        double h = (double)i * 1000.0 + 0.5;
        AtmState s = hot.State(h);
        badP = badP + (Same(s.p, ATM_STD.Pressure(h)) ? 0 : 1);
        double rs = s.p * s.m / (s.rho * s.t);
        double d = fabs(rs / R_STAR - 1.0);
        if (d > worst)
        {
            worst = d;
        }
    }
    printf("  dT = +15 K: pressure differs %ld times, p M / (rho T R*) - 1 worst %.2e\n", badP, worst);
    Check("Offset: pressure profile kept, state equation holds (1e-7)", badP == 0 && worst < 1.0e-7);
}

// (e) State and MolarMass against the separate calls.
static void TestBits()
{
    double dts[3] = { 0.0, 15.0, -20.0 };
    long bad = 0;
    long n = 0;
    int k = 0;
    for (k = 0; k < 6; k = k + 1)
    {
        EnvAtmosphere a(dts[k % 3], k < 3 ? AtmMode::Gost : AtmMode::Smooth);
        long i = 0;
        for (i = -25000; i <= 2600000; i = i + 7)
        {
            double h = (double)i * 0.1;
            AtmState s = a.State(h);
            bool ok = Same(s.t, a.Temperature(h)) && Same(s.p, a.Pressure(h))
                && Same(s.rho, a.Density(h)) && Same(s.m, a.MolarMass(h));
            bad = bad + (ok ? 0 : 1);
            n = n + 1;
        }
    }
    printf("  State against the separate calls (both modes): %ld heights, %ld differ\n", n, bad);
    Check("State: Temperature, Pressure, Density, MolarMass to the bit", bad == 0);
    double a = ATM_STD.SoundSpeed(50000.0);
    double ar = sqrt(1.4 * R_AIR * ATM_STD.Temperature(50000.0));
    Check("SoundSpeed below 94 km: sqrt(1.4 R T) to the bit", Same(a, ar));
}

// Largest |ln f(h + d) - ln f(h)| / d over [h0, h1], step d (f > 0).
static double MaxLogSlope(const EnvAtmosphere& a, bool dens, double h0, double h1, double d)
{
    double w = 0.0;
    double h = h0;
    double f0 = dens ? a.Density(h) : a.Pressure(h);
    while (h < h1)
    {
        double f1 = dens ? a.Density(h + d) : a.Pressure(h + d);
        double g = fabs(log(f1 / f0)) / d;
        if (g > w)
        {
            w = g;
        }
        f0 = f1;
        h = h + d;
    }
    return w;
}

// Largest |second difference of ln rho| / d^2 over [h0, h1], step d.
static double MaxLogCurv(const EnvAtmosphere& a, double h0, double h1, double d)
{
    double w = 0.0;
    double h = h0;
    while (h < h1)
    {
        double c = (log(a.Density(h + d)) - 2.0 * log(a.Density(h)) + log(a.Density(h - d))) / (d * d);
        if (fabs(c) > w)
        {
            w = fabs(c);
        }
        h = h + d;
    }
    return w;
}

// (f) Smooth mode.
static void TestSmooth()
{
    const EnvAtmosphere& g = ATM_STD;
    const EnvAtmosphere& s = ATM_STD_SMOOTH;
    double knots[3] = { 120000.0, 150000.0, 200000.0 };
    double wsm = EnvAtmosphere::W_SMOOTH;

    //  Outside the windows -- Gost to the bit.
    long bad = 0;
    long i = 0;
    for (i = -19990; i <= 2600000; i = i + 3)
    {
        double h = (double)i * 0.1;
        bool inside = false;
        int k = 0;
        for (k = 0; k < 3; k = k + 1)
        {
            if (fabs(h - knots[k]) < wsm)
            {
                inside = true;
            }
        }
        if (inside)
        {
            continue;
        }
        AtmState x = s.State(h);
        AtmState y = g.State(h);
        bool ok = Same(x.t, y.t) && Same(x.p, y.p) && Same(x.rho, y.rho) && Same(x.m, y.m);
        bad = bad + (ok ? 0 : 1);
    }
    Check("Smooth: outside the windows Gost to the bit", bad == 0);

    //  Steps: the largest log slope over 1 mm across each window, against
    //  the natural slope of the profile there (about 1 / scale height).
    double worstS = 0.0;
    double worstG = 0.0;
    double worstDev = 0.0;
    double worstCurvS = 0.0;
    double worstCurvG = 0.0;
    int k = 0;
    for (k = 0; k < 3; k = k + 1)
    {
        double h0 = knots[k] - wsm - 100.0;
        double h1 = knots[k] + wsm + 100.0;
        double ss = MaxLogSlope(s, true, h0, h1, 1.0e-3);
        double sp = MaxLogSlope(s, false, h0, h1, 1.0e-3);
        double gs = MaxLogSlope(g, true, h0, h1, 1.0e-3);
        double cs = MaxLogCurv(s, h0, h1, 1.0);
        double cg = MaxLogCurv(g, h0, h1, 1.0);
        double dev = 0.0;
        double h = h0;
        while (h < h1)
        {
            double d = fabs(s.Density(h) / g.Density(h) - 1.0);
            if (d > dev)
            {
                dev = d;
            }
            h = h + 1.0;
        }
        printf("  knot %.0f km: |dln rho/dh| Smooth %.2e (p %.2e), Gost %.2e 1/m; "
               "|d2ln rho| Smooth %.2e, Gost %.2e 1/m^2; off Gost %.2e\n",
               knots[k] / 1000.0, ss, sp, gs, cs, cg, dev);
        worstS = ss > sp ? (ss > worstS ? ss : worstS) : (sp > worstS ? sp : worstS);
        worstG = gs > worstG ? gs : worstG;
        worstCurvS = cs > worstCurvS ? cs : worstCurvS;
        worstCurvG = cg > worstCurvG ? cg : worstCurvG;
        worstDev = dev > worstDev ? dev : worstDev;
    }
    Check("Smooth: no steps -- log slope below 3e-4 1/m at the knots", worstS < 3.0e-4);
    Check("Gost: the steps are there (log slope over 1 mm above 1)", worstG > 1.0);
    Check("Smooth: density C2 -- second log difference below 1e-7 1/m^2", worstCurvS < 1.0e-7);
    Check("Smooth: off Gost by less than the step (0.76 %)", worstDev < 0.0076);

    //  Never rises; the blend ends match the pieces.
    long rises = 0;
    double pPrev = 1.0e30;
    double rPrev = 1.0e30;
    for (i = 110000; i <= 210000; i = i + 1)
    {
        AtmState x = s.State((double)i);
        rises = rises + (x.p <= pPrev && x.rho <= rPrev ? 0 : 1);
        pPrev = x.p;
        rPrev = x.rho;
    }
    Check("Smooth: pressure and density never rise (110-210 km, 1 m)", rises == 0);
}

// Fall with drag through the knots, RK4: dz/dt = v, dv/dt = -g + rho v^2 / (2 beta)
// (v < 0), from 210 km at -2000 m/s for 45 s, down to about 115 km: the
// knots and no layer base (below 120 km the density kinks at the layer
// bases in either mode); returns the final height.
static double FallRk4(const EnvAtmosphere& a, double dt)
{
    const double beta = 0.02;   // ballistic coefficient, kg/m^2 (light: drag counts)
    double z = 210000.0;
    double v = -2000.0;
    int n = (int)(45.0 / dt + 0.5);
    int i = 0;
    for (i = 0; i < n; i = i + 1)
    {
        double k1z = v;
        double k1v = -G0_STD + a.Density(z) * v * v / (2.0 * beta);
        double z2 = z + 0.5 * dt * k1z;
        double v2 = v + 0.5 * dt * k1v;
        double k2z = v2;
        double k2v = -G0_STD + a.Density(z2) * v2 * v2 / (2.0 * beta);
        double z3 = z + 0.5 * dt * k2z;
        double v3 = v + 0.5 * dt * k2v;
        double k3z = v3;
        double k3v = -G0_STD + a.Density(z3) * v3 * v3 / (2.0 * beta);
        double z4 = z + dt * k3z;
        double v4 = v + dt * k3v;
        double k4z = v4;
        double k4v = -G0_STD + a.Density(z4) * v4 * v4 / (2.0 * beta);
        z = z + dt * (k1z + 2.0 * k2z + 2.0 * k3z + k4z) / 6.0;
        v = v + dt * (k1v + 2.0 * k2v + 2.0 * k3v + k4v) / 6.0;
    }
    return z;
}

// (g) What the integrator sees: order of RK4 across the knots.
static void TestIntegrator()
{
    double dts[3] = { 0.9, 0.45, 0.225 };
    double ord[2] = { 0.0, 0.0 };
    int m = 0;
    for (m = 0; m < 2; m = m + 1)
    {
        const EnvAtmosphere& a = m == 0 ? ATM_STD : ATM_STD_SMOOTH;
        double ref = FallRk4(a, 45.0 / 3200.0);
        double e[3] = { 0.0, 0.0, 0.0 };
        int k = 0;
        for (k = 0; k < 3; k = k + 1)
        {
            e[k] = fabs(FallRk4(a, dts[k]) - ref);
        }
        ord[m] = log(e[1] / e[2]) / log(2.0);
        printf("  RK4 fall 210 -> %.0f km, %s: height error at dt 0.9/0.45/0.225 s: %.2e %.2e %.2e m, order %.2f\n",
               ref / 1000.0, m == 0 ? "Gost  " : "Smooth", e[0], e[1], e[2], ord[m]);
    }
    Check("Smooth: RK4 keeps order above 3.5 across the knots", ord[1] > 3.5);
}

// Entry point.
int main()
{
    printf("TEST_ATMOSPHERE (gnclib %s)\n", GNCLIB_VERSION);
    TestTable();
    TestSteps();
    TestMonotonic();
    TestOffset();
    TestBits();
    TestSmooth();
    TestIntegrator();
    printf("TEST_ATMOSPHERE %s\n", g_bad == 0 ? "PASS" : "FAIL");
    return g_bad == 0 ? 0 : 1;
}

//======================================================================
//  lib_atm_upper.h
//  EnvAtmUpper - density of the upper atmosphere by GOST R 25645.166-2004
//  ("Atmosphere of the Earth, upper. Density model for ballistic support
//  of flights of artificial Earth satellites"), 120 ... 1500 km, with its
//  appendix A for the heights below 120 km. Explicit include (not lib.h).
//
//  MODEL (section 5 of the standard), h -- height above the ellipsoid, km:
//    rho = rho_n K0 (1 + K1 + K2 + K3 + K4),                         (1)
//    rho_n = rho_0 exp(a0 + a1 h + ... + a6 h^6), rho_0 = 1.58868e-8,
//    K0 = 1 + (l0 + l1 h + ... + l4 h^4) (F81 - F0) / F0,
//    K1 = (c0 + ... + c4 h^4) cos^(n0 + n1 h + n2 h^2)(phi / 2),
//      cos phi = (z sin dS + cos dS (x cos beta + y sin beta)) / r,
//      beta = aS - S* - w t + phi1,
//    K2 = (d0 + ... + d4 h^4) A(d),
//    K3 = (b0 + ... + b4 h^4) (F10.7 - F81) / (F81 + |F10.7 - F81|),
//    K4 = (e0 + ... + e4 h^4) (e5 + e6 Kp + e7 Kp^2 + e8 Kp^3)
//         (3-hour indices: et5 ... et8 and the modified kpp).
//  Coefficients: table 2 (first height range) and table 3 (second one,
//  each group from its own lower bound: a from 500 km, b, c, e, l from
//  600 ... 1000 km depending on F0, d never); the level F0 (75, 100,
//  125, 150, 175, 200, 250) nearest to F81. The second range is taken
//  ABOVE its bound: tables 4-11 of the standard give the bound itself
//  by the first range (checked). Below 120 km -- appendix A.4,
//  rho = a0i exp(k1i (h - hi) + k2i (h - hi)^2), table A.2 (static).
//  Above 1500 km (outside the standard) -- the value at 1500 km with the
//  night profile's own logarithmic slope there, so density keeps falling
//  continuously.
//
//  INDICES AND TIME. EnvAtmUpper takes F10.7, F81 and Kp (or kpp) ALREADY
//  taken at the delayed times of 5.7 (F: t - 1.7 d; daily Kp: t - 0.6 d;
//  3-hour kpp: t - 0.25 d) -- UpperAtmIndicesAt does that from series;
//  and the Sun and Earth rotation of the moment (UpperAtmTimeOf from the
//  MJD, UT). The standard does not give the Sun's coordinates; here the
//  low-precision formulas of the Astronomical Almanac (about 0.01 deg)
//  and the sidereal time of the IAU 1982 expression are used.
//
//  CHECKED (tests/test_atm_upper.cpp): rho_n, K'0, K'1, K'2, K'3, K'4,
//  K''4 against tables 4-11 of the standard, every row (70 heights x 7
//  levels, 22 Kp values x 7 levels x 2 kinds).
//
//  TABLES 2-3, CORRECTIONS. Three printed values disagree with the
//  standard's own tables 6 and 7, and are corrected here:
//    c0, second range, F0 = 250: printed -147.859; with it K'1 misses
//      table 6 by 0.031 at every height 1000-1500 km; -147.828 meets it
//      within the rounding (a constant offset: the constant term);
//    d4, F0 = 200 (both tables): printed -4.27908e-13; K'2 misses table
//      7 by up to 0.015 (growing as h^4); -4.24908e-13 meets it;
//    c2, second range, F0 = 250: printed "-671937e-4" (decimal comma
//      lost), read as -6.71937e-4 -- then table 6 is met.
//  (Table 7 itself prints K'2 = 4.466 at 780 km, F0 = 125, between 2.442
//  and 2.487: a misprint of the table, the model gives 2.466.)
//
//  MODES (AtmMode of lib_atmosphere.h). Gost -- the standard as it is:
//  density steps where it switches formulas -- at 120 km (appendix A
//  against formula (1): up to about 15 % with the semiannual term), at
//  the second-range bound of each coefficient group (up to a few %), and
//  in time when F81 crosses the middle between two levels F0. Smooth --
//  for integrators: each switch is blended with the quintic weight
//  s = x^3 (10 - 15 x + 6 x^2) over +-W_120 = 5 km at 120 km (in the
//  logarithm of density), +-W_GROUP = 10 km at a group bound (the
//  polynomial of the group; for a -- the exponent), +-W_F = 2.5 units of
//  F81 at the middle between two levels (in the logarithm of density);
//  C2 in height and in F81, Gost to the bit outside the windows.
//
//  C++17, no STL, no heap, no exceptions.
//======================================================================

#ifndef LIB_ATM_UPPER_H
#define LIB_ATM_UPPER_H

#include <math.h>

#include "lib_atmosphere.h"
#include "lib_frames.h"
#include "lib_geodesy.h"

namespace gnc
{

//  Solar and geomagnetic activity for formula (1), taken at the delayed
//  times (5.7): F10.7 and F81 at t - 1.7 d, Kp at t - 0.6 d (daily) or
//  kpp at t - 0.25 d (3-hour).
struct UpperAtmIndices
{
    double f107;    // daily F10.7, 1e-22 W/(m^2 Hz)
    double f81;     // weighted 81-day F81 (3.2), 1e-22 W/(m^2 Hz)
    double kp;      // daily Kp, or the modified 3-hour kpp when kp3h
    bool kp3h;      // true: kp is the 3-hour kpp (coefficients et5 ... et8)
};

//  The Sun and the Earth's rotation at the moment (UpperAtmTimeOf).
struct UpperAtmTime
{
    double ra;      // right ascension of the Sun, rad
    double dec;     // declination of the Sun, rad
    double sStar;   // Greenwich sidereal time at 0h UT of the date, rad
    double tUt;     // universal time since 0h UT of the date, s
    double day;     // d: days from the start of the year (1 January 0h UT = 1)
};

//  The parts of formula (1) at one point (EnvAtmUpper::Parts); above
//  120 km only (below, rho = rhoN = the appendix A value, K = 0).
struct UpperAtmParts
{
    double rhoN;    // night density rho_n, kg/m^3
    double k0;      // K0
    double k1;      // K1
    double k2;      // K2
    double k3;      // K3
    double k4;      // K4
    int level;      // index of F0: 0 ... 6 = 75, 100, 125, 150, 175, 200, 250
};

class EnvAtmUpper
{
public:

    static constexpr int N_LEV = 7;
    static constexpr double H_LOW = 120.0;          // formula (1) from, km
    static constexpr double H_TOP = 1500.0;         // formula (1) to, km
    static constexpr double RHO_0 = 1.58868e-8;     // night density at 120 km, kg/m^3
    static constexpr double OMEGA_3 = 7.292115e-5;  // Earth rotation rate of the standard, rad/s
    static constexpr double DT_F = 1.7;             // delay of the solar indices, days (5.7)
    static constexpr double DT_KP_DAY = 0.6;        // delay of daily Kp, days
    static constexpr double DT_KP_3H = 0.25;        // delay of 3-hour kpp, days
    static constexpr double W_120 = 5.0;            // Smooth: half-width at 120 km, km
    static constexpr double W_GROUP = 10.0;         // Smooth: half-width at a group bound, km
    static constexpr double W_F = 2.5;              // Smooth: half-width in F81 at a level middle

    constexpr explicit EnvAtmUpper(AtmMode mode = AtmMode::Gost)
        : m_mode(mode)
    {
    }

    constexpr void SetMode(AtmMode mode)
    {
        m_mode = mode;
    }

    constexpr AtmMode Mode() const
    {
        return m_mode;
    }

    //  Density, kg/m^3, at a point.
    //    p  -- the point, ECEF (Greenwich) axes, m
    //    h  -- its height above the ellipsoid, m (the standard: the general
    //          terrestrial ellipsoid -- PZ-90; Density(p, tm, ix) finds it)
    //    tm -- the Sun and the Earth rotation of the moment
    //    ix -- indices at the delayed times
    double Density(const CfEcefPoint& p, double h, const UpperAtmTime& tm, const UpperAtmIndices& ix) const
    {
        double hk = h / 1000.0;
        if (m_mode == AtmMode::Smooth)
        {
            double d = hk - H_LOW;
            if (d > -W_120 && d < W_120)
            {
                double s = Weight(d, W_120);
                double lo = log(DensityLowKm(hk));
                double hi = log(DensityF(p, hk, tm, ix));
                return exp(lo + s * (hi - lo));
            }
        }
        if (hk < H_LOW)
        {
            return DensityLowKm(hk);
        }
        return DensityF(p, hk, tm, ix);
    }

    //  The same, the height found above the PZ-90.11 ellipsoid.
    double Density(const CfEcefPoint& p, const UpperAtmTime& tm, const UpperAtmIndices& ix) const
    {
        return Density(p, EcpToAlt(p, ELL_PZ9011), tm, ix);
    }

    //  Appendix A.4: density below 120 km, kg/m^3 (no variations).
    //    h -- height above the ellipsoid, m
    static double DensityLow(double h)
    {
        return DensityLowKm(h / 1000.0);
    }

    //  The parts of formula (1) at one point (Gost; level by F81).
    //  Heights below 120 km give the formula as well (extrapolated).
    UpperAtmParts Parts(const CfEcefPoint& p, double h, const UpperAtmTime& tm, const UpperAtmIndices& ix) const
    {
        UpperAtmParts q { };
        q.level = LevelOf(ix.f81);
        Eval(q.level, p, h / 1000.0, tm, ix, false, &q);
        return q;
    }

    //  Index of the level F0 nearest to F81 (5.5); a tie goes to the lower.
    static int LevelOf(double f81)
    {
        int best = 0;
        int i = 0;
        for (i = 1; i < N_LEV; i = i + 1)
        {
            if (fabs(f81 - F0(i)) < fabs(f81 - F0(best)))
            {
                best = i;
            }
        }
        return best;
    }

    //  F0 of a level, 1e-22 W/(m^2 Hz).
    static double F0(int i)
    {
        static const double F[N_LEV] = { 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 250.0 };
        return F[i];
    }

    //  The height polynomials of tables 2-3 by name, as tables 4-9 give
    //  them: rho_n (kg/m^3), K'0 (l), K'1 (c), K'2 (d), K'3 (b), K'4 (e);
    //  and K''4 of tables 10-11. Gost (no blending).
    //    level -- 0 ... 6; hKm -- height, km
    static double NightDensity(int level, double hKm)
    {
        return RHO_0 * exp(Poly(GA, level, hKm, false));
    }
    static double KP0(int level, double hKm)
    {
        return Poly(GL, level, hKm, false);
    }
    static double KP1(int level, double hKm)
    {
        return Poly(GC, level, hKm, false);
    }
    static double KP2(int level, double hKm)
    {
        return Poly(GD, level, hKm, false);
    }
    static double KP3(int level, double hKm)
    {
        return Poly(GB, level, hKm, false);
    }
    static double KP4(int level, double hKm)
    {
        return Poly(GE, level, hKm, false);
    }
    static double KPP4(int level, double kp, bool kp3h)
    {
        const double* e = Ekp(level, kp3h);
        return e[0] + kp * (e[1] + kp * (e[2] + kp * e[3]));
    }

    //  A(d), table 1: the semiannual effect.
    //    d -- days from the start of the year
    static double SemiAnnual(double d)
    {
        static const double A[9] =
        {
            -2.53418e-2, -2.44075e-3, 3.08389e-6, 2.90115e-6, -4.99606e-8,
            3.36327e-10, -1.0966e-12, 1.73227e-15, -1.06271e-18
        };
        double s = A[8];
        int k = 0;
        for (k = 7; k >= 0; k = k - 1)
        {
            s = s * d + A[k];
        }
        return s;
    }

private:

    //  Coefficient groups.
    static constexpr int GA = 0;
    static constexpr int GB = 1;
    static constexpr int GC = 2;
    static constexpr int GD = 3;
    static constexpr int GE = 4;
    static constexpr int GL = 5;

    //  Quintic blend weight at d in (-w, w): 0 at -w, 1 at +w, s' = s'' = 0 there.
    static double Weight(double d, double w)
    {
        double x = (d + w) / (2.0 * w);
        return x * x * x * (10.0 - 15.0 * x + 6.0 * x * x);
    }

    //  Lower bound of the second height range of a group, km (table 3).
    static double Bound(int g, int level)
    {
        static const double B[6][N_LEV] =
        {
            { 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0 },        // a
            { 600.0, 660.0, 760.0, 800.0, 860.0, 900.0, 1000.0 },       // b
            { 640.0, 700.0, 760.0, 820.0, 860.0, 920.0, 980.0 },        // c
            { 1.0e9, 1.0e9, 1.0e9, 1.0e9, 1.0e9, 1.0e9, 1.0e9 },        // d: one range (table 3 = table 2)
            { 600.0, 700.0, 780.0, 800.0, 800.0, 900.0, 760.0 },        // e
            { 640.0, 660.0, 740.0, 800.0, 860.0, 900.0, 900.0 }         // l
        };
        return B[g][level];
    }

    //  Coefficients of a group: range 0/1, level; c[0] ... c[n-1].
    static const double* Coef(int g, int range, int level, int* n)
    {
        //  a0 ... a6: [range][level][power]
        static const double A[2][N_LEV][7] =
        {
            {
                { 26.8629, -0.451674, 0.00290397, -1.06953e-05, 2.21598e-08, -2.42941e-11, 1.09926e-14 },   // F0 = 75
                { 27.4598, -0.463668, 0.002974, -1.0753e-05, 2.17059e-08, -2.30249e-11, 1.00123e-14 },   // F0 = 100
                { 28.6395, -0.490987, 0.00320649, -1.1681e-05, 2.36847e-08, -2.51809e-11, 1.09536e-14 },   // F0 = 125
                { 29.6418, -0.514957, 0.00341926, -1.25785e-05, 2.5727e-08, -2.75874e-11, 1.21091e-14 },   // F0 = 150
                { 30.1671, -0.527837, 0.00353211, -1.30227e-05, 2.66455e-08, -2.85432e-11, 1.25009e-14 },   // F0 = 175
                { 29.7578, -0.517915, 0.00342699, -1.24137e-05, 2.48209e-08, -2.58413e-11, 1.09383e-14 },   // F0 = 200
                { 30.7854, -0.545695, 0.00370328, -1.37072e-05, 2.80614e-08, -3.00184e-11, 1.31142e-14 },   // F0 = 250
            },
            {
                { 17.8781, -0.132025, 0.000227717, -2.2543e-07, 1.33574e-10, -4.50458e-14, 6.72086e-18 },   // F0 = 75
                { -2.54909, 0.0140064, -0.00016946, 3.27196e-07, -2.8763e-10, 1.22625e-13, -2.05736e-17 },   // F0 = 100
                { -13.9599, 0.0844951, -0.000328875, 5.05918e-07, -3.92299e-10, 1.52279e-13, -2.35576e-17 },   // F0 = 125
                { -23.3079, 0.135141, -0.000420802, 5.73717e-07, -4.03238e-10, 1.42846e-13, -2.01726e-17 },   // F0 = 150
                { -14.7264, 0.0713256, -0.000228015, 2.8487e-07, -1.74383e-10, 5.08071e-14, -5.34955e-18 },   // F0 = 175
                { -4.912, 0.0108326, -8.10546e-05, 1.15712e-07, -8.13296e-11, 3.04913e-14, -4.94989e-18 },   // F0 = 200
                { -5.40952, 0.00550749, -3.78851e-05, 2.4808e-08, 4.92183e-12, -8.65011e-15, 1.9849e-18 },   // F0 = 250
            }
        };
        //  b0 ... b4: [range][level][power]
        static const double B[2][N_LEV][5] =
        {
            {
                { 0.0687894, -0.00284077, 1.83922e-05, 9.19605e-09, -4.16873e-11 },   // F0 = 75
                { 0.15073, -0.00400889, 2.43937e-05, -9.92772e-09, -1.82239e-11 },   // F0 = 100
                { 0.0479451, -0.00239453, 1.70335e-05, -1.31626e-09, -1.74032e-11 },   // F0 = 125
                { 0.0223448, -0.0019798, 1.54101e-05, -2.3543e-09, -1.24994e-11 },   // F0 = 150
                { -0.00326391, -0.00159869, 1.40443e-05, -3.02287e-09, -9.2016e-12 },   // F0 = 175
                { -0.0514749, -0.000921059, 1.15147e-05, -1.22901e-09, -8.13104e-12 },   // F0 = 200
                { -0.107255, -0.000174343, 9.02759e-06, -3.16512e-10, -6.14e-12 },   // F0 = 250
            },
            {
                { 23.1584, -0.0802147, 0.000105824, -6.15036e-08, 1.32453e-11 },   // F0 = 75
                { 33.2732, -0.111099, 0.000141421, -7.94952e-08, 1.65836e-11 },   // F0 = 100
                { 39.1961, -0.12352, 0.000149015, -7.9705e-08, 1.58772e-11 },   // F0 = 125
                { 43.2469, -0.126973, 0.000142637, -7.09985e-08, 1.31646e-11 },   // F0 = 150
                { 49.5738, -0.138613, 0.000147851, -6.96361e-08, 1.21595e-11 },   // F0 = 175
                { 11.278, 0.00143478, -3.69846e-05, 3.58318e-08, -9.91225e-12 },   // F0 = 200
                { -52.6184, 0.214689, -0.000294882, 1.71171e-07, -3.60582e-11 },   // F0 = 250
            }
        };
        //  c0 ... c4: [range][level][power]
        static const double C[2][N_LEV][5] =
        {
            {
                { -1.04825, 0.0166305, -9.24263e-05, 2.72382e-07, -2.41355e-10 },   // F0 = 75
                { -0.93106, 0.0141537, -7.29862e-05, 2.00294e-07, -1.62006e-10 },   // F0 = 100
                { -0.820867, 0.0119916, -5.79835e-05, 1.50707e-07, -1.13026e-10 },   // F0 = 125
                { -0.744047, 0.0104743, -4.78544e-05, 1.18513e-07, -8.31498e-11 },   // F0 = 150
                { -0.722471, 0.00980317, -4.25245e-05, 9.95544e-08, -6.55175e-11 },   // F0 = 175
                { -0.687482, 0.00916594, -3.80932e-05, 8.51275e-08, -5.29972e-11 },   // F0 = 200
                { -0.739984, 0.00952854, -3.62727e-05, 7.3887e-08, -4.23907e-11 },   // F0 = 250
            },
            {
                { 50.5034, -0.170541, 0.000217232, -1.21902e-07, 2.54037e-11 },   // F0 = 75
                { 61.624, -0.192967, 0.000228061, -1.18715e-07, 2.29638e-11 },   // F0 = 100
                { 53.2623, -0.144342, 0.00014659, -6.46443e-08, 1.04227e-11 },   // F0 = 125
                { 18.2236, -0.00840024, -3.88e-05, 4.31384e-08, -1.23832e-11 },   // F0 = 150
                { -31.8442, 0.168327, -0.000262603, 1.65454e-07, -3.69355e-11 },   // F0 = 175
                { -48.7208, 0.222996, -0.000321884, 1.91495e-07, -4.08067e-11 },   // F0 = 200
                { -147.828, 0.531652, -0.000671937, 3.64787e-07, -7.26268e-11 },   // F0 = 250
            }
        };
        //  d0 ... d4: [range][level][power]
        static const double D[2][N_LEV][5] =
        {
            {
                { -0.351899, 0.00577056, 9.95819e-07, -7.25324e-09, 2.9759e-12 },   // F0 = 75
                { -0.047813, 0.00380813, 4.22771e-06, -8.66826e-09, 3.06712e-12 },   // F0 = 100
                { 0.20981, 0.00262881, 4.24379e-06, -6.67328e-09, 2.13496e-12 },   // F0 = 125
                { 0.265174, 0.00275836, 2.08668e-06, -3.69543e-09, 1.11862e-12 },   // F0 = 150
                { 0.23047, 0.00338331, -5.52305e-07, -8.23607e-10, 2.21349e-13 },   // F0 = 175
                { 0.170074, 0.00406131, -2.82114e-06, 1.38369e-09, -4.24908e-13 },   // F0 = 200
                { 0.088141, 0.00468253, -4.24609e-06, 2.53509e-09, -7.29031e-13 },   // F0 = 250
            },
            {
                { -0.351899, 0.00577056, 9.95819e-07, -7.25324e-09, 2.9759e-12 },   // F0 = 75
                { -0.047813, 0.00380813, 4.22771e-06, -8.66826e-09, 3.06712e-12 },   // F0 = 100
                { 0.20981, 0.00262881, 4.24379e-06, -6.67328e-09, 2.13496e-12 },   // F0 = 125
                { 0.265174, 0.00275836, 2.08668e-06, -3.69543e-09, 1.11862e-12 },   // F0 = 150
                { 0.23047, 0.00338331, -5.52305e-07, -8.23607e-10, 2.21349e-13 },   // F0 = 175
                { 0.170074, 0.00406131, -2.82114e-06, 1.38369e-09, -4.24908e-13 },   // F0 = 200
                { 0.088141, 0.00468253, -4.24609e-06, 2.53509e-09, -7.29031e-13 },   // F0 = 250
            }
        };
        //  e0 ... e4: [range][level][power]
        static const double E[2][N_LEV][5] =
        {
            {
                { -0.731596, 0.00597345, -5.82037e-06, 6.84634e-08, -9.50483e-11 },   // F0 = 75
                { -0.752175, 0.00565925, 1.8082e-06, 3.33822e-08, -5.13965e-11 },   // F0 = 100
                { -0.570476, 0.00295802, 1.68896e-05, -4.7475e-09, -1.72711e-11 },   // F0 = 125
                { -0.949573, 0.00813121, -3.87813e-06, 2.37694e-08, -2.77469e-11 },   // F0 = 150
                { -0.967598, 0.00841991, -3.585e-06, 1.74801e-08, -1.96221e-11 },   // F0 = 175
                { -1.02278, 0.00923633, -6.10128e-06, 1.78211e-08, -1.70073e-11 },   // F0 = 200
                { -0.757903, 0.00606068, 7.85296e-06, -9.74891e-09, 1.58377e-12 },   // F0 = 250
            },
            {
                { 38.6199, -0.132147, 0.000175411, -1.02417e-07, 2.21446e-11 },   // F0 = 75
                { 51.249, -0.167373, 0.000211832, -1.18221e-07, 2.45055e-11 },   // F0 = 100
                { 68.4746, -0.215659, 0.000262273, -1.40972e-07, 2.82285e-11 },   // F0 = 125
                { 58.422, -0.166664, 0.000185486, -9.12345e-08, 1.67118e-11 },   // F0 = 150
                { 7.20188, 0.0216109, -6.52882e-05, 5.37077e-08, -1.4095e-11 },   // F0 = 175
                { 21.5948, -0.0202239, -1.72029e-05, 2.83017e-08, -8.94486e-12 },   // F0 = 200
                { -88.4076, 0.338518, -0.000445581, 2.51729e-07, -5.203e-11 },   // F0 = 250
            }
        };
        //  l0 ... l4: [range][level][power]
        static const double L[2][N_LEV][5] =
        {
            {
                { -0.407768, 0.00148506, 1.25357e-05, 3.77311e-08, -7.78953e-11 },   // F0 = 75
                { -0.902739, 0.00826803, -1.25448e-05, 6.12853e-08, -7.07966e-11 },   // F0 = 100
                { -0.733037, 0.00523396, 6.35667e-06, 1.09065e-08, -2.61427e-11 },   // F0 = 125
                { -1.31444, 0.0133124, -2.55585e-05, 5.43981e-08, -4.33784e-11 },   // F0 = 150
                { -1.20026, 0.0114087, -1.47324e-05, 2.7804e-08, -2.2632e-11 },   // F0 = 175
                { -1.52158, 0.015704, -3.02859e-05, 4.57668e-08, -2.82926e-11 },   // F0 = 200
                { -1.67664, 0.0177194, -3.69498e-05, 5.09134e-08, -2.82878e-11 },   // F0 = 250
            },
            {
                { 48.6536, -0.170291, 0.000226242, -1.32032e-07, 2.85193e-11 },   // F0 = 75
                { 54.4867, -0.178298, 0.000222725, -1.227e-07, 2.51316e-11 },   // F0 = 100
                { 60.1267, -0.183144, 0.000212481, -1.08497e-07, 2.0571e-11 },   // F0 = 125
                { 47.0996, -0.12526, 0.000126352, -5.51584e-08, 8.75272e-12 },   // F0 = 150
                { 50.6174, -0.129047, 0.000124842, -5.24993e-08, 8.08272e-12 },   // F0 = 175
                { 8.01942, 0.0185302, -6.14733e-05, 4.97674e-08, -1.26162e-11 },   // F0 = 200
                { -15.5728, 0.0936704, -0.000149036, 9.42151e-08, -2.0961e-11 },   // F0 = 250
            }
        };
        if (g == GA)
        {
            *n = 7;
            return A[range][level];
        }
        *n = 5;
        if (g == GB)
        {
            return B[range][level];
        }
        if (g == GC)
        {
            return C[range][level];
        }
        if (g == GD)
        {
            return D[range][level];
        }
        if (g == GE)
        {
            return E[range][level];
        }
        return L[range][level];
    }

    //  Polynomial of one range at h, km.
    static double PolyRange(int g, int range, int level, double h)
    {
        int n = 0;
        const double* c = Coef(g, range, level, &n);
        double s = c[n - 1];
        int k = 0;
        for (k = n - 2; k >= 0; k = k - 1)
        {
            s = s * h + c[k];
        }
        return s;
    }

    //  Polynomial of a group at h, km: first range up to the bound
    //  (inclusive), second above; smooth -- blended across the bound.
    static double Poly(int g, int level, double h, bool smooth)
    {
        double b = Bound(g, level);
        if (smooth)
        {
            double d = h - b;
            if (d > -W_GROUP && d < W_GROUP)
            {
                double s = Weight(d, W_GROUP);
                double p0 = PolyRange(g, 0, level, h);
                double p1 = PolyRange(g, 1, level, h);
                return p0 + s * (p1 - p0);
            }
        }
        if (h > b)
        {
            return PolyRange(g, 1, level, h);
        }
        return PolyRange(g, 0, level, h);
    }

    //  e5 ... e8 (daily Kp) or et5 ... et8 (3-hour kpp) of a level.
    static const double* Ekp(int level, bool kp3h)
    {
        static const double E[2][N_LEV][4] =
        {
            {
                { -0.20670, 9.7533e-2, -1.1817e-2, 1.6145e-3 },
                { -0.16971, 7.9830e-2, -9.4393e-3, 1.2622e-3 },
                { -0.14671, 6.8808e-2, -7.9836e-3, 1.0535e-3 },
                { -0.13150, 6.1603e-2, -7.0866e-3, 9.2813e-4 },
                { -0.120916, 5.6538e-2, -6.4324e-3, 8.3723e-4 },
                { -0.11363, 5.3178e-2, -6.0436e-3, 7.7982e-4 },
                { -0.10444, 4.8551e-2, -5.3567e-3, 6.8809e-4 }
            },
            {
                { -0.2061, 9.4449e-2, -8.7953e-3, 8.8385e-4 },
                { -0.169279, 7.7599e-2, -7.1375e-3, 6.9025e-4 },
                { -0.146377, 6.7052e-2, -6.0951e-3, 5.7456e-4 },
                { -0.13121, 6.0105e-2, -5.4388e-3, 5.0585e-4 },
                { -0.12067, 5.5232e-2, -4.9580e-3, 4.5512e-4 },
                { -0.113399, 5.1994e-2, -4.6876e-3, 4.2548e-4 },
                { -0.104243, 4.7573e-2, -4.1711e-3, 3.7068e-4 }
            }
        };
        return E[kp3h ? 1 : 0][level];
    }

    //  phi1 of a level, rad (the lag of the density maximum).
    static double Phi1(int level)
    {
        static const double P[N_LEV] = { 0.5411, 0.5515, 0.5585, 0.5585, 0.5585, 0.5585, 0.5585 };
        return P[level];
    }

    //  Formula (1) at one level, h in km (any h; above H_TOP see Density).
    //  Fills the parts when q is given.
    double Eval(int level, const CfEcefPoint& p, double h, const UpperAtmTime& tm,
                const UpperAtmIndices& ix, bool smooth, UpperAtmParts* q) const
    {
        double f0 = F0(level);
        double rhoN = RHO_0 * exp(Poly(GA, level, h, smooth));
        double k0 = 1.0 + Poly(GL, level, h, smooth) * (ix.f81 - f0) / f0;
        //  Diurnal effect: the angle from the density maximum.
        double beta = tm.ra - tm.sStar - OMEGA_3 * tm.tUt + Phi1(level);
        double r = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        double cosPhi = 0.0;
        if (r > 0.0)
        {
            cosPhi = (p.z * sin(tm.dec) + cos(tm.dec) * (p.x * cos(beta) + p.y * sin(beta))) / r;
        }
        double c2 = 0.5 * (1.0 + cosPhi);
        if (c2 < 0.0)
        {
            c2 = 0.0;
        }
        double cosHalf = sqrt(c2);
        double nExp = 2.058 + h * (5.887e-3 - 4.012e-6 * h);
        double k1 = 0.0;
        if (cosHalf > 0.0)
        {
            k1 = Poly(GC, level, h, smooth) * pow(cosHalf, nExp);
        }
        double k2 = Poly(GD, level, h, smooth) * SemiAnnual(tm.day);
        double df = ix.f107 - ix.f81;
        double k3 = 0.0;
        double den = ix.f81 + fabs(df);
        if (den > 0.0)
        {
            k3 = Poly(GB, level, h, smooth) * df / den;
        }
        double k4 = Poly(GE, level, h, smooth) * KPP4(level, ix.kp, ix.kp3h);
        if (q != nullptr)
        {
            q->rhoN = rhoN;
            q->k0 = k0;
            q->k1 = k1;
            q->k2 = k2;
            q->k3 = k3;
            q->k4 = k4;
        }
        return rhoN * k0 * (1.0 + k1 + k2 + k3 + k4);
    }

    //  Formula (1) with the level chosen (or blended) by F81 and the
    //  continuation above H_TOP. h in km.
    double DensityF(const CfEcefPoint& p, double h, const UpperAtmTime& tm, const UpperAtmIndices& ix) const
    {
        bool smooth = m_mode == AtmMode::Smooth;
        int lev = LevelOf(ix.f81);
        if (smooth)
        {
            int k = 0;
            for (k = 0; k < N_LEV - 1; k = k + 1)
            {
                double mid = 0.5 * (F0(k) + F0(k + 1));
                double d = ix.f81 - mid;
                if (d > -W_F && d < W_F)
                {
                    double s = Weight(d, W_F);
                    double lo = log(DensityLevel(k, p, h, tm, ix, true));
                    double hi = log(DensityLevel(k + 1, p, h, tm, ix, true));
                    return exp(lo + s * (hi - lo));
                }
            }
        }
        return DensityLevel(lev, p, h, tm, ix, smooth);
    }

    //  One level with the continuation above H_TOP: the value at H_TOP
    //  (point moved radially), times exp(slope (h - H_TOP)), slope -- the
    //  logarithmic derivative of rho_n at H_TOP.
    double DensityLevel(int level, const CfEcefPoint& p, double h, const UpperAtmTime& tm,
                        const UpperAtmIndices& ix, bool smooth) const
    {
        if (h <= H_TOP)
        {
            return Eval(level, p, h, tm, ix, smooth, nullptr);
        }
        double r = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        double k = (r - (h - H_TOP) * 1000.0) / r;
        CfEcefPoint pt { };
        pt.x = p.x * k;
        pt.y = p.y * k;
        pt.z = p.z * k;
        double top = Eval(level, pt, H_TOP, tm, ix, smooth, nullptr);
        int n = 0;
        const double* a = Coef(GA, 1, level, &n);
        double slope = 0.0;
        int j = 0;
        for (j = n - 1; j >= 1; j = j - 1)
        {
            slope = slope * H_TOP + (double)j * a[j];
        }
        return top * exp(slope * (h - H_TOP));
    }

    //  Appendix A.4, table A.2; h in km (below 0 -- the first layer).
    static double DensityLowKm(double h)
    {
        static const double HB[4] = { 0.0, 20.0, 60.0, 100.0 };
        static const double A0[4] = { 1.228, 9.013e-2, 3.104e-4, 3.66e-7 };
        static const double K1[4] = { -9.0764e-2, -0.16739, -0.137, -0.18553 };
        static const double K2[4] = { -2.0452e-3, 6.2669e-4, -7.8653e-4, 1.5397e-3 };
        int i = 0;
        if (h >= 20.0)
        {
            i = 1;
        }
        if (h >= 60.0)
        {
            i = 2;
        }
        if (h >= 100.0)
        {
            i = 3;
        }
        double d = h - HB[i];
        return A0[i] * exp(K1[i] * d + K2[i] * d * d);
    }

    AtmMode m_mode;
};

//------------------------------------------------------------------------
//  Time and indices.
//------------------------------------------------------------------------

//  Days from 1 March of year 0 (proleptic Gregorian) to a civil date;
//  internal (H. Hinnant's days_from_civil, shifted).
inline long UpperAtmDaysFromCivil(int y, int m, int d)
{
    long yy = (long)y - (m <= 2 ? 1 : 0);
    long era = (yy >= 0 ? yy : yy - 399) / 400;
    long yoe = yy - era * 400;
    long mp = (long)((m + 9) % 12);
    long doy = (153 * mp + 2) / 5 + (long)d - 1;
    long doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return era * 146097 + doe;
}

//  MJD (0h UT) of a civil date.
inline double MjdOfDate(int year, int month, int day)
{
    //  days_from_civil gives 0 at 0000-03-01; MJD 0 = 1858-11-17.
    return (double)(UpperAtmDaysFromCivil(year, month, day) - UpperAtmDaysFromCivil(1858, 11, 17));
}

//  Civil year of an MJD (UT).
inline int YearOfMjd(double mjd)
{
    long z = (long)floor(mjd) + UpperAtmDaysFromCivil(1858, 11, 17);
    long era = (z >= 0 ? z : z - 146096) / 146097;
    long doe = z - era * 146097;
    long yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    long y = yoe + era * 400;
    long doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    long mp = (5 * doy + 2) / 153;
    long m = mp < 10 ? mp + 3 : mp - 9;
    return (int)(m <= 2 ? y + 1 : y);
}

//  The Sun and the Earth rotation for formula (1) at a moment.
//    mjd -- modified Julian date, UT (UTC to within a second is enough)
//  Sun: the low-precision formulas of the Astronomical Almanac (about
//  0.01 deg, 1950-2050); sidereal time: GMST = 6.697374558 h +
//  0.06570982441908 h D0 at 0h UT, D0 -- days from J2000.0.
inline UpperAtmTime UpperAtmTimeOf(double mjd)
{
    const double deg = PI / 180.0;
    UpperAtmTime tm { };
    double n = mjd - 51544.5;
    double l = 280.460 + 0.9856474 * n;
    double g = (357.528 + 0.9856003 * n) * deg;
    double lam = (l + 1.915 * sin(g) + 0.020 * sin(2.0 * g)) * deg;
    double eps = (23.439 - 0.0000004 * n) * deg;
    tm.ra = atan2(cos(eps) * sin(lam), cos(lam));
    if (tm.ra < 0.0)
    {
        tm.ra = tm.ra + 2.0 * PI;
    }
    tm.dec = asin(sin(eps) * sin(lam));
    double mjd0 = floor(mjd);
    double gmst0 = 6.697374558 + 0.06570982441908 * (mjd0 - 51544.5);
    gmst0 = fmod(gmst0, 24.0);
    if (gmst0 < 0.0)
    {
        gmst0 = gmst0 + 24.0;
    }
    tm.sStar = gmst0 * (PI / 12.0);
    tm.tUt = (mjd - mjd0) * 86400.0;
    int y = YearOfMjd(mjd);
    tm.day = 1.0 + (mjd - MjdOfDate(y, 1, 1));
    return tm;
}

//  F81 (3.2) of day i of a daily F10.7 series: weights W = 1 + 0.5 j / 80,
//  j = -80 ... 0. Needs i >= 80.
inline double F81Of(const double* f107, int i)
{
    double sw = 0.0;
    double sf = 0.0;
    int j = 0;
    for (j = -80; j <= 0; j = j + 1)
    {
        double w = 1.0 + 0.5 * (double)j / 80.0;
        sw = sw + w;
        sf = sf + w * f107[i + j];
    }
    return sf / sw;
}

//  Modified 3-hour indices kpp (section 4): kpp_0 = kp_0,
//  kpp_j = kp_j - r (kp_j - kpp_{j-1}), r = 0.3 when the index rises,
//  0.7 when it falls.
inline void KppOf(const double* kp, int n, double* kpp)
{
    int j = 0;
    for (j = 0; j < n; j = j + 1)
    {
        if (j == 0)
        {
            kpp[j] = kp[j];
            continue;
        }
        double d = kp[j] - kpp[j - 1];
        double r = 0.0;
        if (d > 0.0)
        {
            r = 0.3;
        }
        if (d < 0.0)
        {
            r = 0.7;
        }
        kpp[j] = kp[j] - r * d;
    }
}

//  Kp from Ap (or kp from ap), table A.1, linear between the rows.
inline double KpOfAp(double ap)
{
    static const double AP[28] =
    {
        0.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 9.0, 12.0, 15.0, 18.0, 22.0, 27.0, 32.0,
        39.0, 48.0, 56.0, 67.0, 80.0, 94.0, 111.0, 132.0, 154.0, 179.0, 207.0, 236.0, 300.0, 400.0
    };
    if (ap <= AP[0])
    {
        return 0.0;
    }
    if (ap >= AP[27])
    {
        return 9.0;
    }
    int i = 0;
    while (ap > AP[i + 1])
    {
        i = i + 1;
    }
    double f = (ap - AP[i]) / (AP[i + 1] - AP[i]);
    return ((double)i + f) / 3.0;
}

//  Index series for UpperAtmIndicesAt (owned by the caller).
struct UpperAtmSeries
{
    double mjd0;            // MJD of 0h UT of the first day
    int nDays;              // days in the series
    const double* f107;     // daily F10.7, nDays values
    const double* kpDay;    // daily Kp, nDays values (referred to 12:00 UT), or nullptr
    const double* kpp3h;    // modified 3-hour kpp (KppOf), 8 nDays values
                            // (referred to the middle of each interval), or nullptr
};

//  Value of a series sampled at t0 + k step, linear between samples;
//  false outside.
inline bool UpperAtmSample(const double* v, int n, double t0, double step, double t, double* out)
{
    double x = (t - t0) / step;
    if (x < 0.0 || x > (double)(n - 1))
    {
        return false;
    }
    int k = (int)floor(x);
    if (k > n - 2)
    {
        k = n - 2;
    }
    double f = x - (double)k;
    *out = v[k] + f * (v[k + 1] - v[k]);
    return true;
}

//  Indices for formula (1) at a moment, appendix A.2: each at its delayed
//  time, linear between its reference times -- F10.7 at 20:00 UT (17:00
//  UT before 1991-06-01), daily Kp at 12:00 UT, 3-hour kpp at the middle
//  of its interval. Uses kpp3h when given, else kpDay.
//    mjd -- the moment, MJD UT
//  Returns false when the series do not cover the moment (F81 needs 80
//  days of F10.7 before the delayed time).
inline bool UpperAtmIndicesAt(const UpperAtmSeries& s, double mjd, UpperAtmIndices* ix)
{
    if (s.f107 == nullptr || s.nDays < 82)
    {
        return false;
    }
    double tF = mjd - EnvAtmUpper::DT_F;
    double hourF = tF < 48408.0 ? 17.0 : 20.0;     // 48408 = 1991-06-01
    double x = tF - (s.mjd0 + hourF / 24.0);
    int k = (int)floor(x);
    if (k < 80 || k + 1 > s.nDays - 1)
    {
        return false;
    }
    double f = x - (double)k;
    ix->f107 = s.f107[k] + f * (s.f107[k + 1] - s.f107[k]);
    double a = F81Of(s.f107, k);
    double b = F81Of(s.f107, k + 1);
    ix->f81 = a + f * (b - a);
    if (s.kpp3h != nullptr)
    {
        ix->kp3h = true;
        return UpperAtmSample(s.kpp3h, 8 * s.nDays, s.mjd0 + 1.5 / 24.0, 0.125,
                              mjd - EnvAtmUpper::DT_KP_3H, &ix->kp);
    }
    if (s.kpDay != nullptr)
    {
        ix->kp3h = false;
        return UpperAtmSample(s.kpDay, s.nDays, s.mjd0 + 0.5, 1.0, mjd - EnvAtmUpper::DT_KP_DAY, &ix->kp);
    }
    return false;
}

}   // namespace gnc

#endif  // LIB_ATM_UPPER_H

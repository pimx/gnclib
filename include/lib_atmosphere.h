//======================================================================
//  lib_atmosphere.h
//  EnvAtmosphere - standard atmosphere GOST 4401-81 from -2 km to 250 km,
//  with a temperature offset from the standard (hot/cold day: temperature
//  shifted, pressure profile kept, density recomputed from the state
//  equation).
//
//  PROFILE (1.8.0) -- GOST 4401-81 throughout, its formulas and constants:
//    -2 ... 120 km -- molar temperature T_M linear in the GEOPOTENTIAL
//      height H (layers below), pressure by the hydrostatic equation with
//      T_M: p = p_b (T_M / T_Mb)^(-g0 / (R beta)), or the exponent in an
//      isothermal layer. Above 94 km (geometric) the molar mass M falls
//      (MolarMass below); kinetic temperature T = T_M M / M_C, density
//      rho = p / (R T_M) = p M / (R* T).
//    120 ... 250 km -- kinetic temperature linear in the GEOMETRIC height
//      (GOST table 6), molar mass a cubic in h (table 3), particle
//      concentration a 4th-degree polynomial in h (table 7, eq. 9);
//      p = n k T, rho = n M / N_A (the same as p M / (R* T)).
//    above 250 km -- no air: pressure and density 0, temperature held.
//
//  THE STANDARD'S OWN STEPS. GOST passes from the barometric formula to
//  the concentration polynomial at 120 km, and from one polynomial to the
//  next at 150 and 200 km, without joining them: pressure and density step
//  down by 0.75 % at 120 km (temperature up by 0.003 K), by 0.08 % at
//  150 km and by 0.75 % at 200 km. GOST table 4 gives each knot the value
//  of the interval below it. Mode Gost (the default) reproduces the
//  standard and keeps the steps; a knot belongs to the interval above it.
//
//  SMOOTH MODE (1.9.0, AtmMode::Smooth, ATM_STD_SMOOTH) -- for integrators.
//  Within +-W_SMOOTH = 2 km of each knot both pieces are evaluated (each
//  a little past its interval) and blended with the quintic weight
//  s = x^3 (10 - 15 x + 6 x^2): temperature and molar mass linearly,
//  pressure in its logarithm, density from the state equation. The
//  profile is then C2 in height across the knots (the standard's own
//  temperature kinks -- layer bases, 140/160/200 km -- stay, as in any
//  standard atmosphere). Outside the windows the numbers are Gost to the
//  bit; inside, pressure and density move from one piece to the other
//  and differ from Gost by up to 0.53 % (from GOST table 4 rows, up to
//  0.42 %, at 119 km) -- less than the 0.75 % the standard itself steps.
//  MEASURED (tests/test_atmosphere.cpp): RK4, a fall with drag from 210
//  to 110 km, height error at dt 0.9 / 0.45 / 0.225 s -- Gost 1.0e-2 /
//  7.3e-3 / 3.2e-3 m (order 1, set by the steps), Smooth 6.1e-3 /
//  1.3e-4 / 3.9e-7 m.
//
//  94 ... 97 km MOLAR MASS. The square-root formula of GOST section 4 is
//  not legible in the copies available to the owner; the one used
//  (MolarMassOf) meets the 3-decimal M of table 4 within 0.0009 and its
//  97 km end within 7.4e-5 kg/kmol. Below 120 km pressure and density do
//  not depend on M (they go by T_M); kinetic temperature does, within
//  0.005 K of the table.
//
//  TABLE 7, two constant terms. The text of the standard available to
//  the owner prints c0 = 0.10163937e4 (150-200 km) and 0.7631575e3
//  (200-250 km); with them the polynomial misses the standard's own
//  table 4 by up to 1.3e-4 and 1.3e-5. The constants used here,
//  1016.393796 and 763.1575226, are fitted to table 4 (every row given,
//  151-199 and 201-250 km); the misses fall to 3.5e-6 and 1.3e-6 --
//  the rounding of the 6-digit table. The other coefficients are
//  exactly as printed (120-150 km meets table 4 within 3e-6 as printed).
//
//  HEIGHT. Every function takes the GEOMETRIC height, m -- the height above
//  the WGS84 ellipsoid. The layers up to 120 km are geopotential and are
//  read through H = r h / (r + h), r = R_GEOPOT = 6 356 767 m (GOST).
//  Below -2 km geopotential the height is held at -2 km.
//
//  SOUND SPEED AND VISCOSITY. GOST tabulates them to 80 km; the same
//  formulas are used above, the sound speed with the molar mass:
//  a = sqrt(kappa R* T / M).
//
//  CHECK (tests/test_atmosphere.cpp): against GOST 4401-81 table 4, sea
//  level and 80 ... 250 km, 75 rows -- temperature within 0.0052 K,
//  pressure and density within 9.3e-5 (above 120 km within 4e-6), molar
//  mass within its printed digits; sound speed within 0.001 m/s.
//
//  History: 1.0.2 put the geometric height into the geopotential layers
//  without conversion; 1.0.3 converted it and ended at 86 km; 1.0.4 took
//  GOST to 94 km and US-76 (scaled) to 150 km; 1.8.0 -- GOST to 250 km,
//  GOST constants (R = 287.05287, r = 6 356 767 m), -2 km lower bound;
//  1.9.0 -- AtmMode::Smooth (Gost mode the same numbers to the bit).
//
//  The class is constexpr-constructible: a standard-day instance can be
//  a constant-initialized static (ATM_STD below). C++17, no STL, no heap.
//======================================================================

#ifndef LIB_ATMOSPHERE_H
#define LIB_ATMOSPHERE_H

#include <math.h>

#include "lib_frames.h"

namespace gnc
{

inline constexpr double G0_STD = 9.80665;           // standard gravity, m/s^2
inline constexpr double R_STAR = 8314.32;           // universal gas constant, J/(kmol K), GOST
inline constexpr double M_AIR0 = 28.964420;         // molar mass of air below 94 km, kg/kmol, GOST
inline constexpr double N_AVOGADRO = 602.257e24;    // Avogadro number, 1/kmol, GOST
inline constexpr double R_AIR = 287.05287;          // gas constant of air, J/(kg K), GOST
inline constexpr double R_GEOPOT = 6356767.0;       // conventional Earth radius, m, GOST
inline constexpr double K_BOLTZ = R_STAR / N_AVOGADRO;  // Boltzmann constant, J/K

//  The air at one height (1.4.0): what Temperature, Pressure, Density and
//  MolarMass return there (EnvAtmosphere::State).
struct AtmState
{
    double t;     // temperature, K (with the offset)
    double p;     // static pressure, Pa
    double rho;   // density, kg/m^3
    double m;     // molar mass, kg/kmol (1.8.0)
};

//  1.9.0: how the pieces of GOST 4401-81 meet at 120, 150 and 200 km.
//    Gost   -- as the standard: pressure and density step there (0.75 %,
//              0.08 %, 0.75 %);
//    Smooth -- the two pieces blended across +-W_SMOOTH around each knot,
//              C2 in height (see the header); exactly Gost outside.
enum class AtmMode : int
{
    Gost = 0,
    Smooth = 1
};

class EnvAtmosphere
{
public:

    static constexpr double P_SL = 101325.0;    // sea-level pressure, Pa
    static constexpr double Z_MAX = 250000.0;   // top of the model, geometric, m
    static constexpr double W_SMOOTH = 2000.0;  // Smooth: half-width of the blend at a knot, m

    //  dTemp - temperature offset from the standard, K (0 - standard day);
    //  mode  - Gost (the standard as printed) or Smooth (for integrators).
    constexpr explicit EnvAtmosphere(double dTemp = 0.0, AtmMode mode = AtmMode::Gost)
        : m_dt(dTemp)
        , m_mode(mode)
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

    constexpr void SetDeltaTemp(double dTemp)
    {
        m_dt = dTemp;
    }

    constexpr double DeltaTemp() const
    {
        return m_dt;
    }

    //  Air temperature, K, with the offset; not below 1e-3 K.
    //    h -- geometric height, m
    double Temperature(double h) const
    {
        Point q { };
        Eval(h, false, &q);
        return OffsetT(q.tIsa);
    }

    //  Static pressure, Pa; 0 above Z_MAX. Does NOT depend on dTemp
    //  (non-standard-day convention: pressure profile kept, temperature
    //  shifted).
    //    h -- geometric height, m
    double Pressure(double h) const
    {
        Point q { };
        Eval(h, true, &q);
        return q.p;
    }

    //  Density, kg/m^3, from the state equation rho = p M / (R* T) with
    //  the shifted temperature; 0 above Z_MAX.
    //    h -- geometric height, m
    double Density(double h) const
    {
        Point q { };
        Eval(h, true, &q);
        return Rho(q.p, q.m, OffsetT(q.tIsa));
    }

    //  Molar mass of air, kg/kmol: M_AIR0 up to 94 km, falling above (GOST).
    //    h -- geometric height, m
    double MolarMass(double h) const
    {
        Point q { };
        Eval(h, false, &q);
        return q.m;
    }

    //  1.4.0: temperature, pressure, density (and 1.8.0: molar mass) at
    //  once -- the same numbers as the separate calls, to the bit.
    //    h -- geometric height, m
    AtmState State(double h) const
    {
        Point q { };
        Eval(h, true, &q);
        AtmState r { };
        r.t = OffsetT(q.tIsa);
        r.p = q.p;
        r.rho = Rho(q.p, q.m, r.t);
        r.m = q.m;
        return r;
    }

    //  Pressure relative to sea level, clamped to [0, 1] (engine model).
    double PressureRatio(double h) const
    {
        double pr = Pressure(h) / P_SL;
        if (pr > 1.0)
        {
            pr = 1.0;
        }
        if (pr < 0.0)
        {
            pr = 0.0;
        }
        return pr;
    }

    //  Speed of sound, m/s: a = sqrt(kappa R* T / M), kappa = 1.4
    //  (= sqrt(1.4 R_AIR T) below 94 km, to the bit).
    double SoundSpeed(double h) const
    {
        Point q { };
        Eval(h, false, &q);
        return sqrt(1.4 * R_AIR * OffsetT(q.tIsa) * (M_AIR0 / q.m));
    }

    //  Mach number; 0 on degenerate speed of sound.
    double Mach(double h, double v) const
    {
        double a = SoundSpeed(h);
        if (a > 1.0e-6)
        {
            return v / a;
        }
        return 0.0;
    }

    //  Dynamic pressure q = 1/2 rho V^2, Pa.
    double DynamicPressure(double h, double v) const
    {
        return 0.5 * Density(h) * v * v;
    }

    //  Dynamic viscosity (Sutherland, GOST: beta_s = 1.458e-6, S = 110.4 K), Pa*s.
    double DynamicViscosity(double h) const
    {
        double t = Temperature(h);
        return 1.458e-6 * t * sqrt(t) / (t + 110.4);
    }

    //  Kinematic viscosity nu = mu / rho, m^2/s; 0 at zero density.
    double KinematicViscosity(double h) const
    {
        double rho = Density(h);
        if (rho <= 0.0)
        {
            return 0.0;
        }
        return DynamicViscosity(h) / rho;
    }

    //  Reynolds number Re = V L / nu; 0 on degenerate viscosity.
    double Reynolds(double h, double v, double len) const
    {
        double nu = KinematicViscosity(h);
        if (nu <= 0.0)
        {
            return 0.0;
        }
        return v * len / nu;
    }

    //  Geometric -> geopotential height.
    static constexpr double GeomToGeopot(double h)
    {
        return R_GEOPOT * h / (R_GEOPOT + h);
    }

    //  Geopotential -> geometric height.
    static constexpr double GeopotToGeom(double hp)
    {
        return R_GEOPOT * hp / (R_GEOPOT - hp);
    }

private:

    //  One height, standard day: kinetic temperature, pressure, molar mass.
    struct Point
    {
        double tIsa;
        double p;
        double m;
    };

    //  Molar-temperature layer: base geopotential height, base molar
    //  temperature, gradient dT_M/dH, base pressure (computed from the
    //  layers below with the GOST constants).
    struct Layer
    {
        double hb;
        double tb;
        double lb;
        double pb;
    };

    static constexpr int N_LAY = 10;
    static constexpr double Z_MIN = -1999.3709472292348;    // GOST lower bound -2000 m geopotential, as geometric, m
    static constexpr double Z_M = 94000.0;          // molar mass constant below, geometric, m
    static constexpr int N_KNOT = 3;                // 120, 150, 200 km (Knot)

    //  GOST 4401-81 table 5 (T_M by geopotential height to 120 km geometric).
    static const Layer& Lay(int i)
    {
        static const Layer L[N_LAY] =
        {
            { 0.0,      288.15, -0.0065, 101325.0 },
            { 11000.0,  216.65,  0.0,    22632.040095007793 },
            { 20000.0,  216.65,  0.0010, 5474.877424281046 },
            { 32000.0,  228.65,  0.0028, 868.0157766202153 },
            { 47000.0,  270.65,  0.0,    110.9057733673104 },
            { 51000.0,  270.65, -0.0028, 66.93852812117996 },
            { 71000.0,  214.65, -0.0020, 3.956392160396618 },
            { 85000.0,  186.65,  0.0,    0.3634093782341348 },
            { 94000.0,  186.65,  0.0030, 0.06998127390940209 },
            { 102450.0, 212.00,  0.0110, 0.01641157496321042 }
        };
        return L[i];
    }

    //  Layer whose base is at or below H (layer 0 below its base).
    static int LayerIndex(double hp)
    {
        int i = N_LAY - 1;
        int k = 0;
        for (k = 0; k < N_LAY - 1; k = k + 1)
        {
            if (hp < Lay(k + 1).hb)
            {
                i = k;
                break;
            }
        }
        return i;
    }

    //  Molar mass, kg/kmol, GOST section 4 (geometric height h, m):
    //  constant to 94 km; 94-97 km
    //    M = 28.82 + 0.158 sqrt(1 - 7.5e-8 (h - 94000)^2) - 2.479e-4 sqrt(97000 - h);
    //  linear 97-97.5 km (-0.00012 /m) and 97.5-120 km (-0.0001511 /m);
    //  120-250 km the cubic of table 3. upper -- the cubic at any h (the
    //  pieces above 120 km), else the formulas below 120 km at any h.
    static double MolarMassOf(double h, bool upper)
    {
        if (upper)
        {
            return MolarMassCubic(h);
        }
        if (h <= Z_M)
        {
            return M_AIR0;
        }
        if (h <= 97000.0)
        {
            double d = h - 94000.0;
            return 28.82 + 0.158 * sqrt(1.0 - 7.5e-8 * d * d) - 2.479e-4 * sqrt(97000.0 - h);
        }
        if (h <= 97500.0)
        {
            return 28.91 - 0.00012 * (h - 97000.0);
        }
        return 28.85 - 0.0001511 * (h - 97500.0);
    }

    //  Molar mass 120-250 km, kg/kmol, GOST table 3, interval 1.
    static double MolarMassCubic(double h)
    {
        double b = -1.85675e-15;
        b = b * h + 12.08693e-10;
        b = b * h - 29.71210e-5;
        b = b * h + 46.9083;
        return b;
    }

    //  Kinetic temperature above 120 km, K, GOST table 6 (geometric h, m).
    static double TempUpper(double h)
    {
        if (h < 140000.0)
        {
            return 334.42 + 0.011259 * (h - 120000.0);
        }
        if (h < 160000.0)
        {
            return 559.60 + 0.006800 * (h - 140000.0);
        }
        if (h < 200000.0)
        {
            return 695.60 + 0.003970 * (h - 160000.0);
        }
        return 854.40 + 0.001750 * (h - 200000.0);
    }

    //  Particle concentration above 120 km, 1/m^3, GOST table 7, eq. (9):
    //  n = (c0 + c1 h + c2 h^2 + c3 h^3 + c4 h^4) 10^e (geometric h, m).
    //  c0 of 150-200 and 200-250 km fitted to table 4 (see the header).
    //    i -- interval: 0 - 120-150 km, 1 - 150-200 km, 2 - 200-250 km
    static double ConcUpper(double h, int i)
    {
        static const double C[3][5] =
        {
            { 0.210005867e4,  -0.5618444757e-1, 0.5663986231e-6, -0.2547466858e-11, 0.4309844119e-17 },
            { 0.1016393796e4, -0.2119530830e-1, 0.1671627815e-6, -0.5894237068e-12, 0.7826684089e-18 },
            { 0.7631575226e3, -0.1150600844e-1, 0.6612598428e-7, -0.1708736137e-12, 0.1669823114e-18 }
        };
        static const double E[3] = { 1.0e17, 1.0e16, 1.0e15 };
        double s = C[i][4];
        s = s * h + C[i][3];
        s = s * h + C[i][2];
        s = s * h + C[i][1];
        s = s * h + C[i][0];
        return s * E[i];
    }

    //  1.9.0: the knots where GOST does not join its pieces, geometric, m.
    static double Knot(int k)
    {
        static const double K[N_KNOT] = { 120000.0, 150000.0, 200000.0 };
        return K[k];
    }

    //  Piece of the profile at a height: 0 - barometric (below 120 km);
    //  1, 2, 3 - the concentration intervals 120-150, 150-200, 200-250 km.
    static int PieceOf(double h)
    {
        int i = 0;
        int k = 0;
        for (k = 0; k < N_KNOT; k = k + 1)
        {
            if (h >= Knot(k))
            {
                i = k + 1;
            }
        }
        return i;
    }

    //  One piece of the profile at a height, standard day; a piece is
    //  evaluated as is a little outside its interval too (Smooth blend).
    //  The pressure is computed only when asked (wantP).
    //    h     -- geometric height, m (not below Z_MIN)
    //    piece -- see PieceOf
    static void EvalPiece(double h, int piece, bool wantP, Point* q)
    {
        q->p = 0.0;
        if (piece > 0)
        {
            double z = h;
            if (z > Z_MAX)
            {
                z = Z_MAX;
            }
            q->tIsa = TempUpper(z);
            q->m = MolarMassOf(z, true);
            if (wantP && h <= Z_MAX)
            {
                q->p = ConcUpper(z, piece - 1) * K_BOLTZ * q->tIsa;
            }
            return;
        }
        double hp = GeomToGeopot(h);
        int i = LayerIndex(hp);
        const Layer& L = Lay(i);
        double dh = hp - L.hb;
        double tm = L.tb + L.lb * dh;
        q->m = MolarMassOf(h, false);
        q->tIsa = tm;
        if (q->m != M_AIR0)
        {
            q->tIsa = tm * (q->m / M_AIR0);
        }
        if (!wantP)
        {
            return;
        }
        if (L.lb != 0.0)
        {
            q->p = L.pb * pow(tm / L.tb, -G0_STD / (R_AIR * L.lb));
        }
        else
        {
            q->p = L.pb * exp(-G0_STD * dh / (R_AIR * L.tb));
        }
    }

    //  Smooth: the pieces below and above knot k blended at height h,
    //  |h - Knot(k)| < W_SMOOTH. Weight s(x) = x^3 (10 - 15 x + 6 x^2),
    //  x = (h - Knot(k) + W) / (2 W): s, s', s'' are 0 / 1, 0, 0 at the
    //  window ends, so the profile is C2 there. Temperature and molar mass
    //  are blended linearly, pressure in its logarithm (the density then
    //  follows from the state equation).
    static void Blend(double h, int k, bool wantP, Point* q)
    {
        Point a { };
        Point b { };
        EvalPiece(h, k, wantP, &a);
        EvalPiece(h, k + 1, wantP, &b);
        double x = (h - Knot(k) + W_SMOOTH) / (2.0 * W_SMOOTH);
        double s = x * x * x * (10.0 - 15.0 * x + 6.0 * x * x);
        q->tIsa = a.tIsa + s * (b.tIsa - a.tIsa);
        q->m = a.m + s * (b.m - a.m);
        q->p = 0.0;
        if (wantP)
        {
            q->p = a.p * exp(s * log(b.p / a.p));
        }
    }

    //  The profile at one height (standard day), in the mode of the object.
    //    h -- geometric height, m
    void Eval(double h, bool wantP, Point* q) const
    {
        if (h < Z_MIN)
        {
            h = Z_MIN;
        }
        if (m_mode == AtmMode::Smooth)
        {
            int k = 0;
            for (k = 0; k < N_KNOT; k = k + 1)
            {
                double d = h - Knot(k);
                if (d > -W_SMOOTH && d < W_SMOOTH)
                {
                    Blend(h, k, wantP, q);
                    return;
                }
            }
        }
        EvalPiece(h, PieceOf(h), wantP, q);
    }

    //  Standard temperature + offset, not below 1e-3 K.
    double OffsetT(double tIsa) const
    {
        double t = tIsa + m_dt;
        if (t < 1.0e-3)
        {
            t = 1.0e-3;
        }
        return t;
    }

    //  State equation: rho = p M / (R* T), written as p (M / M_AIR0) / (R T)
    //  so that below 94 km it is p / (R T) to the bit.
    static double Rho(double p, double m, double t)
    {
        if (m == M_AIR0)
        {
            return p / (R_AIR * t);
        }
        return p * (m / M_AIR0) / (R_AIR * t);
    }

    double m_dt;
    AtmMode m_mode;
};

//  Constant-initialized standard-day atmosphere.
inline constexpr EnvAtmosphere ATM_STD = EnvAtmosphere(0.0);

//  1.9.0: the same, blended at the knots of the standard (for integrators).
inline constexpr EnvAtmosphere ATM_STD_SMOOTH = EnvAtmosphere(0.0, AtmMode::Smooth);

}   // namespace gnc

#endif  // LIB_ATMOSPHERE_H

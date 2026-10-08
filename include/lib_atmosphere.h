//======================================================================
//  lib_atmosphere.h
//  EnvAtmosphere - standard atmosphere from sea level to 150 km, with a
//  temperature offset from the standard (hot/cold day: temperature shifted,
//  pressure profile kept, density recomputed from the state equation).
//
//  PROFILE (1.0.4)
//    0 ... 94 km GEOPOTENTIAL -- GOST 4401-81, its layers exactly: the same as
//      US-76 up to 71 km; the -2 K/km layer ends at 85 km (US-76: 84.852 km),
//      then 186.65 K isothermal up to 94 km.
//    94 km geopotential (95.411 km geometric) ... 150 km geometric --
//      TEMPORARILY US-76 (owner: to be replaced by the GOST 4401-81 table):
//      its tabulated temperature, pressure and density every 1 km, pressure
//      and density scaled so they join GOST continuously at 94 km (x0.990 and
//      x1.010), the temperature blended from 186.65 K into US-76 with a 5 km
//      decay. Pressure and density are log-linear between the rows.
//    above 150 km -- no air: pressure and density 0.
//
//  HEIGHT. Every function takes the GEOMETRIC height, m -- the height above
//  the WGS84 ellipsoid. The layers up to 94 km are geopotential and are read
//  through H = R0 h / (R0 + h), R0 = R_GEOPOT = 6 356 766 m; the table above
//  is geometric already.
//
//  1.0.2 put the geometric height into the geopotential layers without
//  conversion (pressure 5 % low at 51.6 km); 1.0.3 converted it and ended at
//  86 km.
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

inline constexpr double G0_STD = 9.80665;       // standard gravity, m/s^2
inline constexpr double R_AIR = 287.0528;       // gas constant of air, J/(kg K)
inline constexpr double R_GEOPOT = 6356766.0;   // geopotential radius, m

//  The air at one height (1.4.0): what Temperature, Pressure and Density
//  return there (EnvAtmosphere::State).
struct AtmState
{
    double t;     // temperature, K (with the offset)
    double p;     // static pressure, Pa
    double rho;   // density, kg/m^3
};

class EnvAtmosphere
{
public:

    static constexpr double P_SL = 101325.0;    // sea-level pressure, Pa

    //  dTemp - temperature offset from ISA, K (0 - standard day).
    constexpr explicit EnvAtmosphere(double dTemp = 0.0)
        : m_lay()
        , m_dt(dTemp)
    {
        //  GOST 4401-81 layers (the same as US-76 below 71 km): base height,
        //  base temperature, lapse rate, base pressure. Fields BY NAME.
        m_lay[0].hb = 0.0;
        m_lay[0].tb = 288.15;
        m_lay[0].lb = -0.0065;
        m_lay[0].pb = 101325.0;
        m_lay[1].hb = 11000.0;
        m_lay[1].tb = 216.65;
        m_lay[1].lb = 0.0;
        m_lay[1].pb = 22632.06;
        m_lay[2].hb = 20000.0;
        m_lay[2].tb = 216.65;
        m_lay[2].lb = 0.001;
        m_lay[2].pb = 5474.889;
        m_lay[3].hb = 32000.0;
        m_lay[3].tb = 228.65;
        m_lay[3].lb = 0.0028;
        m_lay[3].pb = 868.0187;
        m_lay[4].hb = 47000.0;
        m_lay[4].tb = 270.65;
        m_lay[4].lb = 0.0;
        m_lay[4].pb = 110.9063;
        m_lay[5].hb = 51000.0;
        m_lay[5].tb = 270.65;
        m_lay[5].lb = -0.0028;
        m_lay[5].pb = 66.93887;
        m_lay[6].hb = 71000.0;
        m_lay[6].tb = 214.65;
        m_lay[6].lb = -0.002;
        m_lay[6].pb = 3.956420;
        //  GOST 4401-81: the -2 K/km layer ends at 85 km, then isothermal.
        m_lay[7].hb = 85000.0;
        m_lay[7].tb = 186.65;
        m_lay[7].lb = 0.0;
        m_lay[7].pb = 0.363411935;
    }

    constexpr void SetDeltaTemp(double dTemp)
    {
        m_dt = dTemp;
    }

    constexpr double DeltaTemp() const
    {
        return m_dt;
    }

    //  Air temperature, K; clamped below 0 m; above 150 km held at its value there.
    //  Includes the ISA offset.
    double Temperature(double h) const
    {
        double t = TemperatureIsa(h) + m_dt;
        if (t < 1.0e-3)
        {
            t = 1.0e-3;
        }
        return t;
    }

    //  Static pressure, Pa: hydrostatic integration inside a layer,
    //  power law with lapse, exponent in isothermal layers; above the
    //  GOST layers the table; 0 above 150 km. Does NOT depend on dTemp (non-standard-day convention:
    //  pressure profile kept, temperature shifted).
    //    h -- geometric height, m; converted to geopotential inside
    double Pressure(double h) const
    {
        if (h <= 0.0)
        {
            h = 0.0;
        }
        if (h >= Z_MAX)
        {
            return 0.0;
        }
        if (h >= Z_TOP)
        {
            double f = 0.0;
            int k = UpIndex(h, &f);
            return K_P * UpLog(UpP(k), UpP(k + 1), f);
        }
        h = GeomToGeopot(h);
        int i = LayerIndex(h);
        double hb = m_lay[i].hb;
        double tb = m_lay[i].tb;
        double lb = m_lay[i].lb;
        double pb = m_lay[i].pb;
        double dh = h - hb;
        if (fabs(lb) > 1.0e-12)
        {
            double t = tb + lb * dh;
            if (t < 1.0e-3)
            {
                t = 1.0e-3;
            }
            return pb * pow(t / tb, -G0_STD / (R_AIR * lb));
        }
        return pb * exp(-G0_STD * dh / (R_AIR * tb));
    }

    //  Density, kg/m^3, from the state equation rho = P / (R T);
    //  0 above 150 km. Depends on dTemp through the temperature.
    double Density(double h) const
    {
        if (h <= 0.0)
        {
            h = 0.0;
        }
        if (h >= Z_MAX)
        {
            return 0.0;
        }
        if (h >= Z_TOP)
        {
            //  The table density (the molar mass falls up there), scaled to
            //  GOST at the joint, and moved with the temperature offset as
            //  the state equation moves it.
            double f = 0.0;
            int k = UpIndex(h, &f);
            return K_R * UpLog(UpRho(k), UpRho(k + 1), f) * TemperatureIsa(h) / Temperature(h);
        }
        double t = Temperature(h);
        double p = Pressure(h);
        return p / (R_AIR * t);
    }

    //  1.4.0: temperature, pressure and density at once, with one height
    //  conversion and one layer search -- the same numbers as Temperature,
    //  Pressure and Density, to the bit (the three calls convert the height
    //  and search the layers four times).
    //    h -- geometric height, m
    AtmState State(double h) const
    {
        AtmState r { };
        if (h <= 0.0)
        {
            h = 0.0;
        }
        if (h >= Z_TOP)
        {
            double z = h < Z_MAX ? h : Z_MAX;
            double f = 0.0;
            int k = UpIndex(z, &f);
            double tu = UpT(k) + f * (UpT(k + 1) - UpT(k));
            double tIsa = tu + (T_TOP - UpT(0)) * exp(-(z - Z_TOP) / 5000.0);
            r.t = tIsa + m_dt;
            if (r.t < 1.0e-3)
            {
                r.t = 1.0e-3;
            }
            if (h >= Z_MAX)
            {
                return r;
            }
            r.p = K_P * UpLog(UpP(k), UpP(k + 1), f);
            r.rho = K_R * UpLog(UpRho(k), UpRho(k + 1), f) * tIsa / r.t;
            return r;
        }
        double hg = GeomToGeopot(h);
        int i = LayerIndex(hg);
        double tb = m_lay[i].tb;
        double lb = m_lay[i].lb;
        double dh = hg - m_lay[i].hb;
        double tl = tb + lb * dh;
        if (tl < 1.0e-3)
        {
            tl = 1.0e-3;
        }
        r.t = tl + m_dt;
        if (r.t < 1.0e-3)
        {
            r.t = 1.0e-3;
        }
        if (fabs(lb) > 1.0e-12)
        {
            r.p = m_lay[i].pb * pow(tl / tb, -G0_STD / (R_AIR * lb));
        }
        else
        {
            r.p = m_lay[i].pb * exp(-G0_STD * dh / (R_AIR * tb));
        }
        r.rho = r.p / (R_AIR * r.t);
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

    //  Speed of sound, m/s: a = sqrt(gamma R T), gamma = 1.4.
    double SoundSpeed(double h) const
    {
        return sqrt(1.4 * R_AIR * Temperature(h));
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

    //  Dynamic viscosity (Sutherland), Pa*s.
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
    static double GeomToGeopot(double h)
    {
        return R_GEOPOT * h / (R_GEOPOT + h);
    }

    //  Geopotential -> geometric height.
    static double GeopotToGeom(double hp)
    {
        return R_GEOPOT * hp / (R_GEOPOT - hp);
    }

private:

    struct Layer
    {
        double hb;
        double tb;
        double lb;
        double pb;
    };

    static constexpr int N_LAY = 8;

    //  GOST 4401-81 ends its layers at 94 km geopotential; the table above
    //  begins there, at the same geometric height.
    static constexpr double Z_TOP = 95410.878;        // geometric height of 94 km geopotential, m
    static constexpr double Z_MAX = 150000.0;   // no air above, m
    static constexpr double T_TOP = 186.65;     // GOST temperature at the top of its layers, K
    static constexpr double K_P = 0.990088207;           // US-76 pressure scaled to GOST at the joint
    static constexpr double K_R = 1.009883930;           // US-76 density scaled to GOST at the joint
    static constexpr int N_UP = 56;

    //  US-76 above the joint, geometric height every 1 km (TEMPORARY, owner).
    static double UpZ(int i)
    {
        static const double Z[N_UP] =
        {
            95410.878, 96000.000, 97000.000, 98000.000, 99000.000, 100000.000,
            101000.000, 102000.000, 103000.000, 104000.000, 105000.000, 106000.000,
            107000.000, 108000.000, 109000.000, 110000.000, 111000.000, 112000.000,
            113000.000, 114000.000, 115000.000, 116000.000, 117000.000, 118000.000,
            119000.000, 120000.000, 121000.000, 122000.000, 123000.000, 124000.000,
            125000.000, 126000.000, 127000.000, 128000.000, 129000.000, 130000.000,
            131000.000, 132000.000, 133000.000, 134000.000, 135000.000, 136000.000,
            137000.000, 138000.000, 139000.000, 140000.000, 141000.000, 142000.000,
            143000.000, 144000.000, 145000.000, 146000.000, 147000.000, 148000.000,
            149000.000, 150000.000
        };
        return Z[i];
    }
    static double UpT(int i)
    {
        static const double T[N_UP] =
        {
            188.7575, 189.3050, 190.4035, 191.7234, 193.2773, 195.0813,
            197.1559, 199.5274, 202.2305, 205.3117, 208.8352, 212.8939,
            217.6306, 223.2860, 230.3312, 239.9997, 252.0000, 264.0000,
            276.0000, 288.0000, 300.0000, 312.0000, 324.0000, 336.0000,
            348.0000, 360.0000, 371.8864, 383.5485, 394.9905, 406.2167,
            417.2313, 428.0383, 438.6416, 449.0452, 459.2528, 469.2683,
            479.0953, 488.7374, 498.1981, 507.4810, 516.5893, 525.5265,
            534.2957, 542.9003, 551.3433, 559.6279, 567.7570, 575.7336,
            583.5606, 591.2409, 598.7773, 606.1725, 613.4291, 620.5499,
            627.5374, 634.3941
        };
        return T[i];
    }
    static double UpP(int i)
    {
        static const double P[N_UP] =
        {
            7.068236e-02, 6.376231e-02, 5.357131e-02, 4.505906e-02, 3.794895e-02, 3.200574e-02,
            2.719660e-02, 2.314841e-02, 1.974340e-02, 1.688073e-02, 1.447451e-02, 1.245185e-02,
            1.075117e-02, 9.320610e-03, 8.116589e-03, 7.102788e-03, 6.262352e-03, 5.555474e-03,
            4.956891e-03, 4.446640e-03, 4.008837e-03, 3.630765e-03, 3.302183e-03, 3.014787e-03,
            2.761814e-03, 2.537377e-03, 2.340366e-03, 2.164232e-03, 2.006248e-03, 1.864094e-03,
            1.735789e-03, 1.619640e-03, 1.514193e-03, 1.418196e-03, 1.330567e-03, 1.250371e-03,
            1.176795e-03, 1.109131e-03, 1.046763e-03, 9.891506e-04, 9.358211e-04, 8.863591e-04,
            8.403991e-04, 7.976186e-04, 7.577327e-04, 7.204894e-04, 6.856648e-04, 6.530602e-04,
            6.224989e-04, 5.938229e-04, 5.668915e-04, 5.415785e-04, 5.177710e-04, 4.953675e-04,
            4.742768e-04, 4.541520e-04
        };
        return P[i];
    }
    static double UpRho(int i)
    {
        static const double D[N_UP] =
        {
            1.293372e-06, 1.162025e-06, 9.685670e-07, 8.071057e-07, 6.725008e-07, 5.601843e-07,
            4.695725e-07, 3.934839e-07, 3.298591e-07, 2.767593e-07, 2.324417e-07, 1.953891e-07,
            1.643117e-07, 1.381326e-07, 1.159663e-07, 9.706754e-08, 8.113322e-08, 6.839327e-08,
            5.811486e-08, 4.974959e-08, 4.288343e-08, 3.720118e-08, 3.246087e-08, 2.847535e-08,
            2.509881e-08, 2.220555e-08, 1.976862e-08, 1.767169e-08, 1.585874e-08, 1.428415e-08,
            1.291058e-08, 1.170727e-08, 1.064885e-08, 9.714208e-09, 8.885781e-09, 8.148848e-09,
            7.491040e-09, 6.901922e-09, 6.372656e-09, 5.895733e-09, 5.464748e-09, 5.074223e-09,
            4.719453e-09, 4.396390e-09, 4.101536e-09, 3.831862e-09, 3.584738e-09, 3.357870e-09,
            3.149257e-09, 2.957144e-09, 2.779991e-09, 2.616440e-09, 2.465295e-09, 2.325495e-09,
            2.196099e-09, 2.075208e-09
        };
        return D[i];
    }

    //  Row below a geometric height of the upper table, and the fraction.
    static int UpIndex(double z, double* f)
    {
        //  1.4.0: the rows go every kilometre from UpZ(1); the row is
        //  guessed from that and then moved by the same test as the search
        //  from the bottom, so it is the same row (was a search over up to
        //  56 rows at every call above 95 km).
        int i = z >= UpZ(1) ? (int)((z - UpZ(1)) / 1000.0) + 1 : 0;
        i = i < 0 ? 0 : (i > N_UP - 2 ? N_UP - 2 : i);
        while (i > 0 && z < UpZ(i))
        {
            i = i - 1;
        }
        while (i < N_UP - 2 && z >= UpZ(i + 1))
        {
            i = i + 1;
        }
        double span = UpZ(i + 1) - UpZ(i);
        *f = span > 1.0e-9 ? (z - UpZ(i)) / span : 0.0;
        *f = *f < 0.0 ? 0.0 : (*f > 1.0 ? 1.0 : *f);
        return i;
    }

    //  Upper table, log-linear.
    static double UpLog(double a, double b, double f)
    {
        return a > 0.0 && b > 0.0 ? a * exp(f * log(b / a)) : a + f * (b - a);
    }

    //  Layer whose base is at or below h (the last layer above all bases).
    int LayerIndex(double h) const
    {
        int i = N_LAY - 1;
        int k = 0;
        for (k = 0; k < N_LAY - 1; k = k + 1)
        {
            if (h < m_lay[k + 1].hb)
            {
                i = k;
                break;
            }
        }
        return i;
    }

    //  ISA temperature without the offset (internal, for the pressure).
    //    h -- geometric height, m; converted to geopotential inside
    double TemperatureIsa(double h) const
    {
        if (h <= 0.0)
        {
            h = 0.0;
        }
        if (h >= Z_TOP)
        {
            double z = h < Z_MAX ? h : Z_MAX;
            double f = 0.0;
            int k = UpIndex(z, &f);
            double tu = UpT(k) + f * (UpT(k + 1) - UpT(k));
            //  Blend from the GOST 186.65 K into US-76, 5 km decay.
            return tu + (T_TOP - UpT(0)) * exp(-(z - Z_TOP) / 5000.0);
        }
        h = GeomToGeopot(h);
        int i = LayerIndex(h);
        double t = m_lay[i].tb + m_lay[i].lb * (h - m_lay[i].hb);
        if (t < 1.0e-3)
        {
            t = 1.0e-3;
        }
        return t;
    }

    Layer m_lay[N_LAY];
    double m_dt;
};

//  Constant-initialized standard-day atmosphere.
inline constexpr EnvAtmosphere ATM_STD = EnvAtmosphere(0.0);

}   // namespace gnc

#endif  // LIB_ATMOSPHERE_H

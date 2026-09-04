//======================================================================
//  lib_atmosphere.h
//  EnvAtmosphere - layered US Standard Atmosphere 1976, valid to ~86 km, with a
//  temperature offset from ISA (hot/cold day: temperature shifted,
//  layered pressure kept, density recomputed from the state equation).
//  Height is treated as geometric (error < 0.5% below 50 km); exact
//  geometric <-> geopotential conversions are provided.
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

class EnvAtmosphere
{
public:

    static constexpr double P_SL = 101325.0;    // sea-level pressure, Pa

    //  dTemp - temperature offset from ISA, K (0 - standard day).
    constexpr explicit EnvAtmosphere(double dTemp = 0.0)
        : m_lay()
        , m_dt(dTemp)
    {
        //  US-76 layers: base height, base temperature, lapse rate,
        //  base pressure. Fields are assigned BY NAME.
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
    }

    constexpr void SetDeltaTemp(double dTemp)
    {
        m_dt = dTemp;
    }

    constexpr double DeltaTemp() const
    {
        return m_dt;
    }

    //  Air temperature, K; clamped below 0 m; constant 186 K above 86 km.
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
    //  power law with lapse, exponent in isothermal layers; 0 above
    //  86 km. Does NOT depend on dTemp (non-standard-day convention:
    //  pressure profile kept, temperature shifted).
    double Pressure(double h) const
    {
        if (h <= 0.0)
        {
            h = 0.0;
        }
        if (h >= 86000.0)
        {
            return 0.0;
        }
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
    //  0 above 86 km. Depends on dTemp through the temperature.
    double Density(double h) const
    {
        if (h <= 0.0)
        {
            h = 0.0;
        }
        if (h >= 86000.0)
        {
            return 0.0;
        }
        double t = Temperature(h);
        double p = Pressure(h);
        return p / (R_AIR * t);
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

    static constexpr int N_LAY = 7;

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
    double TemperatureIsa(double h) const
    {
        if (h <= 0.0)
        {
            h = 0.0;
        }
        if (h >= 86000.0)
        {
            return 186.0;
        }
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

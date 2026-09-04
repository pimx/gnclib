//======================================================================
//  lib_wind.h
//  EnvWind - wind in the local axes: constant mean plus gusts. Gusts are a
//  first-order Gauss-Markov process per axis:
//      dw/dt = -w/tau + sigma sqrt(2/tau) n(t)
//  (a simplified Dryden analogue, adequate at low speeds; vertical
//  gusts attenuated by 0.3). Exact discretization of the OU process:
//      a = exp(-dt/tau), s = sigma sqrt(1 - a^2).
//  Static height profiles of the mean wind: power law and log law.
//  C++17, no STL, no heap. EnvWind is constexpr-constructible.
//======================================================================

#ifndef LIB_WIND_H
#define LIB_WIND_H

#include <math.h>

#include "lib_frames.h"
#include "lib_random.h"

namespace gnc
{

class EnvWind
{
public:

    constexpr EnvWind()
        : m_mean()
        , m_gust()
        , m_sigma(0.0)
        , m_tau(1.0)
        , m_rng(123456789ULL)
        , m_time(0.0)
    {
    }

    //  mean  - mean wind in the local axes, m/s
    //  sigma - gust standard deviation, m/s
    //  tau   - gust time constant, s (clamped from below)
    //  seed  - RNG seed, t0 - start time, s
    void Setup(const CfLocalVector& mean, double sigma, double tau,
               unsigned long long seed, double t0)
    {
        m_mean = mean;
        m_sigma = sigma;
        m_tau = tau;
        if (m_tau < 1.0e-3)
        {
            m_tau = 1.0e-3;
        }
        CfLocalVector z { };
        m_gust = z;
        m_rng.Reseed(seed);
        m_time = t0;
    }

    //  Advance the gust process to time t, s; returns the current wind
    //  velocity in the local axes.
    CfLocalVector Update(double t)
    {
        double dt = t - m_time;
        m_time = t;
        if (dt <= 0.0)
        {
            return Velocity();
        }
        double a = exp(-dt / m_tau);
        double s = m_sigma * sqrt(1.0 - a * a);
        m_gust.e = a * m_gust.e + s * m_rng.Gauss();
        m_gust.n = a * m_gust.n + s * m_rng.Gauss();
        m_gust.u = a * m_gust.u + s * (0.3 * m_rng.Gauss());
        return Velocity();
    }

    //  Current wind velocity without advancing, m/s.
    CfLocalVector Velocity() const
    {
        return Add(m_mean, m_gust);
    }

private:

    CfLocalVector m_mean;
    CfLocalVector m_gust;
    double m_sigma;
    double m_tau;
    RNG m_rng;
    double m_time;
};

//  Power-law mean wind profile: V(h) = vRef (h/hRef)^alpha;
//  alpha ~0.14 for open terrain. h and hRef clamped from below at 0.1 m.
inline double WindPowerLaw(double vRef, double hRef, double h, double alpha)
{
    if (h < 0.1)
    {
        h = 0.1;
    }
    if (hRef < 0.1)
    {
        hRef = 0.1;
    }
    return vRef * pow(h / hRef, alpha);
}

//  Logarithmic surface-layer profile: V(h) = vRef ln(h/z0)/ln(hRef/z0);
//  z0 - roughness length, m. Zero at or below z0.
inline double WindLogLaw(double vRef, double hRef, double h, double z0)
{
    if (z0 < 1.0e-4)
    {
        z0 = 1.0e-4;
    }
    if (h <= z0)
    {
        return 0.0;
    }
    if (hRef <= z0)
    {
        hRef = z0 * 2.0;
    }
    return vRef * log(h / z0) / log(hRef / z0);
}

}   // namespace gnc

#endif  // LIB_WIND_H

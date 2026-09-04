//======================================================================
//  lib_gravity.h
//  Gravity models over the named ECEF types. Field parameters (mu, j2,
//  om) are taken from the CfGeoid structure, so PZ-90.11 and WGS-84 are
//  interchangeable. C++17, no STL, no heap.
//======================================================================

#ifndef LIB_GRAVITY_H
#define LIB_GRAVITY_H

#include <math.h>

#include "lib_frames.h"
#include "lib_geodesy.h"

namespace gnc
{

//  Gravitational acceleration in ECEF: central field + J2,
//  WITHOUT the centrifugal term. Zero inside a 1 m guard sphere.
inline CfEcefVector GravitationJ2(const CfEcefPoint& p, const CfGeoid& el)
{
    double r2 = p.x * p.x + p.y * p.y + p.z * p.z;
    double r = sqrt(r2);
    CfEcefVector g { };
    if (r < 1.0)
    {
        return g;
    }
    double f = -el.mu / (r2 * r);
    double zr2 = (p.z * p.z) / r2;
    double j2c = 1.5 * el.j2 * (el.a * el.a) / r2;
    g.x = f * p.x * (1.0 - j2c * (5.0 * zr2 - 1.0));
    g.y = f * p.y * (1.0 - j2c * (5.0 * zr2 - 1.0));
    g.z = f * p.z * (1.0 - j2c * (5.0 * zr2 - 3.0));
    return g;
}

//  Centrifugal acceleration of a point fixed in ECEF, taken as the
//  apparent force in the rotating frame: -w x (w x r) = om^2 (x, y, 0).
inline CfEcefVector CentrifugalEcf(const CfEcefPoint& p, const CfGeoid& el)
{
    double w2 = el.om * el.om;
    CfEcefVector a { };
    a.x = w2 * p.x;
    a.y = w2 * p.y;
    a.z = 0.0;
    return a;
}

//  Full gravity of the rotating Earth (gravitation J2 + centrifugal) -
//  what a static accelerometer measures.
inline CfEcefVector GravityEcf(const CfEcefPoint& p, const CfGeoid& el)
{
    return Add(GravitationJ2(p, el), CentrifugalEcf(p, el));
}

//  Somigliana normal gravity on the WGS-84 ellipsoid plus a second-order
//  free-air height correction, m/s^2. Scalar along the normal.
//  The constants are WGS-84 SPECIFIC; for PZ-90.11 use the magnitude of
//  GravityEcf, which follows the ellipsoid passed in.
inline double GravityNormalWgs84(double lat, double h)
{
    double s2 = sin(lat);
    s2 = s2 * s2;
    double g0 = 9.7803253359 * (1.0 + 0.00193185265241 * s2) /
                sqrt(1.0 - 0.00669437999013 * s2);
    return g0 - (3.0877e-6 - 4.4e-9 * s2) * h + 7.2e-13 * h * h;
}

//  Gravity vector of the site in the local axes: GravityEcf at the
//  origin of the local frame, rotated into CfLocalVector (~(0, 0, -9.81)).
//  The site is static, the vector is constant - compute once at task
//  initialization and store where needed.
inline CfLocalVector SiteGravityLoc(const GEO& g)
{
    return g.E2L(GravityEcf(g.OriginEcp(), g.Ellipsoid()));
}

}   // namespace gnc

#endif  // LIB_GRAVITY_H

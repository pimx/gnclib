//============================================================================
// lib_earth.h -- rotating Earth for a site-fixed (Earth-fixed) frame
// (gnclib 1.1.0). The local frame of a GEO site is ECEF with the origin
// moved to the site and the axes rotated, so integrating in it with the
// terms below is integrating in ECEF.
//   EarthGravityLoc -- gravitation (J2) + centrifugal at a local point;
//   EarthRateLoc    -- Earth angular rate in the local axes;
//   CoriolisAccLoc  -- Coriolis acceleration -2 Omega x v;
//   SpecificForceLoc -- what accelerometers sense for the relative
//                       acceleration a_rel: a_rel + 2 Omega x v - g(r).
// C++17, no STL, no heap.
//============================================================================
#ifndef GNCLIB_LIB_EARTH_H
#define GNCLIB_LIB_EARTH_H

#include "lib_geodesy.h"
#include "lib_gravity.h"

namespace gnc
{

inline CfLocalVector EarthGravityLoc(const GEO& f, const CfLocalPoint& p)
{
    return f.E2L(GravityEcf(f.L2E(p), f.Ellipsoid()));
}

inline CfLocalVector EarthRateLoc(const GEO& f)
{
    double lat = f.Origin().lat;
    CfLocalVector w { };
    w.e = 0.0;
    w.n = OMEGA_E * cos(lat);
    w.u = OMEGA_E * sin(lat);
    return w;
}

inline CfLocalVector CoriolisAccLoc(const CfLocalVector& omegaL, const CfLocalVector& v)
{
    CfLocalVector c { };
    c.e = -2.0 * (omegaL.n * v.u - omegaL.u * v.n);
    c.n = -2.0 * (omegaL.u * v.e - omegaL.e * v.u);
    c.u = -2.0 * (omegaL.e * v.n - omegaL.n * v.e);
    return c;
}

inline CfLocalVector SpecificForceLoc(const GEO& f, const CfLocalPoint& p,
                                      const CfLocalVector& v, const CfLocalVector& aRel)
{
    CfLocalVector g = EarthGravityLoc(f, p);
    CfLocalVector c = CoriolisAccLoc(EarthRateLoc(f), v);
    CfLocalVector r { };
    r.e = aRel.e - c.e - g.e;
    r.n = aRel.n - c.n - g.n;
    r.u = aRel.u - c.u - g.u;
    return r;
}

}   // namespace gnc

#endif // GNCLIB_LIB_EARTH_H

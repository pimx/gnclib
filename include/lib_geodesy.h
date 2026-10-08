//======================================================================
//  lib_geodesy.h
//  GEO - operator between ECEF (x,y,z) and CfLocalVector (n,e,u).
//  Holds the geodetic origin of the local frame, therefore it
//  distinguishes free vectors (rotation only) from points
//  (rotation plus translation of the origin).
//  C++17, no STL, no heap. GEO is a literal type.
//======================================================================

#ifndef LIB_GEODESY_H
#define LIB_GEODESY_H

#include <math.h>

#include "lib_frames.h"

namespace gnc
{

//----------------------------------------------------------------------
//  GCM - direction cosine matrix with NAMED elements.
//  Element  A_b  is the cosine between ECEF axis A and CfLocalVector axis b.
//
//      ECEF.A = sum over b of ( A_b * CfLocalVector.b )      (direct, L2E)
//      CfLocalVector.b  = sum over A of ( A_b * ECEF.A )     (transpose, E2L)
//----------------------------------------------------------------------
struct GCM
{
    double x_n;
    double x_e;
    double x_u;
    double y_n;
    double y_e;
    double y_u;
    double z_n;
    double z_e;
    double z_u;
};

//  Origin on the equator at the prime meridian: lat = 0, lon = 0.
//      east  = ( 0, 1, 0 )
//      north = ( 0, 0, 1 )
//      up    = ( 1, 0, 0 )
constexpr GCM MakeGcmNull()
{
    GCM c { };
    c.x_u = 1.0;
    c.y_e = 1.0;
    c.z_n = 1.0;
    return c;
}

constexpr CfEcefPoint MakeEcpNull()
{
    CfEcefPoint p { };
    p.x = 6378136.0;
    return p;
}

//----------------------------------------------------------------------
//  Geodetic <-> ECEF, closed form (Bowring) for the inverse problem.
//----------------------------------------------------------------------

//  Radius of curvature in the prime vertical N(lat), m.
inline double PrimeVerticalRadius(double lat, const CfGeoid& el)
{
    double e2 = el.f * (2.0 - el.f);
    double sb = sin(lat);
    return el.a / sqrt(1.0 - e2 * sb * sb);
}

//  Radius of curvature in the meridian M(lat), m.
inline double MeridianRadius(double lat, const CfGeoid& el)
{
    double e2 = el.f * (2.0 - el.f);
    double sb = sin(lat);
    double w2 = 1.0 - e2 * sb * sb;
    return el.a * (1.0 - e2) / (w2 * sqrt(w2));
}

//  Geocentric latitude of a point ON the ellipsoid surface.
inline double GeocentricLat(double lat, const CfGeoid& el)
{
    double e2 = el.f * (2.0 - el.f);
    return atan((1.0 - e2) * tan(lat));
}

inline CfEcefPoint GeoToEcp(const CfGeodetic& g, const CfGeoid& el)
{
    double e2 = el.f * (2.0 - el.f);
    double sb = sin(g.lat);
    double cb = cos(g.lat);
    double sl = sin(g.lon);
    double cl = cos(g.lon);
    double rn = el.a / sqrt(1.0 - e2 * sb * sb);
    CfEcefPoint p { };
    p.x = (rn + g.alt) * cb * cl;
    p.y = (rn + g.alt) * cb * sl;
    p.z = (rn * (1.0 - e2) + g.alt) * sb;
    return p;
}

inline CfGeodetic EcpToGeo(const CfEcefPoint& p, const CfGeoid& el)
{
    double e2 = el.f * (2.0 - el.f);
    double bb = el.a * (1.0 - el.f);
    double ep2 = e2 / (1.0 - e2);
    double rp = sqrt(p.x * p.x + p.y * p.y);
    CfGeodetic g { };
    g.lon = atan2(p.y, p.x);
    if (rp < 1.0e-6)
    {
        //  Polar case: latitude is defined by the sign of z only.
        double sg = 1.0;
        if (p.z < 0.0)
        {
            sg = -1.0;
        }
        g.lat = sg * 1.57079632679489661923;
        g.alt = fabs(p.z) - bb;
        return g;
    }
    double th = atan2(p.z * el.a, rp * bb);
    double st = sin(th);
    double ct = cos(th);
    g.lat = atan2(p.z + ep2 * bb * st * st * st, rp - e2 * el.a * ct * ct * ct);
    //  Two refinement steps: Bowring alone loses accuracy at high
    //  altitude, the fixed point below converges to machine precision.
    int it = 0;
    double sb = 0.0;
    double rn = 0.0;
    for (it = 0; it < 2; it = it + 1)
    {
        sb = sin(g.lat);
        rn = el.a / sqrt(1.0 - e2 * sb * sb);
        g.alt = rp / cos(g.lat) - rn;
        g.lat = atan2(p.z, rp * (1.0 - e2 * rn / (rn + g.alt)));
    }
    sb = sin(g.lat);
    rn = el.a / sqrt(1.0 - e2 * sb * sb);
    g.alt = rp / cos(g.lat) - rn;
    return g;
}

//  1.5.0: the geodetic latitude as its sine and cosine, and the height.
struct CfGeoSc
{
    double sinLat;  // sine of the geodetic latitude
    double cosLat;  // cosine of the geodetic latitude (not negative)
    double alt;     // height above the ellipsoid, m
};

//  1.5.0: the latitude and the height of EcpToGeo without trigonometry.
//  The iteration of EcpToGeo (Bowring, two refinements, the final height)
//  with the latitude carried as its sine and cosine -- the direction (num,
//  den) normalised -- instead of atan2, sin and cos, and without the
//  longitude; the latitude returned is the one of the final height, as in
//  EcpToGeo. On the axis: the sine is the sign of z, the height |z| - b.
inline CfGeoSc EcpToGeoSc(const CfEcefPoint& p, const CfGeoid& el)
{
    CfGeoSc g { };
    double e2 = el.f * (2.0 - el.f);
    double bb = el.a * (1.0 - el.f);
    double ep2 = e2 / (1.0 - e2);
    double rp = sqrt(p.x * p.x + p.y * p.y);
    if (rp < 1.0e-6)
    {
        g.sinLat = p.z < 0.0 ? -1.0 : 1.0;
        g.cosLat = 0.0;
        g.alt = fabs(p.z) - bb;
        return g;
    }
    double ta = p.z * el.a;
    double tb = rp * bb;
    double nt = sqrt(ta * ta + tb * tb);
    double st = ta / nt;
    double ct = tb / nt;
    double num = p.z + ep2 * bb * st * st * st;
    double den = rp - e2 * el.a * ct * ct * ct;
    int it = 0;
    for (it = 0; it < 3; it = it + 1)
    {
        double nl = sqrt(num * num + den * den);
        g.sinLat = num / nl;
        g.cosLat = den / nl;
        double rn = el.a / sqrt(1.0 - e2 * g.sinLat * g.sinLat);
        g.alt = rp / g.cosLat - rn;
        num = p.z;
        den = rp * (1.0 - e2 * rn / (rn + g.alt));
    }
    return g;
}

//  1.4.0: the geodetic height only: EcpToGeoSc (1.5.0; the same steps as
//  in 1.4.0, the same numbers). MEASURED against EcpToGeo().alt on
//  PZ-90.11, latitude 0 ... 90 deg by 0.1 deg, height -2 ... 398 km by
//  1 km: at most 7.5e-9 m (rounding); about seven times cheaper.
inline double EcpToAlt(const CfEcefPoint& p, const CfGeoid& el)
{
    return EcpToGeoSc(p, el).alt;
}

//  1.5.0: the distance from the centre of the point at a geodetic latitude
//  (its sine and cosine) and a height: the length of GeoToEcp, which the
//  longitude does not change, without trigonometry.
//    sinLat, cosLat -- latitude
//    alt            -- height above the ellipsoid, m
//    el             -- ellipsoid
inline double GeoRadius(double sinLat, double cosLat, double alt, const CfGeoid& el)
{
    double e2 = el.f * (2.0 - el.f);
    double rn = el.a / sqrt(1.0 - e2 * sinLat * sinLat);
    double rho = (rn + alt) * cosLat;
    double z = (rn * (1.0 - e2) + alt) * sinLat;
    return sqrt(rho * rho + z * z);
}

//----------------------------------------------------------------------
//  GEO - the ECEF <-> CfLocalVector operator.
//----------------------------------------------------------------------
class GEO
{
public:

    //  Constant-initialized default: PZ-90.11, origin lat = 0, lon = 0,
    //  alt = 0. Call SetOrigin() at run time before use.
    constexpr GEO()
        : el(MakeEllPz9011())
        , org()
        , orp(MakeEcpNull())
        , m(MakeGcmNull())
    {
    }

    void SetEllipsoid(const CfGeoid& e)
    {
        el = e;
        SetOrigin(org);
    }

    void SetOrigin(const CfGeodetic& g)
    {
        org = g;
        orp = GeoToEcp(g, el);
        double sb = sin(g.lat);
        double cb = cos(g.lat);
        double sl = sin(g.lon);
        double cl = cos(g.lon);
        //  Columns of the matrix are the CfLocalVector axes written in ECEF:
        //      north = ( -sb*cl, -sb*sl,  cb )
        //      east  = ( -sl,     cl,     0  )
        //      up    = (  cb*cl,  cb*sl,  sb )
        m.x_n = -sb * cl;
        m.x_e = -sl;
        m.x_u = cb * cl;
        m.y_n = -sb * sl;
        m.y_e = cl;
        m.y_u = cb * sl;
        m.z_n = cb;
        m.z_e = 0.0;
        m.z_u = sb;
    }

    //-- transforms of free vectors (rotation only) -------------------
    //  One implementation, typed entry points per quantity kind:
    //  geometric vector, angular rate, angular acceleration.

    CfEcefVector L2E(const CfLocalVector& l) const
    {
        return MapL2E<CfEcefVector>(l);
    }

    CfEcefOmega L2E(const CfLocalOmega& l) const
    {
        return MapL2E<CfEcefOmega>(l);
    }

    CfEcefOmegaAcc L2E(const CfLocalOmegaAcc& l) const
    {
        return MapL2E<CfEcefOmegaAcc>(l);
    }

    CfLocalVector E2L(const CfEcefVector& v) const
    {
        return MapE2L<CfLocalVector>(v);
    }

    CfLocalOmega E2L(const CfEcefOmega& v) const
    {
        return MapE2L<CfLocalOmega>(v);
    }

    CfLocalOmegaAcc E2L(const CfEcefOmegaAcc& v) const
    {
        return MapE2L<CfLocalOmegaAcc>(v);
    }

    //-- transforms of points (rotation and translation) --------------

    CfEcefPoint L2E(const CfLocalPoint& l) const
    {
        CfEcefPoint p { };
        p.x = orp.x + m.x_n * l.n + m.x_e * l.e + m.x_u * l.u;
        p.y = orp.y + m.y_n * l.n + m.y_e * l.e + m.y_u * l.u;
        p.z = orp.z + m.z_n * l.n + m.z_e * l.e + m.z_u * l.u;
        return p;
    }

    CfLocalPoint E2L(const CfEcefPoint& p) const
    {
        double dx = p.x - orp.x;
        double dy = p.y - orp.y;
        double dz = p.z - orp.z;
        CfLocalPoint l { };
        l.n = m.x_n * dx + m.y_n * dy + m.z_n * dz;
        l.e = m.x_e * dx + m.y_e * dy + m.z_e * dz;
        l.u = m.x_u * dx + m.y_u * dy + m.z_u * dz;
        return l;
    }

    //-- derived quantities ------------------------------------------

    //  Earth rotation vector in the local axes (transport rate source
    //  for Coriolis terms). Equals el.om * ( cos B, 0, sin B ).
    CfLocalOmega EarthRate() const
    {
        CfLocalOmega w { };
        w.n = el.om * cos(org.lat);
        w.e = 0.0;
        w.u = el.om * sin(org.lat);
        return w;
    }

    //  The same vector in the ECEF axes: (0, 0, el.om).
    CfEcefOmega EarthRateEcw() const
    {
        CfEcefOmega w { };
        w.z = el.om;
        return w;
    }

    //  Geodetic coordinates of a point given in the local frame.
    CfGeodetic GeoOf(const CfLocalPoint& l) const
    {
        return EcpToGeo(L2E(l), el);
    }

    //  Geodetic height of a point given in the local frame (1.4.0,
    //  EcpToAlt: GeoOf(l).alt without the latitude and longitude).
    double AltOf(const CfLocalPoint& l) const
    {
        return EcpToAlt(L2E(l), el);
    }

    constexpr const CfGeodetic& Origin() const
    {
        return org;
    }

    constexpr const CfEcefPoint& OriginEcp() const
    {
        return orp;
    }

    constexpr const GCM& Mat() const
    {
        return m;
    }

    constexpr const CfGeoid& Ellipsoid() const
    {
        return el;
    }

private:

    //  TE and TL must expose the component names of their families.
    template <class TE, class TL>
    TE MapL2E(const TL& l) const
    {
        TE v { };
        v.x = m.x_n * l.n + m.x_e * l.e + m.x_u * l.u;
        v.y = m.y_n * l.n + m.y_e * l.e + m.y_u * l.u;
        v.z = m.z_n * l.n + m.z_e * l.e + m.z_u * l.u;
        return v;
    }

    template <class TL, class TE>
    TL MapE2L(const TE& v) const
    {
        TL l { };
        l.n = m.x_n * v.x + m.y_n * v.y + m.z_n * v.z;
        l.e = m.x_e * v.x + m.y_e * v.y + m.z_e * v.z;
        l.u = m.x_u * v.x + m.y_u * v.y + m.z_u * v.z;
        return l;
    }

    CfGeoid el;
    CfGeodetic org;
    CfEcefPoint orp;
    GCM m;
};

//  Constant-initialized operator, usable as a static.
inline constexpr GEO GEO_NULL = GEO();

}   // namespace gnc

#endif  // LIB_GEODESY_H

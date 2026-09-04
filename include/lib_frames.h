//======================================================================
//  lib_frames.h
//  Frame structures with NAMED components. No indexed vectors anywhere.
//  C++17, no STL, no heap, no exceptions, no virtual, no dynamic init.
//
//  Naming rule:
//      first two letters - frame family, third letter - kind of quantity
//          C - geometric triple (free vector: direction, velocity, force)
//          P - point (radius-vector, position)
//          W - angular rate
//          M - moment
//          A - angular acceleration
//  Component names are fixed per family and never repeat across families,
//  so a value of a foreign frame cannot be passed by mistake.
//======================================================================

#ifndef LIB_FRAMES_H
#define LIB_FRAMES_H

#include <math.h>

namespace gnc
{

//----------------------------------------------------------------------
//  CfLocalVector family - local frame, components n (north), e (east), u (up).
//  Right-handed triple in the order (u, e, n): u x e = n.
//----------------------------------------------------------------------

//  Free vector in local axes (direction, velocity, acceleration, force).
struct CfLocalVector
{
    double n;
    double e;
    double u;
};

//  Point, offset from the origin of the local frame.
struct CfLocalPoint
{
    double n;
    double e;
    double u;
};

//  Angular rate projected on the local axes, rad/s.
struct CfLocalOmega
{
    double n;
    double e;
    double u;
};

//  Moment projected on the local axes, N*m.
struct CfLocalMoment
{
    double n;
    double e;
    double u;
};

//  Angular acceleration projected on the local axes, rad/s^2.
struct CfLocalOmegaAcc
{
    double n;
    double e;
    double u;
};

//----------------------------------------------------------------------
//  CfBodyVector family - body frame, origin at a fixed structural point (not CoM),
//  components bu (forward / nose), be (port), bn = bu x be.
//  At t = 0: bu == u, be == e, bn == n.
//----------------------------------------------------------------------

//  Free vector in body axes.
struct CfBodyVector
{
    double bn;
    double be;
    double bu;
};

//  Angular rate of CfBodyVector relative to CfLocalVector, projected on the body axes, rad/s.
struct CfBodyOmega
{
    double bn;
    double be;
    double bu;
};

//  Moment about the body axes, N*m.
struct CfBodyMoment
{
    double bn;
    double be;
    double bu;
};

//  Angular acceleration of CfBodyVector relative to CfLocalVector, body axes, rad/s^2.
//  The derivative of the angular rate is the same vector in both
//  frames (omega x omega = 0), so it converts as a plain vector.
struct CfBodyOmegaAcc
{
    double bn;
    double be;
    double bu;
};

//----------------------------------------------------------------------
//  ECEF family - earth-centred earth-fixed frame, components x, y, z.
//      x - equator, prime meridian
//      y - equator, 90 deg east
//      z - rotation axis, north
//  Right-handed in the order (x, y, z).
//----------------------------------------------------------------------

//  Free vector in ECEF axes.
struct CfEcefVector
{
    double x;
    double y;
    double z;
};

//  Point, radius-vector from the centre of the Earth.
struct CfEcefPoint
{
    double x;
    double y;
    double z;
};

//  Angular rate projected on the ECEF axes, rad/s.
struct CfEcefOmega
{
    double x;
    double y;
    double z;
};

//  Angular acceleration projected on the ECEF axes, rad/s^2.
struct CfEcefOmegaAcc
{
    double x;
    double y;
    double z;
};

//----------------------------------------------------------------------
//  Geodetic coordinates and reference ellipsoid.
//----------------------------------------------------------------------

struct CfGeodetic
{
    double lat;     // geodetic latitude, rad
    double lon;     // longitude, rad
    double alt;     // height above the ellipsoid, m
};

struct CfGeoid
{
    double a;       // semi-major axis, m
    double f;       // flattening
    double mu;      // geocentric gravitational constant GM, m^3/s^2
    double j2;      // second zonal harmonic
    double om;      // rotation rate of the Earth, rad/s
};

//----------------------------------------------------------------------
//  Minimal algebra on named components. Only the operations the models
//  need; every formula is written by field names, no indices.
//  Cross products are TYPED: the kind of the result follows the physics
//  (rate x arm = velocity, arm x force = moment), so a wrong combination
//  does not compile.
//----------------------------------------------------------------------

inline constexpr double PI = 3.14159265358979323846;

constexpr CfLocalVector Add(const CfLocalVector& a, const CfLocalVector& b)
{
    CfLocalVector v { };
    v.n = a.n + b.n;
    v.e = a.e + b.e;
    v.u = a.u + b.u;
    return v;
}

constexpr CfBodyVector Add(const CfBodyVector& a, const CfBodyVector& b)
{
    CfBodyVector v { };
    v.bn = a.bn + b.bn;
    v.be = a.be + b.be;
    v.bu = a.bu + b.bu;
    return v;
}

constexpr CfBodyMoment Add(const CfBodyMoment& a, const CfBodyMoment& b)
{
    CfBodyMoment v { };
    v.bn = a.bn + b.bn;
    v.be = a.be + b.be;
    v.bu = a.bu + b.bu;
    return v;
}

constexpr CfEcefVector Add(const CfEcefVector& a, const CfEcefVector& b)
{
    CfEcefVector v { };
    v.x = a.x + b.x;
    v.y = a.y + b.y;
    v.z = a.z + b.z;
    return v;
}

constexpr CfLocalVector Sub(const CfLocalVector& a, const CfLocalVector& b)
{
    CfLocalVector v { };
    v.n = a.n - b.n;
    v.e = a.e - b.e;
    v.u = a.u - b.u;
    return v;
}

constexpr CfBodyVector Sub(const CfBodyVector& a, const CfBodyVector& b)
{
    CfBodyVector v { };
    v.bn = a.bn - b.bn;
    v.be = a.be - b.be;
    v.bu = a.bu - b.bu;
    return v;
}

constexpr CfLocalVector Scale(const CfLocalVector& a, double s)
{
    CfLocalVector v { };
    v.n = a.n * s;
    v.e = a.e * s;
    v.u = a.u * s;
    return v;
}

constexpr CfBodyVector Scale(const CfBodyVector& a, double s)
{
    CfBodyVector v { };
    v.bn = a.bn * s;
    v.be = a.be * s;
    v.bu = a.bu * s;
    return v;
}

constexpr CfEcefVector Scale(const CfEcefVector& a, double s)
{
    CfEcefVector v { };
    v.x = a.x * s;
    v.y = a.y * s;
    v.z = a.z * s;
    return v;
}

constexpr double Norm2(const CfLocalVector& a)
{
    return a.n * a.n + a.e * a.e + a.u * a.u;
}

constexpr double Norm2(const CfBodyVector& a)
{
    return a.bn * a.bn + a.be * a.be + a.bu * a.bu;
}

constexpr double Norm2(const CfEcefVector& a)
{
    return a.x * a.x + a.y * a.y + a.z * a.z;
}

inline double Norm(const CfLocalVector& a)
{
    return sqrt(Norm2(a));
}

inline double Norm(const CfBodyVector& a)
{
    return sqrt(Norm2(a));
}

inline double Norm(const CfEcefVector& a)
{
    return sqrt(Norm2(a));
}

//  Velocity of a body point from rotation: v = w x r.
//  Right-handed order (bu, be, bn): bu x be = bn.
constexpr CfBodyVector OmegaCrossR(const CfBodyOmega& w, const CfBodyVector& r)
{
    CfBodyVector v { };
    v.bu = w.be * r.bn - w.bn * r.be;
    v.be = w.bn * r.bu - w.bu * r.bn;
    v.bn = w.bu * r.be - w.be * r.bu;
    return v;
}

//  Moment of a force about the origin of CfBodyVector: m = r x f.
constexpr CfBodyMoment Torque(const CfBodyVector& r, const CfBodyVector& f)
{
    CfBodyMoment m { };
    m.bu = r.be * f.bn - r.bn * f.be;
    m.be = r.bn * f.bu - r.bu * f.bn;
    m.bn = r.bu * f.be - r.be * f.bu;
    return m;
}

//----------------------------------------------------------------------
//  Constant-initialized objects. Every one of them is built by a
//  constexpr function that assigns FIELDS BY NAME, so no positional
//  aggregate initialization appears anywhere in the module and no
//  dynamic initialization is emitted at load time.
//----------------------------------------------------------------------

constexpr CfLocalVector MakeAxisN()
{
    CfLocalVector v { };
    v.n = 1.0;
    return v;
}

constexpr CfLocalVector MakeAxisE()
{
    CfLocalVector v { };
    v.e = 1.0;
    return v;
}

constexpr CfLocalVector MakeAxisU()
{
    CfLocalVector v { };
    v.u = 1.0;
    return v;
}

constexpr CfEcefVector MakeAxisX()
{
    CfEcefVector v { };
    v.x = 1.0;
    return v;
}

constexpr CfEcefVector MakeAxisY()
{
    CfEcefVector v { };
    v.y = 1.0;
    return v;
}

constexpr CfEcefVector MakeAxisZ()
{
    CfEcefVector v { };
    v.z = 1.0;
    return v;
}

//  PZ-90.11. Check the values against the governing document before use.
constexpr CfGeoid MakeEllPz9011()
{
    CfGeoid e { };
    e.a = 6378136.0;
    e.f = 1.0 / 298.25784;
    e.mu = 398600.4418e9;
    e.j2 = 1.08262575e-3;
    e.om = 7.292115e-5;
    return e;
}

//  WGS-84.
constexpr CfGeoid MakeEllWgs84()
{
    CfGeoid e { };
    e.a = 6378137.0;
    e.f = 1.0 / 298.257223563;
    e.mu = 3.986004418e14;
    e.j2 = 1.08262668355e-3;
    e.om = 7.2921159e-5;
    return e;
}

inline constexpr CfLocalVector AXIS_N = MakeAxisN();
inline constexpr CfLocalVector AXIS_E = MakeAxisE();
inline constexpr CfLocalVector AXIS_U = MakeAxisU();

inline constexpr CfEcefVector AXIS_X = MakeAxisX();
inline constexpr CfEcefVector AXIS_Y = MakeAxisY();
inline constexpr CfEcefVector AXIS_Z = MakeAxisZ();

inline constexpr CfGeoid ELL_PZ9011 = MakeEllPz9011();
inline constexpr CfGeoid ELL_WGS84 = MakeEllWgs84();

//  Earth rotation rate, rad/s (PZ-90.11).
inline constexpr double OMEGA_E = 7.292115e-5;

}   // namespace gnc

#endif  // LIB_FRAMES_H

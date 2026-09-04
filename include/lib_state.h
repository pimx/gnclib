//======================================================================
//  lib_state.h
//  VehState - state of the object at one instant of time.
//
//  Canonical storage (the single source of truth):
//      GEO  - origin of the local frame: geodetic coordinates, its
//             radius-vector in ECEF and the ECEF<->CfLocalVector matrix
//      t    - time, s
//      CfLocalPoint  - position of the object in the local frame, m
//      CfLocalVector  - velocity, m/s, and acceleration, m/s^2, local axes
//      CfRotation  - attitude operator CfBodyVector <-> CfLocalVector
//      CfBodyOmega  - angular rate, body axes, rad/s
//      CfBodyOmegaAcc  - angular acceleration (epsilon), body axes, rad/s^2
//      CfBodyVector  - arm from the CfBodyVector origin (the tail end) to the current
//             CoM, m; drifts with propellant burn, update together
//             with the mass summary (VehMassProps::CmBdy)
//
//  Every other representation is DERIVED by the typed operators on
//  request; nothing is cached, so the state cannot become internally
//  inconsistent. Getters are named by the type they return.
//
//  Kinematic conventions:
//      - the position, velocity and acceleration of the state refer
//        to the CENTRE OF MASS (dynamics and telemetry are
//        CoM-centric); arms of body points are measured from the CfBodyVector
//        origin - the TAIL end - and are reduced through the stored
//        CoM arm;
//      - the local frame is fixed to the rotating Earth, so all linear
//        quantities are relative to the rotating Earth; the ECEF
//        getters are AXIS PROJECTIONS of the same vectors, not
//        inertial quantities (no Coriolis or centrifugal terms added);
//      - the angular rate is that of CfBodyVector relative to CfLocalVector; relative to
//        ECEF it is the same vector (both frames are Earth-fixed);
//        relative to the inertial frame add GEO::EarthRate();
//      - the angular acceleration converts as a plain vector: the
//        derivative of omega is the same in both frames (w x w = 0);
//      - the position of the object in its own body frame is not a
//        meaningful quantity, therefore there is no body getter for
//        the position; positions of body POINTS are given by
//        BodyPointLop / BodyPointEcp with an explicit arm.
//
//  VehState is a literal type: a default state (zero site of GEO, identity
//  attitude, zeros) can be a constant-initialized static (STA_INIT).
//  C++17, no STL, no heap, no exceptions.
//======================================================================

#ifndef LIB_STATE_H
#define LIB_STATE_H

#include "lib_frames.h"
#include "lib_rotation.h"
#include "lib_geodesy.h"

namespace gnc
{

class VehState
{
public:

    constexpr VehState()
        : m_geo()
        , m_t(0.0)
        , m_pos()
        , m_vel()
        , m_acc()
        , m_att()
        , m_omg()
        , m_eps()
        , m_cm()
    {
    }

    //-- setters -----------------------------------------------------

    //  The site is configured outside (GEO::SetOrigin) and copied in.
    constexpr void SetGeo(const GEO& g)
    {
        m_geo = g;
    }

    constexpr void SetTime(double t)
    {
        m_t = t;
    }

    constexpr void SetPos(const CfLocalPoint& p)
    {
        m_pos = p;
    }

    constexpr void SetVel(const CfLocalVector& v)
    {
        m_vel = v;
    }

    constexpr void SetAcc(const CfLocalVector& a)
    {
        m_acc = a;
    }

    constexpr void SetAtt(const CfRotation& r)
    {
        m_att = r;
    }

    constexpr void SetOmg(const CfBodyOmega& w)
    {
        m_omg = w;
    }

    constexpr void SetEps(const CfBodyOmegaAcc& e)
    {
        m_eps = e;
    }

    //  Bulk update of the translational part of the step.
    constexpr void SetLinear(const CfLocalPoint& p, const CfLocalVector& v, const CfLocalVector& a)
    {
        m_pos = p;
        m_vel = v;
        m_acc = a;
    }

    //  Arm from the CfBodyVector origin (tail) to the current CoM; update
    //  together with the mass summary (VehMassProps::CmBdy).
    constexpr void SetCmArm(const CfBodyVector& c)
    {
        m_cm = c;
    }

    //  Bulk update of the rotational part of the step.
    constexpr void SetAngular(const CfRotation& r, const CfBodyOmega& w, const CfBodyOmegaAcc& e)
    {
        m_att = r;
        m_omg = w;
        m_eps = e;
    }

    //-- canonical storage, by reference -----------------------------

    constexpr const GEO& Geo() const
    {
        return m_geo;
    }

    constexpr double Time() const
    {
        return m_t;
    }

    constexpr const CfLocalPoint& PosLop() const
    {
        return m_pos;
    }

    constexpr const CfLocalVector& VelLoc() const
    {
        return m_vel;
    }

    constexpr const CfLocalVector& AccLoc() const
    {
        return m_acc;
    }

    constexpr const CfRotation& Att() const
    {
        return m_att;
    }

    constexpr const CfBodyOmega& OmgBdw() const
    {
        return m_omg;
    }

    constexpr const CfBodyOmegaAcc& EpsBda() const
    {
        return m_eps;
    }

    constexpr const CfBodyVector& CmArm() const
    {
        return m_cm;
    }

    //-- position: Earth frames --------------------------------------

    CfEcefPoint PosEcp() const
    {
        return m_geo.L2E(m_pos);
    }

    CfGeodetic PosGdt() const
    {
        return m_geo.GeoOf(m_pos);
    }

    //-- velocity ----------------------------------------------------

    CfBodyVector VelBdy() const
    {
        return m_att.L2B(m_vel);
    }

    CfEcefVector VelEcf() const
    {
        return m_geo.L2E(m_vel);
    }

    //-- acceleration ------------------------------------------------

    CfBodyVector AccBdy() const
    {
        return m_att.L2B(m_acc);
    }

    CfEcefVector AccEcf() const
    {
        return m_geo.L2E(m_acc);
    }

    //-- angular rate ------------------------------------------------

    CfLocalOmega OmgLow() const
    {
        return m_att.B2L(m_omg);
    }

    CfEcefOmega OmgEcw() const
    {
        return m_geo.L2E(OmgLow());
    }

    //-- angular acceleration ----------------------------------------

    CfLocalOmegaAcc EpsLoa() const
    {
        return m_att.B2L(m_eps);
    }

    CfEcefOmegaAcc EpsEca() const
    {
        return m_geo.L2E(EpsLoa());
    }

    //-- attitude in the Earth frame ---------------------------------

    //  Body vector expressed in the ECEF axes and back; together with
    //  NoseDirEcf() this gives the full attitude information relative
    //  to the Earth without introducing a separate CfBodyVector<->ECEF operator.
    CfEcefVector B2E(const CfBodyVector& b) const
    {
        return m_geo.L2E(m_att.B2L(b));
    }

    CfBodyVector E2B(const CfEcefVector& v) const
    {
        return m_att.L2B(m_geo.E2L(v));
    }

    CfLocalVector NoseDirLoc() const
    {
        return m_att.NoseDir();
    }

    CfEcefVector NoseDirEcf() const
    {
        return m_geo.L2E(m_att.NoseDir());
    }

    //-- body points -------------------------------------------------
    //  "arm" is the offset of the point from the CfBodyVector origin (the tail
    //  end), m; it is reduced to the CoM through the stored CoM arm.

    CfLocalPoint BodyPointLop(const CfBodyVector& arm) const
    {
        CfLocalVector d = m_att.B2L(Sub(arm, m_cm));
        CfLocalPoint p { };
        p.n = m_pos.n + d.n;
        p.e = m_pos.e + d.e;
        p.u = m_pos.u + d.u;
        return p;
    }

    CfEcefPoint BodyPointEcp(const CfBodyVector& arm) const
    {
        return m_geo.L2E(BodyPointLop(arm));
    }

    //  Velocity of a body point in the local axes:
    //  v_cm + R ( w x (arm - cm) ).
    CfLocalVector BodyPointVelLoc(const CfBodyVector& arm) const
    {
        CfLocalVector d = m_att.B2L(OmegaCrossR(m_omg, Sub(arm, m_cm)));
        return Add(m_vel, d);
    }

    //  The CfBodyVector origin (tail) itself: BodyPointLop of a zero arm.
    CfLocalPoint TailLop() const
    {
        CfBodyVector z { };
        return BodyPointLop(z);
    }

private:

    GEO m_geo;
    double m_t;
    CfLocalPoint m_pos;
    CfLocalVector m_vel;
    CfLocalVector m_acc;
    CfRotation m_att;
    CfBodyOmega m_omg;
    CfBodyOmegaAcc m_eps;
    CfBodyVector m_cm;
};

//  Constant-initialized default state: GEO null site, zero kinematics,
//  identity attitude.
inline constexpr VehState STA_INIT = VehState();

}   // namespace gnc

#endif  // LIB_STATE_H

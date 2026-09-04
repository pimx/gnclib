//======================================================================
//  lib_dynamics.h
//  Rigid-body integration of the VehState state over a step dt under an
//  external force F and moment M.
//
//  INR - mass properties: mass and the symmetric inertia TENSOR about
//  the CoM in the body axes. The off-diagonal fields are TENSOR
//  elements J_ab (equal to MINUS the centrifugal moments of inertia
//  I_ab of the GOST convention). Mass is constant within a step; for
//  propellant flow update INR between steps.
//
//  Integrate(state, dt, F, M, inr) -> new state:
//      - translation: F is the TOTAL external force on the CoM
//        (including gravity - nothing is added implicitly), constant
//        in the LOCAL axes over the step; with a = F/m constant the
//        integration is EXACT: v' = v + a dt, p' = p + v dt + a dt^2/2;
//      - rotation: Euler's equations J eps = M - w x (J w) with M
//        constant in the BODY axes, integrated by RK4 over the step
//        (handles the gyroscopic term); the attitude is advanced by
//        the exact exponential map of the trapezoidally averaged rate
//        0.5 (w + w') - a second-order scheme with an exactly unit
//        quaternion. Subdivide the step externally for higher accuracy;
//      - the local frame is taken quasi-inertial for the rotational
//        dynamics (Earth-rate terms ~1e-4 rad/s are neglected); the
//        translational Coriolis/centrifugal terms are NOT added - put
//        them into F if the model requires them;
//      - the returned state carries the accelerations of ITS time:
//        AccLoc = F/m, EpsBda = J^-1 (M - w' x J w') at the new rate.
//
//  The moment M is taken ABOUT THE CoM. A force applied at a body
//  point (thrust at the gimbal, leg contact) is reduced to the CoM by
//  MomentAboutCm(): the arms are measured from the CfBodyVector origin - the
//  TAIL end - and the current CoM arm comes from the state
//  (VehState::CmArm) or the mass model (VehMassProps::CmBdy).
//
//  An overload takes the force in the BODY axes (thrust): it is
//  converted to the local axes at the INITIAL attitude, so it is
//  treated as constant in the local axes within the step, and it is
//  assumed to act THROUGH the CoM (produces no moment); the error is
//  O(|w| dt) - subdivide the step when the vehicle turns fast.
//
//  C++17, no STL, no heap, no exceptions.
//======================================================================

#ifndef LIB_DYNAMICS_H
#define LIB_DYNAMICS_H

#include "lib_frames.h"
#include "lib_rotation.h"
#include "lib_state.h"

namespace gnc
{

//----------------------------------------------------------------------
//  Mass properties. Inertia tensor about the CoM, body axes, kg*m^2:
//      row bn: [ j_bnbn, j_bnbe, j_bnbu ]
//      row be: [ j_bnbe, j_bebe, j_bebu ]
//      row bu: [ j_bnbu, j_bebu, j_bubu ]
//----------------------------------------------------------------------
struct INR
{
    double m;           // mass, kg
    double j_bnbn;
    double j_bebe;
    double j_bubu;
    double j_bnbe;      // tensor element (= -I_bnbe)
    double j_bnbu;      // tensor element (= -I_bnbu)
    double j_bebu;      // tensor element (= -I_bebu)
};

//  Angular acceleration from Euler's equations:
//      eps = J^-1 ( M - w x (J w) )
//  Closed-form symmetric 3x3 solve (cofactors). Returns zero on a
//  degenerate tensor (|det| below the guard).
inline CfBodyOmegaAcc EpsOf(const INR& j, const CfBodyMoment& mom, const CfBodyOmega& w)
{
    //  Angular momentum h = J w, by field names.
    double hn = j.j_bnbn * w.bn + j.j_bnbe * w.be + j.j_bnbu * w.bu;
    double he = j.j_bnbe * w.bn + j.j_bebe * w.be + j.j_bebu * w.bu;
    double hu = j.j_bnbu * w.bn + j.j_bebu * w.be + j.j_bubu * w.bu;

    //  Gyroscopic term c = w x h (right-handed order bu, be, bn).
    double cu = w.be * hn - w.bn * he;
    double ce = w.bn * hu - w.bu * hn;
    double cn = w.bu * he - w.be * hu;

    double rn = mom.bn - cn;
    double re = mom.be - ce;
    double ru = mom.bu - cu;

    //  Symmetric inverse via cofactors; A..F name the tensor blocks.
    double a = j.j_bnbn;
    double b = j.j_bebe;
    double c = j.j_bubu;
    double f = j.j_bnbe;
    double e = j.j_bnbu;
    double d = j.j_bebu;
    double k11 = b * c - d * d;
    double k12 = e * d - f * c;
    double k13 = f * d - e * b;
    double k22 = a * c - e * e;
    double k23 = e * f - a * d;
    double k33 = a * b - f * f;
    double det = a * k11 + f * k12 + e * k13;

    CfBodyOmegaAcc eps { };
    if (det < 1.0e-30 && det > -1.0e-30)
    {
        return eps;
    }
    double inv = 1.0 / det;
    eps.bn = (k11 * rn + k12 * re + k13 * ru) * inv;
    eps.be = (k12 * rn + k22 * re + k23 * ru) * inv;
    eps.bu = (k13 * rn + k23 * re + k33 * ru) * inv;
    return eps;
}

//  Moment about the CoM of a force f applied at a body point:
//  "atArm" and "cmArm" are measured from the CfBodyVector origin (the tail).
constexpr CfBodyMoment MomentAboutCm(const CfBodyVector& f, const CfBodyVector& atArm, const CfBodyVector& cmArm)
{
    return Torque(Sub(atArm, cmArm), f);
}

//  w + eps * s : the legal mixing of kinds inside an integrator step.
constexpr CfBodyOmega AdvanceRate(const CfBodyOmega& w, const CfBodyOmegaAcc& e, double s)
{
    CfBodyOmega r { };
    r.bn = w.bn + e.bn * s;
    r.be = w.be + e.be * s;
    r.bu = w.bu + e.bu * s;
    return r;
}

//  One integration step. F in the LOCAL axes (total, incl. gravity),
//  M about the CoM in the BODY axes; both constant over the step.
inline VehState Integrate(const VehState& s, double dt, const CfLocalVector& f, const CfBodyMoment& mom,
                     const INR& inr)
{
    VehState r = s;
    r.SetTime(s.Time() + dt);

    //-- translation: exact for constant acceleration -----------------
    double m = inr.m;
    if (m < 1.0e-9)
    {
        m = 1.0e-9;
    }
    CfLocalVector a { };
    a.n = f.n / m;
    a.e = f.e / m;
    a.u = f.u / m;
    const CfLocalVector& v0 = s.VelLoc();
    CfLocalVector v1 { };
    v1.n = v0.n + a.n * dt;
    v1.e = v0.e + a.e * dt;
    v1.u = v0.u + a.u * dt;
    const CfLocalPoint& p0 = s.PosLop();
    double hdt2 = 0.5 * dt * dt;
    CfLocalPoint p1 { };
    p1.n = p0.n + v0.n * dt + a.n * hdt2;
    p1.e = p0.e + v0.e * dt + a.e * hdt2;
    p1.u = p0.u + v0.u * dt + a.u * hdt2;
    r.SetLinear(p1, v1, a);

    //-- rotation: RK4 on Euler's equations ---------------------------
    const CfBodyOmega& w0 = s.OmgBdw();
    CfBodyOmegaAcc k1 = EpsOf(inr, mom, w0);
    CfBodyOmegaAcc k2 = EpsOf(inr, mom, AdvanceRate(w0, k1, 0.5 * dt));
    CfBodyOmegaAcc k3 = EpsOf(inr, mom, AdvanceRate(w0, k2, 0.5 * dt));
    CfBodyOmegaAcc k4 = EpsOf(inr, mom, AdvanceRate(w0, k3, dt));
    CfBodyOmega w1 { };
    double s6 = dt / 6.0;
    w1.bn = w0.bn + s6 * (k1.bn + 2.0 * k2.bn + 2.0 * k3.bn + k4.bn);
    w1.be = w0.be + s6 * (k1.be + 2.0 * k2.be + 2.0 * k3.be + k4.be);
    w1.bu = w0.bu + s6 * (k1.bu + 2.0 * k2.bu + 2.0 * k3.bu + k4.bu);

    //  Attitude: exponential map of the trapezoidally averaged rate.
    CfBodyOmega wm { };
    wm.bn = 0.5 * (w0.bn + w1.bn);
    wm.be = 0.5 * (w0.be + w1.be);
    wm.bu = 0.5 * (w0.bu + w1.bu);
    CfRotation att = s.Att();
    att.Propagate(wm, dt);

    //  Stored accelerations belong to the NEW time.
    r.SetAngular(att, w1, EpsOf(inr, mom, w1));
    return r;
}

//  The same step with the force given in the BODY axes (thrust).
//  Converted at the initial attitude, i.e. treated as constant in the
//  local axes within the step; error O(|w| dt).
inline VehState Integrate(const VehState& s, double dt, const CfBodyVector& fBdy,
                     const CfBodyMoment& mom, const INR& inr)
{
    return Integrate(s, dt, s.Att().B2L(fBdy), mom, inr);
}

}   // namespace gnc

#endif  // LIB_DYNAMICS_H

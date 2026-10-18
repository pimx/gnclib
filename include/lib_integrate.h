//============================================================================
// lib_integrate.h -- RK4 integration of the rigid-body motion with the
// forces re-evaluated at every stage, plus auxiliary states (actuators,
// slosh, bending, ...) integrated by the same RK4 (gnclib 1.1.0).
//
// The vehicle supplies one derivative function:
//     fn(ctx, s, aux, &d): for the state s (position, velocity, attitude,
//     body rate, time) and the auxiliary states aux it returns the total
//     force in the LOCAL axes d.f (incl. gravity, Coriolis if wanted), the
//     moment about the CoM in the BODY axes d.mom, the mass properties
//     d.inr at that state and the derivatives of the auxiliary states
//     d.xdot[0 .. aux.n - 1].
// Stages (classic RK4): translation and body rate by RK4; the attitude by
// RKMK4 (Munthe-Kaas, 1.7.1): the rotation vector th of the step obeys the
// Bortz equation th' = w + th x w / 2 + th x (th x w) / 12 (body rate w,
// right multiplication); its RK4 slopes k_i = dexp^-1(th_i, w_i) carry the
// non-commutativity (coning) term, the stages k = 2, 3, 4 start from the
// initial attitude turned by th_i = c_i dt k_(i-1), the final attitude --
// the exponential map of dt (k1 + 2 k2 + 2 k3 + k4) / 6. Up to 1.7.0 the
// slopes were the bare rates w_i: 2nd order in the attitude for a rate
// that changes its direction (bitwise the same for a rotation about one
// axis). Stored accelerations are the RK4 means over the step.
// The existing Integrate (lib_dynamics.h, forces constant over the step)
// stays as the light option. C++17, no STL, no heap.
//============================================================================
#ifndef GNCLIB_LIB_INTEGRATE_H
#define GNCLIB_LIB_INTEGRATE_H

#include "lib_dynamics.h"

namespace gnc
{

// Auxiliary states integrated together with the rigid body.
struct AuxState
{
    static const int NMAX = 32;
    int n;
    double x[NMAX];
};

// Result of the derivative function at one stage.
struct DynDeriv
{
    CfLocalVector f;              // total force, LOCAL axes, N
    CfBodyMoment mom;             // moment about the CoM, BODY axes, N*m
    INR inr;                      // mass properties at the stage
    double xdot[AuxState::NMAX];  // derivatives of the auxiliary states
};

typedef void (*DynFn)(void* ctx, const VehState& s, const AuxState& aux, DynDeriv* d);

namespace detail
{
// Slope of the rotation vector at a stage (Bortz, dexp^-1 to 2nd order --
// enough for RKMK4): k = w + th x w / 2 + th x (th x w) / 12, th = kPrev tau.
//   kPrev -- slope of the previous stage, rad/s; tau -- stage time, s;
//   w -- body rate at the stage, rad/s
inline CfBodyOmega RotSlope(const CfBodyOmega& kPrev, double tau, const CfBodyOmega& w)
{
    double tn = kPrev.bn * tau;
    double te = kPrev.be * tau;
    double tu = kPrev.bu * tau;
    // a = th x w, right-handed order (bu, be, bn)
    double au = te * w.bn - tn * w.be;
    double ae = tn * w.bu - tu * w.bn;
    double an = tu * w.be - te * w.bu;
    // b = th x a
    double bu = te * an - tn * ae;
    double be = tn * au - tu * an;
    double bn = tu * ae - te * au;
    CfBodyOmega k { };
    k.bn = w.bn + 0.5 * an + bn / 12.0;
    k.be = w.be + 0.5 * ae + be / 12.0;
    k.bu = w.bu + 0.5 * au + bu / 12.0;
    return k;
}

inline CfLocalVector LinAcc(const DynDeriv& d)
{
    double m = (d.inr.m < 1.0e-9) ? 1.0e-9 : d.inr.m;
    CfLocalVector a { };
    a.n = d.f.n / m;
    a.e = d.f.e / m;
    a.u = d.f.u / m;
    return a;
}

// stage state: from s0 advanced by the stage increments
inline VehState Stage(const VehState& s0, double tau, const CfLocalVector& dv,
                      const CfLocalVector& vStage, const CfBodyOmega& wStage,
                      const CfBodyOmega& wProp, const CfLocalVector& aPrev,
                      const CfBodyOmegaAcc& ePrev)
{
    VehState r = s0;
    r.SetTime(s0.Time() + tau);
    const CfLocalPoint& p0 = s0.PosLop();
    CfLocalPoint p { };
    p.n = p0.n + vStage.n * tau;
    p.e = p0.e + vStage.e * tau;
    p.u = p0.u + vStage.u * tau;
    CfLocalVector v { };
    v.n = s0.VelLoc().n + dv.n;
    v.e = s0.VelLoc().e + dv.e;
    v.u = s0.VelLoc().u + dv.u;
    r.SetLinear(p, v, aPrev);
    CfRotation att = s0.Att();
    att.Propagate(wProp, tau);
    r.SetAngular(att, wStage, ePrev);
    return r;
}
}   // namespace detail

//  One RK4 step of dt with the stage forces from fn; aux (may be 0) is
//  advanced in place. Returns the new state.
inline VehState IntegrateRK4(const VehState& s, AuxState* aux, double dt, DynFn fn, void* ctx)
{
    AuxState a0 { };
    if (aux != 0)
    {
        a0 = *aux;
    }
    int na = (a0.n < 0) ? 0 : ((a0.n > AuxState::NMAX) ? AuxState::NMAX : a0.n);
    AuxState ak = a0;
    DynDeriv d1 { };
    DynDeriv d2 { };
    DynDeriv d3 { };
    DynDeriv d4 { };
    const CfLocalVector& v0 = s.VelLoc();
    const CfBodyOmega& w0 = s.OmgBdw();
    int i = 0;
    // stage 1
    fn(ctx, s, a0, &d1);
    CfLocalVector a1 = detail::LinAcc(d1);
    CfBodyOmegaAcc e1 = EpsOf(d1.inr, d1.mom, w0);
    // stage 2 (dt / 2)
    double h = 0.5 * dt;
    CfLocalVector dv2 { a1.n * h, a1.e * h, a1.u * h };
    CfLocalVector v2 { v0.n + dv2.n, v0.e + dv2.e, v0.u + dv2.u };
    CfBodyOmega w2 = AdvanceRate(w0, e1, h);
    for (i = 0; i < na; ++i)
    {
        ak.x[i] = a0.x[i] + d1.xdot[i] * h;
    }
    VehState s2 = detail::Stage(s, h, dv2, v0, w2, w0, a1, e1);
    fn(ctx, s2, ak, &d2);
    CfLocalVector a2 = detail::LinAcc(d2);
    CfBodyOmegaAcc e2 = EpsOf(d2.inr, d2.mom, w2);
    CfBodyOmega k2 = detail::RotSlope(w0, h, w2);
    // stage 3 (dt / 2)
    CfLocalVector dv3 { a2.n * h, a2.e * h, a2.u * h };
    CfBodyOmega w3 = AdvanceRate(w0, e2, h);
    for (i = 0; i < na; ++i)
    {
        ak.x[i] = a0.x[i] + d2.xdot[i] * h;
    }
    VehState s3 = detail::Stage(s, h, dv3, v2, w3, k2, a2, e2);
    fn(ctx, s3, ak, &d3);
    CfLocalVector a3 = detail::LinAcc(d3);
    CfBodyOmegaAcc e3 = EpsOf(d3.inr, d3.mom, w3);
    CfBodyOmega k3 = detail::RotSlope(k2, h, w3);
    // stage 4 (dt)
    CfLocalVector dv4 { a3.n * dt, a3.e * dt, a3.u * dt };
    CfLocalVector v3 { v0.n + dv3.n, v0.e + dv3.e, v0.u + dv3.u };
    CfBodyOmega w4 = AdvanceRate(w0, e3, dt);
    for (i = 0; i < na; ++i)
    {
        ak.x[i] = a0.x[i] + d3.xdot[i] * dt;
    }
    VehState s4 = detail::Stage(s, dt, dv4, v3, w4, k3, a3, e3);
    fn(ctx, s4, ak, &d4);
    CfLocalVector a4 = detail::LinAcc(d4);
    CfBodyOmegaAcc e4 = EpsOf(d4.inr, d4.mom, w4);
    CfBodyOmega k4 = detail::RotSlope(k3, dt, w4);
    CfLocalVector v4 { v0.n + dv4.n, v0.e + dv4.e, v0.u + dv4.u };
    // combine
    double s6 = dt / 6.0;
    VehState r = s;
    r.SetTime(s.Time() + dt);
    const CfLocalPoint& p0 = s.PosLop();
    CfLocalPoint p1 { };
    p1.n = p0.n + s6 * (v0.n + 2.0 * v2.n + 2.0 * v3.n + v4.n);
    p1.e = p0.e + s6 * (v0.e + 2.0 * v2.e + 2.0 * v3.e + v4.e);
    p1.u = p0.u + s6 * (v0.u + 2.0 * v2.u + 2.0 * v3.u + v4.u);
    CfLocalVector v1 { };
    v1.n = v0.n + s6 * (a1.n + 2.0 * a2.n + 2.0 * a3.n + a4.n);
    v1.e = v0.e + s6 * (a1.e + 2.0 * a2.e + 2.0 * a3.e + a4.e);
    v1.u = v0.u + s6 * (a1.u + 2.0 * a2.u + 2.0 * a3.u + a4.u);
    CfLocalVector aEnd { };
    aEnd.n = (a1.n + 2.0 * a2.n + 2.0 * a3.n + a4.n) / 6.0;
    aEnd.e = (a1.e + 2.0 * a2.e + 2.0 * a3.e + a4.e) / 6.0;
    aEnd.u = (a1.u + 2.0 * a2.u + 2.0 * a3.u + a4.u) / 6.0;
    r.SetLinear(p1, v1, aEnd);
    CfBodyOmega w1 { };
    w1.bn = w0.bn + s6 * (e1.bn + 2.0 * e2.bn + 2.0 * e3.bn + e4.bn);
    w1.be = w0.be + s6 * (e1.be + 2.0 * e2.be + 2.0 * e3.be + e4.be);
    w1.bu = w0.bu + s6 * (e1.bu + 2.0 * e2.bu + 2.0 * e3.bu + e4.bu);
    CfBodyOmega km { };
    km.bn = (w0.bn + 2.0 * k2.bn + 2.0 * k3.bn + k4.bn) / 6.0;
    km.be = (w0.be + 2.0 * k2.be + 2.0 * k3.be + k4.be) / 6.0;
    km.bu = (w0.bu + 2.0 * k2.bu + 2.0 * k3.bu + k4.bu) / 6.0;
    CfRotation att = s.Att();
    att.Propagate(km, dt);
    CfBodyOmegaAcc eEnd { };
    eEnd.bn = (e1.bn + 2.0 * e2.bn + 2.0 * e3.bn + e4.bn) / 6.0;
    eEnd.be = (e1.be + 2.0 * e2.be + 2.0 * e3.be + e4.be) / 6.0;
    eEnd.bu = (e1.bu + 2.0 * e2.bu + 2.0 * e3.bu + e4.bu) / 6.0;
    r.SetAngular(att, w1, eEnd);
    if (aux != 0)
    {
        for (i = 0; i < na; ++i)
        {
            aux->x[i] = a0.x[i] + s6 * (d1.xdot[i] + 2.0 * d2.xdot[i] + 2.0 * d3.xdot[i]
                                        + d4.xdot[i]);
        }
    }
    return r;
}

}   // namespace gnc

#endif // GNCLIB_LIB_INTEGRATE_H

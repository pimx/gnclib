//============================================================================
// lib_imu.h -- ideal strapdown IMU readings for a site-fixed (rotating
// Earth) frame (gnclib 1.1.0). Noise and biases are
// added by the caller.
//   gyro  = body rate relative to the Earth + Earth rate (body axes);
//   accel = a_rel + 2 Omega x v - g(r) (body axes) + eps x r + w x (w x r)
//           for the IMU at the lever r from the CoM (body axes).
//   geo = 0: flat non-rotating Earth with the constant gravity g0.
// On the ground (grounded) the body rate and the lever terms are zero.
// C++17, no STL, no heap.
//============================================================================
#ifndef GNCLIB_LIB_IMU_H
#define GNCLIB_LIB_IMU_H

#include "lib_earth.h"
#include "lib_state.h"
#include "lib_frame_ops.h"

namespace gnc
{

inline void ImuIdeal(const VehState& st, const CfLocalVector& accRel, const GEO* geo, double g0,
                     const CfBodyVector& lever, bool grounded, CfBodyOmega* gyro,
                     CfBodyVector* accel)
{
    CfLocalVector gL = ::mkLoc(0.0, 0.0, -g0);
    CfLocalVector omL = ::mkLoc(0.0, 0.0, 0.0);
    if (geo != 0)
    {
        gL = EarthGravityLoc(*geo, st.PosLop());
        omL = EarthRateLoc(*geo);
    }
    CfLocalVector corA = ::scale(::cross(omL, st.VelLoc()), 2.0);
    CfBodyVector f = st.Att().L2B(::sub(::add(accRel, corA), gL));
    if (!grounded)
    {
        CfBodyVector wv = ::omgVec(st.OmgBdw());
        CfBodyVector ev = ::mkBody(st.EpsBda().bu, st.EpsBda().be, st.EpsBda().bn);
        f = ::add(f, ::add(::cross(ev, lever), ::cross(wv, ::cross(wv, lever))));
    }
    CfBodyOmega w = grounded ? ::mkOmg(0.0, 0.0, 0.0) : st.OmgBdw();
    *gyro = ::omgOf(::add(::omgVec(w), st.Att().L2B(omL)));
    *accel = f;
}

}   // namespace gnc

#endif // GNCLIB_LIB_IMU_H

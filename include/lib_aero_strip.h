//============================================================================
// lib_aero_strip.h -- strip aerodynamics of a cylinder with the Mach
// dependence of both drag coefficients and the retro-plume shielding by
// the thrust coefficient C_T (gnclib 1.1.0; from HOP90 V90.33..V90.37).
// Complements EnvAeroforces::CylinderStrip (lib_aerodynamics.h).
//   MachFactorAxialBlunt(M) -- axial Cd factor of a blunt body with the
//     transonic rise: M 0 / 0.8 / 1.0 / 1.2 / 1.5 / 2.0 / 3.0 ->
//     1.00 / 1.06 / 1.44 / 1.56 / 1.38 / 1.25 / 1.13;
//   MachFactorCrossflow(Mc) -- crossflow Cd factor by the crossflow Mach
//     number: Mc 0 / 0.6 / 0.9 / 1.2 / 2.0 -> 1.0 / 1.0 / 1.3 / 1.5 / 1.3;
//   CxRetroPlumeCT(cx, k, ct) -- descending tail first with thrust:
//     cx / (1 + k C_T), C_T = T / (q A_end);
//   CylinderStripCT -- the body split into nStrips along its axis, the
//     crossflow of each strip by its own velocity v + w x r (centre of
//     pressure and rotation damping follow), axial drag on the axis.
// C++17, no STL, no heap.
//============================================================================
#ifndef GNCLIB_LIB_AERO_STRIP_H
#define GNCLIB_LIB_AERO_STRIP_H

#include "lib_frame_ops.h"
#include "lib_dynamics.h"

namespace gnc
{

inline double AeroTable(const double* xs, const double* ys, int n, double x)
{
    if (x <= xs[0])
    {
        return ys[0];
    }
    int i = 0;
    for (i = 1; i < n; ++i)
    {
        if (x <= xs[i])
        {
            double k = (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
            return ys[i - 1] + k * (ys[i] - ys[i - 1]);
        }
    }
    return ys[n - 1];
}

inline double MachFactorAxialBlunt(double mach)
{
    static const double xs[] = { 0.0, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0 };
    static const double ys[] = { 1.00, 1.06, 1.44, 1.56, 1.38, 1.25, 1.13 };
    return AeroTable(xs, ys, 7, mach);
}

inline double MachFactorCrossflow(double machCross)
{
    static const double xs[] = { 0.0, 0.6, 0.9, 1.2, 2.0 };
    static const double ys[] = { 1.00, 1.00, 1.30, 1.50, 1.30 };
    return AeroTable(xs, ys, 5, machCross);
}

//  vAx -- axial air velocity (< 0: flow onto the tail); q -- dynamic
//  pressure; aEnd -- end area; thrust -- total thrust, N
inline double CxRetroPlumeCT(double cx, double k, double vAx, double q, double aEnd,
                             double thrust)
{
    if (vAx < 0.0 && thrust > 0.0 && q > 1.0)
    {
        double ct = thrust / (q * aEnd);
        return cx / (1.0 + k * ct);
    }
    return cx;
}

//  vRel -- velocity relative to the air, BODY axes; omg -- body rate;
//  rho -- density; sound -- speed of sound; thrust -- total thrust (plume);
//  length, radius -- cylinder; cmFromTail -- CoM above the tail end;
//  cdLateral, cdAxial -- subsonic coefficients; retroPlumeK -- shielding
//  (0 off); nStrips -- 1..32. Force and moment about the CoM, BODY axes.
inline void CylinderStripCT(const CfBodyVector& vRel, const CfBodyOmega& omg, double rho,
                            double sound, double thrust, double length, double radius,
                            double cmFromTail, double cdLateral, double cdAxial,
                            double retroPlumeK, int nStrips, CfBodyVector* force,
                            CfBodyMoment* moment)
{
    *force = ::mkBody(0.0, 0.0, 0.0);
    *moment = ::mkMom(0.0, 0.0, 0.0);
    if (rho <= 0.0 || length <= 0.0 || radius <= 0.0 || sound <= 0.0)
    {
        return;
    }
    int ns = (nStrips < 1) ? 1 : ((nStrips > 32) ? 32 : nStrips);
    double q = 0.5 * rho;
    double V = sqrt(vRel.bu * vRel.bu + vRel.be * vRel.be + vRel.bn * vRel.bn);
    double dz = length / ns;
    double aEnd = PI * radius * radius;
    double cax = CxRetroPlumeCT(cdAxial * MachFactorAxialBlunt(V / sound), retroPlumeK,
                                vRel.bu, q, aEnd, thrust);
    CfBodyVector f = ::mkBody(-q * cax * PI * radius * radius * vRel.bu * fabs(vRel.bu), 0.0,
                            0.0);
    CfBodyVector mom = ::mkBody(0.0, 0.0, 0.0);
    CfBodyVector wv = ::omgVec(omg);
    int k = 0;
    for (k = 0; k < ns; ++k)
    {
        CfBodyVector r = ::mkBody((k + 0.5) * dz - cmFromTail, 0.0, 0.0);
        CfBodyVector vk = ::add(vRel, ::cross(wv, r));
        double vl = sqrt(vk.be * vk.be + vk.bn * vk.bn);
        double cl = cdLateral * MachFactorCrossflow(vl / sound) * 2.0 * radius * dz;
        CfBodyVector dF = ::mkBody(0.0, -q * cl * vk.be * vl, -q * cl * vk.bn * vl);
        f = ::add(f, dF);
        mom = ::add(mom, ::momVec(Torque(r, dF)));
    }
    *force = f;
    *moment = ::momOf(mom);
}

}   // namespace gnc

#endif // GNCLIB_LIB_AERO_STRIP_H

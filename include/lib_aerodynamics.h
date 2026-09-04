//======================================================================
//  lib_aerodynamics.h
//  EnvAeroforces - aerodynamic force models for all flight regimes of the vehicle
//  line: constant-Cd drag, landing drag law Cx(M) with retro-plume
//  shielding, strip model of a cylinder, fin authority, ground effect.
//
//  Frames: drag is COLLINEAR with the air velocity, so it is given as
//  two overloads - in body axes (CfBodyVector) and in local axes (CfLocalVector); the
//  result type equals the argument type, frames cannot be mixed.
//  The strip model consumes CfBodyOmega (angular rate) and produces CfBodyVector (force)
//  and CfBodyMoment (moment) - the kinds are enforced by the compiler.
//
//  EnvAeroforces keeps a reference to EnvAtmosphere (owned by the caller); one instance
//  serves any number of bodies. C++17, no STL, no heap.
//======================================================================

#ifndef LIB_AERODYNAMICS_H
#define LIB_AERODYNAMICS_H

#include <math.h>

#include "lib_frames.h"
#include "lib_atmosphere.h"

namespace gnc
{

class EnvAeroforces
{
public:

    constexpr explicit EnvAeroforces(const EnvAtmosphere& atm)
        : m_atm(atm)
    {
    }

    //-- constant-Cd drag --------------------------------------------

    //  Magnitude of the constant-Cd drag, N, at air speed v (scalar
    //  kernel of the overloads below; for point predictions).
    double DragMagnitude(double h, double v, double cd, double aRef) const
    {
        double rho = m_atm.Density(h);
        if (v < 1.0e-6 || rho <= 0.0)
        {
            return 0.0;
        }
        return 0.5 * rho * v * v * cd * aRef;
    }

    //  F = -1/2 rho |v| v Cd A, against the air velocity.
    CfBodyVector DragConstCd(double h, const CfBodyVector& vAir, double cd, double aRef) const
    {
        double v = Norm(vAir);
        double f = DragMagnitude(h, v, cd, aRef);
        CfBodyVector z { };
        if (f <= 0.0)
        {
            return z;
        }
        return Scale(vAir, -f / v);
    }

    CfLocalVector DragConstCd(double h, const CfLocalVector& vAir, double cd, double aRef) const
    {
        double v = Norm(vAir);
        double f = DragMagnitude(h, v, cd, aRef);
        CfLocalVector z { };
        if (f <= 0.0)
        {
            return z;
        }
        return Scale(vAir, -f / v);
    }

    //-- landing drag law --------------------------------------------

    //  Landing Cx(M) of a blunt base: tabulated (transonic hump on a
    //  supersonic plateau), linear interpolation, nearest node outside
    //  the table. The table is a discretization of the former analytic
    //  fit and is to be replaced by tunnel data without code changes.
    static double CxMachLanding(double mach)
    {
        static constexpr double tabM[] =
        {
            0.0, 0.3, 0.6, 0.8, 0.9, 1.0,
            1.1, 1.2, 1.4, 1.6, 2.0, 2.5,
            3.0, 4.0, 5.0, 6.0, 10.0
        };
        static constexpr double tabCx[] =
        {
            1.5212, 1.5102, 1.5506, 1.6460, 1.6870, 1.6964,
            1.6653, 1.6029, 1.4659, 1.3886, 1.3332, 1.2940,
            1.2715, 1.2546, 1.2509, 1.2502, 1.2500
        };
        const int n = (int)(sizeof(tabM) / sizeof(tabM[0]));
        if (mach <= tabM[0])
        {
            return tabCx[0];
        }
        if (mach >= tabM[n - 1])
        {
            return tabCx[n - 1];
        }
        int i = 0;
        while (i < n - 2 && mach > tabM[i + 1])
        {
            i = i + 1;
        }
        double f = (mach - tabM[i]) / (tabM[i + 1] - tabM[i]);
        return tabCx[i] + (tabCx[i + 1] - tabCx[i]) * f;
    }

    //  Ascent Cx(M) of a cone-cylinder with a nose fairing: tabulated
    //  (subsonic plateau, transonic hump, supersonic decay), linear
    //  interpolation, nearest node outside the table. A generic curve
    //  to be replaced by tunnel data without code changes.
    static double CxMachAscent(double mach)
    {
        static constexpr double tabM[] =
        {
            0.0, 0.5, 0.8, 0.9, 1.0, 1.05,
            1.1, 1.2, 1.4, 1.8, 2.5, 3.5,
            5.0, 8.0
        };
        static constexpr double tabCx[] =
        {
            0.280, 0.270, 0.300, 0.360, 0.480, 0.550,
            0.560, 0.540, 0.480, 0.400, 0.330, 0.290,
            0.260, 0.250
        };
        const int n = (int)(sizeof(tabM) / sizeof(tabM[0]));
        if (mach <= tabM[0])
        {
            return tabCx[0];
        }
        if (mach >= tabM[n - 1])
        {
            return tabCx[n - 1];
        }
        int i = 0;
        while (i < n - 2 && mach > tabM[i + 1])
        {
            i = i + 1;
        }
        double f = (mach - tabM[i]) / (tabM[i + 1] - tabM[i]);
        return tabCx[i] + (tabCx[i + 1] - tabCx[i]) * f;
    }

    //  Relative Mach correction of the axial Cx: 1.0 below M = 0.3
    //  (incompressible), above - growth by CxMachLanding normalized to
    //  the subsonic value.
    static double MachFactorAxial(double mach)
    {
        if (mach <= 0.3)
        {
            return 1.0;
        }
        double f = CxMachLanding(mach) / CxMachLanding(0.3);
        if (f > 1.0)
        {
            return f;
        }
        return 1.0;
    }

    //  Retro-plume shielding: Cx_eff = Cx (1 - k thrFrac), floor 5% of
    //  Cx keeps the integration well-posed; k <= 0 or engine off - Cx.
    static double CxRetroPlume(double cx, double k, double thrFrac)
    {
        if (k <= 0.0 || thrFrac <= 0.0)
        {
            return cx;
        }
        if (thrFrac > 1.0)
        {
            thrFrac = 1.0;
        }
        double f = 1.0 - k * thrFrac;
        if (f < 0.05)
        {
            f = 0.05;
        }
        return cx * f;
    }

    //  Prandtl-Glauert compressibility correction for subsonic lift
    //  coefficients: 1 / sqrt(1 - M^2), M clamped to 0.95.
    static double PrandtlGlauert(double mach)
    {
        if (mach < 0.0)
        {
            mach = 0.0;
        }
        if (mach > 0.95)
        {
            mach = 0.95;
        }
        return 1.0 / sqrt(1.0 - mach * mach);
    }

    //  Cheeseman-Bennett ground effect: T_ige / T_oge =
    //  1 / (1 - (R/(4h))^2); clamped by maxFactor; maxFactor at
    //  h <= R/4; diskRadius <= 0 disables the effect (1.0).
    static double GroundEffect(double diskRadius, double h, double maxFactor)
    {
        if (diskRadius <= 0.0)
        {
            return 1.0;
        }
        if (maxFactor < 1.0)
        {
            maxFactor = 1.0;
        }
        double quarter = 0.25 * diskRadius;
        if (h <= quarter)
        {
            return maxFactor;
        }
        double r = quarter / h;
        double f = 1.0 / (1.0 - r * r);
        if (f > maxFactor)
        {
            f = maxFactor;
        }
        if (f < 1.0)
        {
            f = 1.0;
        }
        return f;
    }

    //  Magnitude of the landing drag, N: Cx(M) with retro-plume
    //  reduction (scalar kernel of the overloads below).
    double DragLandingMagnitude(double h, double v, double aRef,
                                double retroPlumeK, double thrFrac) const
    {
        double rho = m_atm.Density(h);
        if (v < 1.0e-6 || rho <= 0.0)
        {
            return 0.0;
        }
        double mach = m_atm.Mach(h, v);
        double cx = CxRetroPlume(CxMachLanding(mach), retroPlumeK, thrFrac);
        return 0.5 * rho * v * v * cx * aRef;
    }

    //  Landing drag force against the air velocity; thrFrac in [0, 1].
    CfBodyVector DragLanding(double h, const CfBodyVector& vAir, double aRef,
                    double retroPlumeK, double thrFrac) const
    {
        double v = Norm(vAir);
        double f = DragLandingMagnitude(h, v, aRef, retroPlumeK, thrFrac);
        CfBodyVector z { };
        if (f <= 0.0)
        {
            return z;
        }
        return Scale(vAir, -f / v);
    }

    CfLocalVector DragLanding(double h, const CfLocalVector& vAir, double aRef,
                    double retroPlumeK, double thrFrac) const
    {
        double v = Norm(vAir);
        double f = DragLandingMagnitude(h, v, aRef, retroPlumeK, thrFrac);
        CfLocalVector z { };
        if (f <= 0.0)
        {
            return z;
        }
        return Scale(vAir, -f / v);
    }

    //  Scalar landing DECELERATION, m/s^2, for point predictions:
    //  a = q Cx(M) A / m (with plume reduction). For planner/onboard.
    double DragDecelLanding(double h, double v, double aRef, double mass,
                            double retroPlumeK, double thrFrac) const
    {
        if (mass < 1.0)
        {
            mass = 1.0;
        }
        return DragLandingMagnitude(h, v, aRef, retroPlumeK, thrFrac) / mass;
    }

    //  Maximum lateral force of a PAIR of fins at full deflection, N:
    //  2 q Cl_delta A_fin delta_max (two effective panels per axis).
    double FinForceMax(double h, double v, double finArea, double finClDelta,
                       double finMaxDeg) const
    {
        double q = m_atm.DynamicPressure(h, v);
        double dmax = finMaxDeg * PI / 180.0;
        return 2.0 * q * finClDelta * finArea * dmax;
    }

    //-- strip model of a cylinder -----------------------------------

    //  Cylinder body in CfBodyVector axes (body axis bu to the nose, be and bn
    //  lateral). Lateral flow over nStrips strips with local velocity
    //  v + w x r (gives lateral force, wind overturning moment and
    //  rotation damping automatically) + axial drag over the midsection
    //  with the Mach correction and retro-plume shielding on descent
    //  (vRel.bu < 0: flow onto the tail).
    //      vRel   - vehicle velocity relative to the air, CfBodyVector, m/s
    //      omg    - angular rate, CfBodyOmega, rad/s
    //      rho    - air density, kg/m^3
    //      mach   - Mach number of the air speed
    //      thrFrac - thrust/maximum in [0, 1] for the plume shielding
    //      length, radius - cylinder geometry, m
    //      cmFromTail - CoM from the tail, m
    //      cdLateral, cdAxial - drag coefficients
    //      retroPlumeK - shielding coefficient (0 - off)
    //      nStrips - number of strips, clamped to [1, 32]
    //  Returns the force (CfBodyVector) and the moment about the CoM (CfBodyMoment).
    void CylinderStrip(const CfBodyVector& vRel, const CfBodyOmega& omg, double rho,
                       double mach, double thrFrac, double length,
                       double radius, double cmFromTail, double cdLateral,
                       double cdAxial, double retroPlumeK, int nStrips,
                       CfBodyVector* force, CfBodyMoment* moment) const
    {
        CfBodyVector f0 { };
        *force = f0;
        CfBodyMoment m0 { };
        *moment = m0;
        if (rho <= 0.0 || length <= 0.0 || radius <= 0.0)
        {
            return;
        }
        if (nStrips < 1)
        {
            nStrips = 1;
        }
        if (nStrips > 32)
        {
            nStrips = 32;
        }

        double dx = length / (double)nStrips;
        double dia = 2.0 * radius;

        //  Lateral flow strip by strip; strips lie along bu, the lateral
        //  part of the local velocity is its (be, bn) components.
        int i = 0;
        for (i = 0; i < nStrips; i = i + 1)
        {
            double xi = ((double)i + 0.5) * dx - cmFromTail;
            CfBodyVector r { };
            r.bu = xi;
            CfBodyVector vloc = Add(vRel, OmegaCrossR(omg, r));
            CfBodyVector vlat { };
            vlat.be = vloc.be;
            vlat.bn = vloc.bn;
            double vm = Norm(vlat);
            if (vm < 1.0e-9)
            {
                continue;
            }
            CfBodyVector df = Scale(vlat, -0.5 * rho * cdLateral * dia * dx * vm);
            *force = Add(*force, df);
            *moment = Add(*moment, Torque(r, df));
        }

        //  Axial drag over the midsection: Mach correction; on descent
        //  (air flow onto the tail, vRel.bu < 0) the retro plume faces
        //  the flow and shields the base.
        double aAx = PI * radius * radius;
        double vx = vRel.bu;
        double cdAx = cdAxial * MachFactorAxial(mach);
        if (vx < 0.0)
        {
            cdAx = CxRetroPlume(cdAx, retroPlumeK, thrFrac);
        }
        force->bu += -0.5 * rho * cdAx * aAx * fabs(vx) * vx;
    }

    //  Bound atmosphere (for direct rho/Mach queries).
    constexpr const EnvAtmosphere& Atmo() const
    {
        return m_atm;
    }

private:

    const EnvAtmosphere& m_atm;
};

//  Constant-initialized aerodynamics bound to the standard atmosphere.
inline constexpr EnvAeroforces AER_STD = EnvAeroforces(ATM_STD);

}   // namespace gnc

#endif  // LIB_AERODYNAMICS_H

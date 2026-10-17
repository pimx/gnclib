//======================================================================
//  lib_aerodynamics.h
//  EnvAeroforces - aerodynamic force models: constant-Cd drag, retro-plume
//  shielding, Prandtl-Glauert correction, ground effect, strip model of a
//  cylinder. Coefficients and their Mach dependence belong to the caller
//  (1.6.0: the tabulated Cx(M) of particular bodies left the library).
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

    //-- corrections ---------------------------------------------------

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

    //-- strip model of a cylinder -----------------------------------

    //  Cylinder body in CfBodyVector axes (body axis bu to the nose, be and bn
    //  lateral). Lateral flow over nStrips strips with local velocity
    //  v + w x r (gives lateral force, wind overturning moment and
    //  rotation damping automatically) + axial drag over the midsection
    //  with retro-plume shielding when the flow comes onto the tail
    //  (vRel.bu < 0). 1.6.0: the axial coefficient is given as it acts --
    //  a Mach correction, if the body needs one, is the caller's (the
    //  argument mach of 1.5.0 is gone).
    //      vRel   - body velocity relative to the air, CfBodyVector, m/s
    //      omg    - angular rate, CfBodyOmega, rad/s
    //      rho    - air density, kg/m^3
    //      thrFrac - thrust/maximum in [0, 1] for the plume shielding
    //      length, radius - cylinder geometry, m
    //      cmFromTail - CoM from the tail, m
    //      cdLateral - lateral drag coefficient
    //      cdAxial - axial drag coefficient as it acts (0 - no axial term)
    //      retroPlumeK - shielding coefficient (0 - off)
    //      nStrips - number of strips, clamped to [1, 32]
    //  Returns the force (CfBodyVector) and the moment about the CoM (CfBodyMoment).
    void CylinderStrip(const CfBodyVector& vRel, const CfBodyOmega& omg, double rho,
                       double thrFrac, double length,
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
        //  part of the local velocity is its (be, bn) components. 1.4.0:
        //  the arm of a strip is r = xi bu, so omega x r = (w.bn xi,
        //  -w.be xi, 0) and r x df = (-xi df.bn, xi df.be, 0) in (be, bn,
        //  bu); the terms of the general products with a zero factor are
        //  left out -- the sums are the same (MEASURED: consumer records
        //  identical to the bit), the loop about a third cheaper.
        double k = -0.5 * rho * cdLateral * dia * dx;
        double fe = 0.0;
        double fn = 0.0;
        double me = 0.0;
        double mn = 0.0;
        int i = 0;
        for (i = 0; i < nStrips; i = i + 1)
        {
            double xi = ((double)i + 0.5) * dx - cmFromTail;
            double ve = vRel.be + omg.bn * xi;
            double vn = vRel.bn - omg.be * xi;
            double vm = sqrt(vn * vn + ve * ve);
            if (vm < 1.0e-9)
            {
                continue;
            }
            double s = k * vm;
            double dfe = ve * s;
            double dfn = vn * s;
            fe = fe + dfe;
            fn = fn + dfn;
            me = me - xi * dfn;
            mn = mn + xi * dfe;
        }
        force->be = fe;
        force->bn = fn;
        moment->be = me;
        moment->bn = mn;

        if (cdAxial == 0.0)
        {
            //  1.4.0: no axial term asked for -- it would add a zero.
            return;
        }
        //  Axial drag over the midsection; with the air flow onto the tail
        //  (vRel.bu < 0) the retro plume faces the flow and shields the
        //  base.
        double aAx = PI * radius * radius;
        double vx = vRel.bu;
        double cdAx = cdAxial;
        if (vx < 0.0)
        {
            cdAx = CxRetroPlume(cdAx, retroPlumeK, thrFrac);
        }
        force->bu += -0.5 * rho * cdAx * aAx * fabs(vx) * vx;
    }

    //  1.4.0: the lateral strip forces twice at once, with the body turning
    //  and with it still (the damping of the rotation is the difference):
    //  the numbers of CylinderStrip(vRel, omg, ...) and CylinderStrip(vRel,
    //  0, ...) with cdAxial 0, in one pass and one square root for the
    //  still body.
    //    vRel, omg, rho, length, radius, cmFromTail, cdLateral, nStrips --
    //                as CylinderStrip
    //    force, moment           -- out: the body turning (bu components 0)
    //    forceStill, momentStill -- out: the body still
    void CylinderStripPair(const CfBodyVector& vRel, const CfBodyOmega& omg, double rho,
                           double length, double radius, double cmFromTail, double cdLateral,
                           int nStrips, CfBodyVector* force, CfBodyMoment* moment,
                           CfBodyVector* forceStill, CfBodyMoment* momentStill) const
    {
        CfBodyVector f0 { };
        CfBodyMoment m0 { };
        *force = f0;
        *moment = m0;
        *forceStill = f0;
        *momentStill = m0;
        if (rho <= 0.0 || length <= 0.0 || radius <= 0.0)
        {
            return;
        }
        nStrips = nStrips < 1 ? 1 : (nStrips > 32 ? 32 : nStrips);
        double dx = length / (double)nStrips;
        double dia = 2.0 * radius;
        double k = -0.5 * rho * cdLateral * dia * dx;
        double se = vRel.be;
        double sn = vRel.bn;
        double sm = sqrt(sn * sn + se * se);
        bool still = sm >= 1.0e-9;
        double ss = k * sm;
        double sfe = se * ss;
        double sfn = sn * ss;
        double fe = 0.0;
        double fn = 0.0;
        double me = 0.0;
        double mn = 0.0;
        double fe0 = 0.0;
        double fn0 = 0.0;
        double me0 = 0.0;
        double mn0 = 0.0;
        int i = 0;
        for (i = 0; i < nStrips; i = i + 1)
        {
            double xi = ((double)i + 0.5) * dx - cmFromTail;
            if (still)
            {
                fe0 = fe0 + sfe;
                fn0 = fn0 + sfn;
                me0 = me0 - xi * sfn;
                mn0 = mn0 + xi * sfe;
            }
            double ve = vRel.be + omg.bn * xi;
            double vn = vRel.bn - omg.be * xi;
            double vm = sqrt(vn * vn + ve * ve);
            if (vm < 1.0e-9)
            {
                continue;
            }
            double s = k * vm;
            double dfe = ve * s;
            double dfn = vn * s;
            fe = fe + dfe;
            fn = fn + dfn;
            me = me - xi * dfn;
            mn = mn + xi * dfe;
        }
        force->be = fe;
        force->bn = fn;
        moment->be = me;
        moment->bn = mn;
        forceStill->be = fe0;
        forceStill->bn = fn0;
        momentStill->be = me0;
        momentStill->bn = mn0;
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

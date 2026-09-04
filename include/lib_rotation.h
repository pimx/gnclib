//======================================================================
//  lib_rotation.h
//  CfRotation - rotation operator between CfLocalVector (n,e,u) and CfBodyVector (bn,be,bu).
//  C++17, no STL, no heap, no exceptions, no virtual.
//  Header-only. CfRotation is a literal type: static instances are
//  constant-initialized, no entry in .init_array.
//======================================================================

#ifndef LIB_ROTATION_H
#define LIB_ROTATION_H

#include <math.h>

#include "lib_frames.h"

namespace gnc
{

//----------------------------------------------------------------------
//  DCM - direction cosine matrix with NAMED elements.
//  Element  X_bY  is the cosine between CfLocalVector axis X and CfBodyVector axis bY,
//  i.e. the contribution of CfBodyVector component bY to CfLocalVector component X.
//
//      CfLocalVector.X  = sum over bY of ( X_bY * CfBodyVector.bY )
//      CfBodyVector.bY = sum over X  of ( X_bY * CfLocalVector.X )     (transpose)
//
//  No row/column index exists, so no index-order mistake is possible.
//----------------------------------------------------------------------
struct DCM
{
    double n_bn;
    double n_be;
    double n_bu;
    double e_bn;
    double e_be;
    double e_bu;
    double u_bn;
    double u_be;
    double u_bu;
};

constexpr DCM MakeDcmIdent()
{
    DCM c { };
    c.n_bn = 1.0;
    c.e_be = 1.0;
    c.u_bu = 1.0;
    return c;
}

//----------------------------------------------------------------------
//  CfRotation - the rotation operator itself.
//
//  Primary storage : unit quaternion (qw; qn, qe, qu).
//  Derived storage : DCM, rebuilt on every state change (no lazy
//                    mutable fields - behaviour is deterministic).
//
//  The vector part of the quaternion is the rotation axis; its
//  components are numerically identical in CfLocalVector and CfBodyVector names,
//  therefore the fields are named after the CfLocalVector axes.
//
//  Sign convention: the operator maps CfBodyVector components into CfLocalVector
//  components (B2L is the "direct" direction).
//----------------------------------------------------------------------
class CfRotation
{
public:

    //  Constant-initialized identity: bu == u, be == e, bn == n.
    constexpr CfRotation()
        : qw(1.0)
        , qn(0.0)
        , qe(0.0)
        , qu(0.0)
        , m(MakeDcmIdent())
    {
    }

    //-- state -------------------------------------------------------

    constexpr void SetIdentity()
    {
        qw = 1.0;
        qn = 0.0;
        qe = 0.0;
        qu = 0.0;
        m = MakeDcmIdent();
    }

    //  Quaternion is normalized on input.
    void SetQuat(double w, double an, double ae, double au)
    {
        qw = w;
        qn = an;
        qe = ae;
        qu = au;
        Normalize();
        BuildDcm();
    }

    //  Rotation by "ang" [rad] about an axis given in CfLocalVector components.
    //  The axis need not be unit length.
    void SetAxisAngle(const CfLocalVector& axis, double ang)
    {
        double len = sqrt(axis.n * axis.n + axis.e * axis.e + axis.u * axis.u);
        if (len < 1.0e-12)
        {
            SetIdentity();
            return;
        }
        double s = sin(0.5 * ang) / len;
        qw = cos(0.5 * ang);
        qn = s * axis.n;
        qe = s * axis.e;
        qu = s * axis.u;
        BuildDcm();
    }

    //  Attitude of a vertically launched vehicle:
    //      azi  - azimuth of the nose, from north toward east [rad]
    //      tlt  - tilt of the nose away from the vertical     [rad]
    //      gam  - roll about the nose axis bu                 [rad]
    //  tlt = 0 gives the initial attitude; there the pair (azi, gam)
    //  is degenerate - the physical gimbal lock of a vertical vehicle,
    //  not a code defect.
    void SetTilt(double azi, double tlt, double gam)
    {
        CfRotation r1;
        CfRotation r2;
        CfRotation r3;
        r1.SetAxisAngle(AXIS_U, -azi);
        r2.SetAxisAngle(AXIS_E, -tlt);
        r3.SetAxisAngle(AXIS_U, gam);
        *this = Mul(Mul(r1, r2), r3);
    }

    //  Restore the operator from an externally supplied DCM
    //  (Shepperd selection of the largest quaternion component).
    void SetDcm(const DCM& c)
    {
        double tr = c.u_bu + c.e_be + c.n_bn;
        double s = 0.0;
        if (tr > 0.0)
        {
            s = 2.0 * sqrt(1.0 + tr);
            qw = 0.25 * s;
            qu = (c.n_be - c.e_bn) / s;
            qe = (c.u_bn - c.n_bu) / s;
            qn = (c.e_bu - c.u_be) / s;
        }
        else if (c.u_bu > c.e_be && c.u_bu > c.n_bn)
        {
            s = 2.0 * sqrt(1.0 + c.u_bu - c.e_be - c.n_bn);
            qw = (c.n_be - c.e_bn) / s;
            qu = 0.25 * s;
            qe = (c.u_be + c.e_bu) / s;
            qn = (c.u_bn + c.n_bu) / s;
        }
        else if (c.e_be > c.n_bn)
        {
            s = 2.0 * sqrt(1.0 + c.e_be - c.u_bu - c.n_bn);
            qw = (c.u_bn - c.n_bu) / s;
            qu = (c.u_be + c.e_bu) / s;
            qe = 0.25 * s;
            qn = (c.e_bn + c.n_be) / s;
        }
        else
        {
            s = 2.0 * sqrt(1.0 + c.n_bn - c.u_bu - c.e_be);
            qw = (c.e_bu - c.u_be) / s;
            qu = (c.u_bn + c.n_bu) / s;
            qe = (c.e_bn + c.n_be) / s;
            qn = 0.25 * s;
        }
        Normalize();
        BuildDcm();
    }

    //-- transforms --------------------------------------------------
    //  One implementation, three type-safe entry points per direction.
    //  The quantity kind is preserved: a moment cannot come out of a
    //  transform fed with an angular rate.

    CfLocalVector B2L(const CfBodyVector& b) const
    {
        return MapB2L<CfLocalVector>(b);
    }

    CfLocalOmega B2L(const CfBodyOmega& b) const
    {
        return MapB2L<CfLocalOmega>(b);
    }

    CfLocalMoment B2L(const CfBodyMoment& b) const
    {
        return MapB2L<CfLocalMoment>(b);
    }

    CfLocalOmegaAcc B2L(const CfBodyOmegaAcc& b) const
    {
        return MapB2L<CfLocalOmegaAcc>(b);
    }

    CfBodyVector L2B(const CfLocalVector& l) const
    {
        return MapL2B<CfBodyVector>(l);
    }

    CfBodyOmega L2B(const CfLocalOmega& l) const
    {
        return MapL2B<CfBodyOmega>(l);
    }

    CfBodyMoment L2B(const CfLocalMoment& l) const
    {
        return MapL2B<CfBodyMoment>(l);
    }

    CfBodyOmegaAcc L2B(const CfLocalOmegaAcc& l) const
    {
        return MapL2B<CfBodyOmegaAcc>(l);
    }

    //-- algebra -----------------------------------------------------

    //  Algebraic inverse: Mul(a, a.Inv()) == identity.
    //  For the transforms themselves it is never needed - B2L and L2B
    //  already give both directions of the same operator.
    CfRotation Inv() const
    {
        CfRotation r;
        r.qw = qw;
        r.qn = -qn;
        r.qe = -qe;
        r.qu = -qu;
        r.BuildDcm();
        return r;
    }

    //  Quaternion product. Semantics:
    //      Mul(a, b).B2L(x) == a.B2L( b.B2L(x) )
    static CfRotation Mul(const CfRotation& a, const CfRotation& b)
    {
        CfRotation r;
        r.qw = a.qw * b.qw - a.qn * b.qn - a.qe * b.qe - a.qu * b.qu;
        r.qu = a.qw * b.qu + b.qw * a.qu + a.qe * b.qn - a.qn * b.qe;
        r.qe = a.qw * b.qe + b.qw * a.qe + a.qn * b.qu - a.qu * b.qn;
        r.qn = a.qw * b.qn + b.qw * a.qn + a.qu * b.qe - a.qe * b.qu;
        r.Normalize();
        r.BuildDcm();
        return r;
    }

    //  Kinematic step by the angular rate of CfBodyVector relative to CfLocalVector,
    //  projected on the body axes. Exact exponential map, so the
    //  quaternion norm does not drift with the number of steps.
    void Propagate(const CfBodyOmega& w, double dt)
    {
        double an = w.bn * dt;
        double ae = w.be * dt;
        double au = w.bu * dt;
        double t2 = an * an + ae * ae + au * au;
        double th = sqrt(t2);
        double cs = 0.0;
        double sc = 0.0;
        if (th < 1.0e-8)
        {
            cs = 1.0 - t2 / 8.0;
            sc = 0.5 - t2 / 48.0;
        }
        else
        {
            cs = cos(0.5 * th);
            sc = sin(0.5 * th) / th;
        }
        CfRotation d;
        d.qw = cs;
        d.qn = sc * an;
        d.qe = sc * ae;
        d.qu = sc * au;
        //  The increment is expressed in the current body axes,
        //  therefore it multiplies from the right.
        *this = Mul(*this, d);
    }

    void Normalize()
    {
        double s = sqrt(qw * qw + qn * qn + qe * qe + qu * qu);
        if (s < 1.0e-12)
        {
            SetIdentity();
            return;
        }
        s = 1.0 / s;
        qw *= s;
        qn *= s;
        qe *= s;
        qu *= s;
    }

    //-- access ------------------------------------------------------

    constexpr const DCM& Mat() const
    {
        return m;
    }

    constexpr double Qw() const
    {
        return qw;
    }

    constexpr double Qn() const
    {
        return qn;
    }

    constexpr double Qe() const
    {
        return qe;
    }

    constexpr double Qu() const
    {
        return qu;
    }

    //  Nose direction in CfLocalVector components (column bu of the DCM).
    constexpr CfLocalVector NoseDir() const
    {
        CfLocalVector l { };
        l.n = m.n_bu;
        l.e = m.e_bu;
        l.u = m.u_bu;
        return l;
    }

    //  Deviation of the quaternion norm from unity (health monitor).
    double NormErr() const
    {
        return sqrt(qw * qw + qn * qn + qe * qe + qu * qu) - 1.0;
    }

private:

    //  TL and TB must expose the component names of their family.
    //  Replace by three explicit functions if the coding standard
    //  forbids templates; the template only guarantees that the
    //  three overloads can never diverge.

    template <class TL, class TB>
    constexpr TL MapB2L(const TB& b) const
    {
        TL l { };
        l.n = m.n_bn * b.bn + m.n_be * b.be + m.n_bu * b.bu;
        l.e = m.e_bn * b.bn + m.e_be * b.be + m.e_bu * b.bu;
        l.u = m.u_bn * b.bn + m.u_be * b.be + m.u_bu * b.bu;
        return l;
    }

    template <class TB, class TL>
    constexpr TB MapL2B(const TL& l) const
    {
        TB b { };
        b.bn = m.n_bn * l.n + m.e_bn * l.e + m.u_bn * l.u;
        b.be = m.n_be * l.n + m.e_be * l.e + m.u_be * l.u;
        b.bu = m.n_bu * l.n + m.e_bu * l.e + m.u_bu * l.u;
        return b;
    }

    void BuildDcm()
    {
        m.u_bu = 1.0 - 2.0 * (qe * qe + qn * qn);
        m.u_be = 2.0 * (qu * qe - qw * qn);
        m.u_bn = 2.0 * (qu * qn + qw * qe);
        m.e_bu = 2.0 * (qu * qe + qw * qn);
        m.e_be = 1.0 - 2.0 * (qu * qu + qn * qn);
        m.e_bn = 2.0 * (qe * qn - qw * qu);
        m.n_bu = 2.0 * (qu * qn - qw * qe);
        m.n_be = 2.0 * (qe * qn + qw * qu);
        m.n_bn = 1.0 - 2.0 * (qu * qu + qe * qe);
    }

    double qw;
    double qn;
    double qe;
    double qu;
    DCM m;
};

//  Constant-initialized identity operator, usable as a static.
inline constexpr CfRotation ROT_IDENT = CfRotation();

}   // namespace gnc

#endif  // LIB_ROTATION_H

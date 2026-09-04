//======================================================================
//  lib_chain.h
//  Composition of the two operators: CfBodyVector <-> CfLocalVector <-> ECEF.
//  Free vectors only. A point in CfBodyVector must first be reduced to a CfLocalPoint
//  by adding the position of the vehicle origin in the local frame,
//  because CfRotation is a pure rotation and carries no translation.
//======================================================================

#ifndef LIB_CHAIN_H
#define LIB_CHAIN_H

#include "lib_rotation.h"
#include "lib_geodesy.h"

namespace gnc
{

inline CfEcefVector B2E(const CfRotation& r, const GEO& g, const CfBodyVector& b)
{
    return g.L2E(r.B2L(b));
}

inline CfBodyVector E2B(const CfRotation& r, const GEO& g, const CfEcefVector& v)
{
    return r.L2B(g.E2L(v));
}

//  Position of a body point: "arm" is the offset of the point from the
//  CfBodyVector origin, "pos" is the position of the CfBodyVector origin in the local frame.
inline CfLocalPoint BodyPointToLop(const CfRotation& r, const CfLocalPoint& pos, const CfBodyVector& arm)
{
    CfLocalVector d = r.B2L(arm);
    CfLocalPoint l { };
    l.n = pos.n + d.n;
    l.e = pos.e + d.e;
    l.u = pos.u + d.u;
    return l;
}

}   // namespace gnc

#endif  // LIB_CHAIN_H

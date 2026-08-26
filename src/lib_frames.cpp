/*============================================================================
 * lib_frames.cpp -- преобразования связанных раскладок (этап 297).
 *============================================================================*/
#include "lib_frames.h"

VecXBody xbodyFromZbody(const VecZBody& zb)
{
    return VecXBody(zb.z, zb.x, zb.y);
}

VecZBody zbodyFromXbody(const VecXBody& xb)
{
    return VecZBody(xb.y, xb.z, xb.x);
}

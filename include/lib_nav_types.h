//============================================================================
// lib_nav_types.h -- measurement and navigation data of the common ESKF
// (gnclib 1.1.0; the layouts of the FLY90 / HOP90 msg_sensors.h, which
// alias these types). C++17, no STL, no heap.
//============================================================================
#ifndef GNCLIB_LIB_NAV_TYPES_H
#define GNCLIB_LIB_NAV_TYPES_H

#include "lib_frames.h"
#include "lib_rotation.h"

namespace gnc
{

typedef double Time; // s (as msg_types.h of the projects)

struct ImuData
{
    CfBodyOmega gyro;   // body rate, rad/s
    CfBodyVector accel; // specific force (body), m/s^2
    ImuData() : gyro { }, accel { }
    {
    }
};

struct GnssData
{
    CfEcefPoint posEcp;  // position, ECEF, m
    CfEcefVector velEcf; // velocity, ECEF, m/s
    int valid;           // 1 -- a fresh measurement at this tick
    GnssData() : posEcp { }, velEcf { }, valid(0)
    {
    }
};

struct AltimeterData
{
    double range; // slant range to the surface, m
    int valid;    // 1 -- valid (within the range of the instrument)
    AltimeterData() : range(0.0), valid(0)
    {
    }
};

struct NavData
{
    CfLocalPoint posLop;  // reference point (site frame), m
    CfLocalVector velLoc; // velocity (site frame), m/s
    CfRotation att;       // attitude BODY <-> LOCAL
    CfBodyOmega omgBdw;   // body rate, rad/s
    double gnssAge;       // time since the last GNSS correction, s
    NavData() : posLop { }, velLoc { }, omgBdw { }, gnssAge(0.0)
    {
    }
};

}   // namespace gnc

#endif // GNCLIB_LIB_NAV_TYPES_H

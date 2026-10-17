//============================================================================
// lib_orbit.h -- орбитальные величины и геометрия площадки на вращающейся
// Земле (1.6.0).
// C++17, без STL, без кучи, без исключений.
//============================================================================
#ifndef GNCLIB_LIB_ORBIT_H
#define GNCLIB_LIB_ORBIT_H

#include <math.h>

#include "lib_frame_ops.h"
#include "lib_frames.h"
#include "lib_geodesy.h"

namespace gnc
{

// Инерциальное движение точки, заданной в осях площадки.
struct InertialState
{
    CfLocalVector vel { }; // инерциальная скорость в осях площадки, м/с
    double speed = 0.0;    // инерциальная скорость, м/с
    double radius = 0.0;   // геоцентрическое расстояние, м
    double fpa = 0.0;      // наклон траектории к плоскости, нормальной радиусу, рад
    double inclination = 0.0; // наклонение плоскости оскулирующей орбиты, рад
};

// Инерциальная скорость, наклон траектории и наклонение по положению и
// скорости относительно вращающейся Земли. Переносная скорость w x r точки
// в ECEF (вращение вокруг z) прибавляется к относительной в осях площадки;
// направление радиуса переводится в оси площадки (наклон -- в одной
// системе); наклонение -- угол между моментом r x v (ECEF) и осью полюса.
//   geo -- система площадки (эллипсоид и вращение); p, v -- положение и
//   скорость относительно Земли в осях площадки
inline InertialState InertialOf(const GEO& geo, const CfLocalPoint& p, const CfLocalVector& v)
{
    InertialState r;
    CfEcefPoint ep = geo.L2E(p);
    double om = geo.Ellipsoid().om;
    CfEcefVector vt { };
    vt.x = -om * ep.y;
    vt.y = om * ep.x;
    r.vel = Add(v, geo.E2L(vt));
    r.speed = Norm(r.vel);
    r.radius = sqrt(ep.x * ep.x + ep.y * ep.y + ep.z * ep.z);
    CfEcefVector er { };
    er.x = ep.x / r.radius;
    er.y = ep.y / r.radius;
    er.z = ep.z / r.radius;
    CfLocalVector up = geo.E2L(er);
    if (r.speed > 1.0e-6)
    {
        double s = (r.vel.n * up.n + r.vel.e * up.e + r.vel.u * up.u) / r.speed;
        r.fpa = asin(s > 1.0 ? 1.0 : (s < -1.0 ? -1.0 : s));
    }
    CfEcefVector w = geo.L2E(r.vel);
    double hx = ep.y * w.z - ep.z * w.y;
    double hy = ep.z * w.x - ep.x * w.z;
    double hz = ep.x * w.y - ep.y * w.x;
    double h = sqrt(hx * hx + hy * hy + hz * hz);
    if (h > 1.0e-9)
    {
        double c = hz / h;
        r.inclination = acos(c > 1.0 ? 1.0 : (c < -1.0 ? -1.0 : c));
    }
    return r;
}

// Круговая скорость на радиусе, м/с.
//   mu -- гравитационный параметр, м^3/с^2; radius -- м
inline double CircularSpeed(double mu, double radius)
{
    return sqrt(mu / radius);
}

// Начальный азимут пуска на заданное наклонение (относительно вращающейся
// Земли). Инерциальный азимут плоскости по сфере: sin A_i = cos(i) /
// cos(широта); он переводится в азимут над вращающейся Землёй вычитанием
// переносной скорости площадки из орбитальной:
//   tan A = (V sin A_i - w R cos(широта)) / (V cos A_i),  V = sqrt(mu / (R + h)).
// Если наклонение меньше широты (корня нет), ближайшее -- точно на восток.
//   lat -- широта площадки, рад; incl -- наклонение, рад
//   orbitAlt -- высота орбиты над экваториальным радиусом, м
//   south -- южная ветвь (A_i = PI - корень); el -- эллипсоид (a, mu, om)
// Возвращает азимут от севера по часовой, рад.
inline double LaunchAzimuthStart(double lat, double incl, double orbitAlt, bool south,
                                 const CfGeoid& el)
{
    double sinA = cos(incl) / cos(lat);
    double a0 = 0.5 * PI;
    if (sinA < 1.0)
    {
        double root = asin(sinA);
        double ai = south ? PI - root : root;
        double v = sqrt(el.mu / (el.a + orbitAlt));
        double vPad = el.om * el.a * cos(lat);
        a0 = atan2(v * sin(ai) - vPad, v * cos(ai));
    }
    return a0;
}

// Точка системы площадки с заданными n, e на заданной высоте над
// эллипсоидом: u следует за кривизной (-0.078 м на 1 км); три поправки u
// вдоль местной вертикали приводят высоту к заданной до 1e-9 м, n и e --
// точно.
//   geo -- система площадки; n, e -- м; alt -- высота над эллипсоидом, м
// Возвращает геодезические координаты точки (alt -- достигнутая высота).
inline CfGeodetic GeoAtLevel(const GEO& geo, double n, double e, double alt)
{
    CfLocalPoint pad { };
    pad.n = n;
    pad.e = e;
    CfGeodetic g = geo.GeoOf(pad);
    for (int i = 0; i < 3; ++i)
    {
        pad.u = pad.u - (g.alt - alt);
        g = geo.GeoOf(pad);
    }
    return g;
}

} // namespace gnc

#endif // GNCLIB_LIB_ORBIT_H

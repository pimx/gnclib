//============================================================================
// lib_launch.h -- стартовые системы координат (1.11.0): ГИСК, ГОСК, НССК,
// ССК и связанные оси по углам тангажа, рыскания и крена.
//
// Системы (справочный документ "Системы координат", пп. 3, 4; ГОСТ 32453):
//   ГОСК -- гринвичская относительная геоцентрическая (= ECEF gnclib);
//   ГИСК -- гринвичская инерциальная: оси ГОСК, зафиксированные в момент
//           старта t0; ГИСК -> ГОСК -- поворот на om (t - t0) вокруг Z (5);
//   ССК  -- стартовая: начало в точке старта, Y -- вверх по нормали к ОЗЭ,
//           X -- по азимуту пуска (горизонтально), Z -- вправо (правая
//           тройка); вращается с Землёй;
//   НССК -- начальная стартовая: оси ССК, зафиксированные в момент старта
//           (инерциальная); ГИСК -> НССК -- матрица (7).
// Система площадки gnc::GEO (север, восток, верх) вращается с Землёй; её
// начало -- точка старта, если GEO::SetOrigin задан точкой старта. Поэтому
// ССК -- поворот осей GEO на азимут, а НССК в осях GEO на момент t --
// поворот осей ССК на -om (t - t0) вокруг оси Земли.
//
// Связанные оси относительно НССК -- последовательность "тангаж (вокруг Z),
// рыскание (вокруг Y'), крен (вокруг X)", матрица НССК -> СвСК (8).
//
// Подключается явно (не через lib.h). C++17, без STL, без кучи.
//============================================================================
#ifndef GNCLIB_LIB_LAUNCH_H
#define GNCLIB_LIB_LAUNCH_H

#include <math.h>

#include "lib_frames.h"
#include "lib_geodesy.h"

namespace gnc
{

// Составляющие вектора в осях стартовой системы (ССК или НССК): X -- по
// азимуту, Y -- вверх, Z -- вправо
struct CfLaunchVector
{
    double x;
    double y;
    double z;
};

// Оси стартовой системы в ECEF на момент старта (= оси НССК в ГИСК) и
// угловая скорость Земли
struct LaunchAxes
{
    CfEcefVector x; // по азимуту пуска
    CfEcefVector y; // вверх по нормали к ОЗЭ
    CfEcefVector z; // вправо
    double om;      // угловая скорость Земли, рад/с (эллипсоид GEO)
};

// Оси стартовой системы по системе площадки и азимуту пуска: строки
// матрицы (7) ГИСК -> НССК.
//   geo -- система площадки (начало -- точка старта); az -- азимут пуска от
//   севера по часовой, рад
inline LaunchAxes LaunchAxesOf(const GEO& geo, double az)
{
    CfLocalVector x { };
    x.n = cos(az);
    x.e = sin(az);
    CfLocalVector y { };
    y.u = 1.0;
    CfLocalVector z { };
    z.n = -sin(az);
    z.e = cos(az);
    LaunchAxes a { };
    a.x = geo.L2E(x);
    a.y = geo.L2E(y);
    a.z = geo.L2E(z);
    a.om = geo.Ellipsoid().om;
    return a;
}

// ГИСК -> ГОСК (ECEF) на момент t (5): поворот на om t вокруг оси Z.
//   v -- вектор в ГИСК; om -- угловая скорость Земли, рад/с; t -- t - t0, с
inline CfEcefVector InertialToEcef(const CfEcefVector& v, double om, double t)
{
    double c = cos(om * t);
    double s = sin(om * t);
    CfEcefVector r { };
    r.x = c * v.x + s * v.y;
    r.y = -s * v.x + c * v.y;
    r.z = v.z;
    return r;
}

// ГОСК (ECEF) -> ГИСК на момент t (обратный поворот (5)).
//   v -- вектор в ECEF; om -- рад/с; t -- t - t0, с
inline CfEcefVector EcefToInertial(const CfEcefVector& v, double om, double t)
{
    double c = cos(om * t);
    double s = sin(om * t);
    CfEcefVector r { };
    r.x = c * v.x - s * v.y;
    r.y = s * v.x + c * v.y;
    r.z = v.z;
    return r;
}

// Скорость относительно Земли в проекциях на ГИСК по абсолютной (6):
// v - om x r, om вдоль оси Z.
//   r, v -- радиус-вектор и абсолютная скорость в ГИСК; om -- рад/с
inline CfEcefVector InertialVelRel(const CfEcefVector& r, const CfEcefVector& v, double om)
{
    CfEcefVector w { };
    w.x = v.x + om * r.y;
    w.y = v.y - om * r.x;
    w.z = v.z;
    return w;
}

// Вектор, заданный в НССК (неподвижный в инерциальном пространстве), в
// местных осях площадки на момент t: НССК -> ГИСК -> ГОСК (5) -> GEO::E2L.
//   geo -- система площадки; a -- оси стартовой системы; v -- вектор в
//   НССК; t -- время от старта, с
inline CfLocalVector NsskToLoc(const GEO& geo, const LaunchAxes& a, const CfLaunchVector& v,
                               double t)
{
    CfEcefVector e { };
    e.x = v.x * a.x.x + v.y * a.y.x + v.z * a.z.x;
    e.y = v.x * a.x.y + v.y * a.y.y + v.z * a.z.y;
    e.z = v.x * a.x.z + v.y * a.y.z + v.z * a.z.z;
    return geo.E2L(InertialToEcef(e, a.om, t));
}

// Вектор, заданный в ССК (неподвижный относительно Земли), в местных осях.
//   geo -- система площадки; a -- оси стартовой системы; v -- вектор в ССК
inline CfLocalVector SskToLoc(const GEO& geo, const LaunchAxes& a, const CfLaunchVector& v)
{
    return NsskToLoc(geo, a, v, 0.0);
}

// Вектор (или точка относительно точки старта) из местных осей в ССК.
//   geo -- система площадки; a -- оси стартовой системы; v -- местный вектор
inline CfLaunchVector LocToSsk(const GEO& geo, const LaunchAxes& a, const CfLocalVector& v)
{
    CfEcefVector e = geo.L2E(v);
    CfLaunchVector r { };
    r.x = e.x * a.x.x + e.y * a.x.y + e.z * a.x.z;
    r.y = e.x * a.y.x + e.y * a.y.y + e.z * a.y.z;
    r.z = e.x * a.z.x + e.y * a.z.y + e.z * a.z.z;
    return r;
}

// Связанные оси X, Y, Z в осях НССК по углам тангажа, рыскания и крена:
// строки матрицы (8) НССК -> СвСК (при psi = gam = 0: X = (cos th, sin th, 0)).
//   th, psi, gam -- углы, рад; x, y, z -- результат
inline void LaunchBodyAxes(double th, double psi, double gam, CfLaunchVector* x,
                           CfLaunchVector* y, CfLaunchVector* z)
{
    double st = sin(th);
    double ct = cos(th);
    double sp = sin(psi);
    double cp = cos(psi);
    double sg = sin(gam);
    double cg = cos(gam);
    x->x = ct * cp;
    x->y = st * cp;
    x->z = -sp;
    y->x = -st * cg + ct * sp * sg;
    y->y = ct * cg + st * sp * sg;
    y->z = cp * sg;
    z->x = st * sg + ct * sp * cg;
    z->y = -ct * sg + st * sp * cg;
    z->z = cp * cg;
}

} // namespace gnc

#endif // GNCLIB_LIB_LAUNCH_H

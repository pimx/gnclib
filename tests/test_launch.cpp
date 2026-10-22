// test_launch.cpp -- gnclib 1.11.0: стартовые системы координат (lib_launch.h)
// против формул справочного документа "Системы координат" (5)-(8), (12).
// Build: g++ -std=c++17 -O2 -Iinclude tests/test_launch.cpp -o build/test_launch
#include "lib.h"
#include "lib_earth.h"
#include "lib_launch.h"
#include <stdio.h>
#include <math.h>

using namespace gnc;

static int g_bad = 0;

// Строка проверки: значение по модулю не больше предела.
static void Check(const char* name, double val, double lim)
{
    int ok = (fabs(val) <= lim) ? 1 : 0;
    printf("%-58s %12.3e  (limit %.1e)  %s\n", name, val, lim, ok ? "OK" : "FAIL");
    g_bad += ok ? 0 : 1;
}

// Модуль разности векторов ECEF.
static double Dif(const CfEcefVector& a, double x, double y, double z)
{
    return sqrt((a.x - x) * (a.x - x) + (a.y - y) * (a.y - y) + (a.z - z) * (a.z - z));
}

// Модуль разности местных векторов.
static double DifL(const CfLocalVector& a, const CfLocalVector& b)
{
    return sqrt((a.e - b.e) * (a.e - b.e) + (a.n - b.n) * (a.n - b.n) + (a.u - b.u) * (a.u - b.u));
}

// Местный вектор v, повёрнутый вокруг единичной оси k на угол ang (Родриг).
static CfLocalVector Turn(const CfLocalVector& v, const CfLocalVector& k, double ang)
{
    double c = cos(ang);
    double s = sin(ang);
    double kv = k.e * v.e + k.n * v.n + k.u * v.u;
    CfLocalVector r { };
    r.e = v.e * c + (k.n * v.u - k.u * v.n) * s + k.e * kv * (1.0 - c);
    r.n = v.n * c + (k.u * v.e - k.e * v.u) * s + k.n * kv * (1.0 - c);
    r.u = v.u * c + (k.e * v.n - k.n * v.e) * s + k.u * kv * (1.0 - c);
    return r;
}

// Проверки (a)..(f) на площадке Плесецк, азимут 14.499029 град.
int main()
{
    const double D2R = PI / 180.0;
    double b = 62.928008 * D2R;
    double l = 40.456137 * D2R;
    double az = 14.499029 * D2R;
    GEO geo;
    CfGeodetic site { };
    site.lat = b;
    site.lon = l;
    geo.SetOrigin(site);
    LaunchAxes a = LaunchAxesOf(geo, az);
    // (a) строки матрицы (7) справочного документа
    {
        double sa = sin(az);
        double ca = cos(az);
        double sb = sin(b);
        double cb = cos(b);
        double sl = sin(l);
        double cl = cos(l);
        double ex = Dif(a.x, -sa * sl - ca * sb * cl, sa * cl - ca * sb * sl, ca * cb);
        double ey = Dif(a.y, cb * cl, sl * cb, sb);
        double ez = Dif(a.z, -ca * sl + sa * sb * cl, ca * cl + sa * sb * sl, -sa * cb);
        Check("(a) axes vs matrix (7): X", ex, 1.0e-15);
        Check("(a) axes vs matrix (7): Y", ey, 1.0e-15);
        Check("(a) axes vs matrix (7): Z", ez, 1.0e-15);
        // правая тройка: X x Y = Z
        double cx = a.x.y * a.y.z - a.x.z * a.y.y;
        double cy = a.x.z * a.y.x - a.x.x * a.y.z;
        double cz = a.x.x * a.y.y - a.x.y * a.y.x;
        Check("(a) right-handed: |X x Y - Z|", Dif(a.z, cx, cy, cz), 1.0e-15);
    }
    // (b) ГИСК <-> ГОСК (5) и скорость относительно Земли (6)
    {
        CfEcefVector v { };
        v.x = 1000.0;
        v.y = -2000.0;
        v.z = 3000.0;
        CfEcefVector w = EcefToInertial(InertialToEcef(v, a.om, 158.23), a.om, 158.23);
        Check("(b) EcefToInertial(InertialToEcef(v)) - v, m", Dif(w, v.x, v.y, v.z), 1.0e-11);
        CfEcefPoint p = geo.OriginEcp();
        CfEcefVector r { };
        r.x = p.x;
        r.y = p.y;
        r.z = p.z;
        CfEcefVector vp { };
        vp.x = -a.om * r.y;
        vp.y = a.om * r.x;
        CfEcefVector vr = InertialVelRel(r, vp, a.om);
        Check("(b) point at rest on the Earth: velocity rel. (6), m/s", Dif(vr, 0.0, 0.0, 0.0), 1.0e-9);
    }
    // (c) ССК <-> местные оси
    {
        CfLaunchVector v { };
        v.x = 300.0;
        v.y = 2700.0;
        v.z = -159.7;
        CfLaunchVector w = LocToSsk(geo, a, SskToLoc(geo, a, v));
        double e = sqrt((w.x - v.x) * (w.x - v.x) + (w.y - v.y) * (w.y - v.y)
                        + (w.z - v.z) * (w.z - v.z));
        Check("(c) LocToSsk(SskToLoc(v)) - v, m", e, 1.0e-10);
        CfLocalVector loc = SskToLoc(geo, a, v);
        CfLocalVector want { };
        want.n = v.x * cos(az) - v.z * sin(az);
        want.e = v.x * sin(az) + v.z * cos(az);
        want.u = v.y;
        Check("(c) SskToLoc: X by the azimuth, Y up, Z right, m", DifL(loc, want), 1.0e-10);
    }
    // (d) направление, неподвижное в НССК, -- в местных осях поворот на
    // -om t вокруг оси Земли
    {
        CfLocalVector k = EarthRateLoc(geo);
        double om = sqrt(k.e * k.e + k.n * k.n + k.u * k.u);
        k.e = k.e / om;
        k.n = k.n / om;
        k.u = k.u / om;
        CfLaunchVector v { };
        v.x = cos(85.5 * D2R);
        v.y = sin(85.5 * D2R);
        double t = 158.23;
        CfLocalVector got = NsskToLoc(geo, a, v, t);
        CfLocalVector want = Turn(SskToLoc(geo, a, v), k, -a.om * t);
        Check("(d) NSSK direction at 158 s vs turn about the Earth axis", DifL(got, want), 1.0e-14);
        Check("(d) NSSK at t = 0 equals SSK", DifL(NsskToLoc(geo, a, v, 0.0),
                                                   SskToLoc(geo, a, v)), 1.0e-16);
    }
    // (e) связанные оси (8): ортонормированность и частный случай psi = gam = 0
    {
        CfLaunchVector x { };
        CfLaunchVector y { };
        CfLaunchVector z { };
        LaunchBodyAxes(0.3, -0.2, 0.7, &x, &y, &z);
        double xy = x.x * y.x + x.y * y.y + x.z * y.z;
        double xz = x.x * z.x + x.y * z.y + x.z * z.z;
        double cz = x.y * y.z - x.z * y.y - z.x;
        Check("(e) body axes orthogonal: X.Y", xy, 1.0e-15);
        Check("(e) body axes orthogonal: X.Z", xz, 1.0e-15);
        Check("(e) body axes right-handed: (X x Y - Z).x", cz, 1.0e-15);
        LaunchBodyAxes(85.5 * D2R, 0.0, 0.0, &x, &y, &z);
        Check("(e) psi = gam = 0: X = (cos th, sin th, 0)",
              fabs(x.x - cos(85.5 * D2R)) + fabs(x.y - sin(85.5 * D2R)) + fabs(x.z), 1.0e-16);
    }
    // (f) точка ССК по ГОСК: (12) через GeoToEcp и матрица (7)
    {
        CfLocalPoint p { };
        p.n = 250.0;
        p.e = -40.0;
        p.u = 1500.0;
        CfEcefPoint ep = geo.L2E(p);
        CfEcefPoint e0 = geo.OriginEcp();
        double dx = ep.x - e0.x;
        double dy = ep.y - e0.y;
        double dz = ep.z - e0.z;
        CfLaunchVector s = LocToSsk(geo, a, vecOf(p));
        double x = dx * a.x.x + dy * a.x.y + dz * a.x.z;
        double y = dx * a.y.x + dy * a.y.y + dz * a.y.z;
        double z = dx * a.z.x + dy * a.z.y + dz * a.z.z;
        Check("(f) SSK point: LocToSsk vs matrix (7) on ECEF difference, m",
              fabs(s.x - x) + fabs(s.y - y) + fabs(s.z - z), 1.0e-8);
    }
    printf(g_bad == 0 ? "TEST_LAUNCH PASS\n" : "TEST_LAUNCH FAIL\n");
    return g_bad;
}

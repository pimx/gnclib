#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

/* lib_earth.cpp -- геодезия WGS-84, гравитация, ENU-фрейм площадки. */
#include "lib_earth.h"

#include <math.h>

/* ------------------------- EarthGeodesy ------------------------- */

/* Радиус кривизны первого вертикала N(lat). */
double EarthGeodesy::primeVerticalRadius(double lat)
{
    double sl = sin(lat);
    return Wgs84::a() / sqrt(1.0 - Wgs84::e2() * sl * sl);
}

/* Радиус кривизны меридиана M(lat). */
double EarthGeodesy::meridianRadius(double lat)
{
    double sl = sin(lat);
    double w2 = 1.0 - Wgs84::e2() * sl * sl;
    return Wgs84::a() * (1.0 - Wgs84::e2()) / (w2 * sqrt(w2));
}

/* Геоцентрическая широта по геодезической (на поверхности эллипсоида):
 * tan(phi_c) = (1 - e^2) tan(phi). */
double EarthGeodesy::geocentricLat(double geodeticLat)
{
    return atan((1.0 - Wgs84::e2()) * tan(geodeticLat));
}

/* Геодезические -> ECEF. */
VecEcef EarthGeodesy::geodeticToEcef(Geodetic g)
{
    double sl = sin(g.lat), cl = cos(g.lat), so = sin(g.lon), co = cos(g.lon);
    double N = primeVerticalRadius(g.lat);
    return VecEcef((N + g.alt) * cl * co, (N + g.alt) * cl * so,
                   (N * (1.0 - Wgs84::e2()) + g.alt) * sl);
}

/* ECEF -> геодезические: итерации Боуринга по широте. */
Geodetic EarthGeodesy::ecefToGeodetic(VecEcef r)
{
    double p = sqrt(r.x * r.x + r.y * r.y);
    double lon = atan2(r.y, r.x);
    double lat = atan2(r.z, p * (1.0 - Wgs84::e2()));
    for (int i = 0; i < 6; i++)
    {
        double N = primeVerticalRadius(lat);
        double alt = p / cos(lat) - N;
        lat = atan2(r.z, p * (1.0 - Wgs84::e2() * N / (N + alt)));
    }
    double N = primeVerticalRadius(lat);
    Geodetic g;
    g.lat = lat;
    g.lon = lon;
    g.alt = p / cos(lat) - N;
    return g;
}

/* Строки DCM -- орты ENU, выраженные в ECEF. */
Mat3 EarthGeodesy::ecefToEnuDcm(double lat, double lon)
{
    double sl = sin(lat), cl = cos(lat), so = sin(lon), co = cos(lon);
    Mat3 R;
    R.m[0][0] = -so;
    R.m[0][1] = co;
    R.m[0][2] = 0.0;
    R.m[1][0] = -sl * co;
    R.m[1][1] = -sl * so;
    R.m[1][2] = cl;
    R.m[2][0] = cl * co;
    R.m[2][1] = cl * so;
    R.m[2][2] = sl;
    return R;
}

/* ------------------------- EarthGravity ------------------------- */

/* Центральное поле + J2, в ECEF (без центробежного члена). */
VecEcef EarthGravity::gravitationEcefJ2(VecEcef r)
{
    double rr = r.norm();
    if (rr < 1.0)
    {
        return VecEcef::zero();
    }
    double f = -Wgs84::mu() / (rr * rr * rr);
    double zr2 = (r.z * r.z) / (rr * rr);
    double j2c = 1.5 * Wgs84::j2() * (Wgs84::a() * Wgs84::a()) / (rr * rr);
    return VecEcef(f * r.x * (1.0 - j2c * (5.0 * zr2 - 1.0)),
                   f * r.y * (1.0 - j2c * (5.0 * zr2 - 1.0)),
                   f * r.z * (1.0 - j2c * (5.0 * zr2 - 3.0)));
}

/* Центробежное ускорение точки, неподвижной в ECEF:
 * a = -w x (w x r) = w^2 * (x, y, 0). */
VecEcef EarthGravity::centrifugalEcef(VecEcef r)
{
    double w2 = Wgs84::omega() * Wgs84::omega();
    return VecEcef(w2 * r.x, w2 * r.y, 0.0);
}

/* Полная сила тяжести: гравитация + центробежная. */
VecEcef EarthGravity::gravityEcef(VecEcef r)
{
    return gravitationEcefJ2(r).add(centrifugalEcef(r));
}

/* Нормальная гравитация Сомильяны (WGS-84) + высотная поправка. */
double EarthGravity::normal(double lat, double h)
{
    double s2 = sin(lat);
    s2 *= s2;
    double g0 = 9.7803253359 * (1.0 + 0.00193185265241 * s2) /
                sqrt(1.0 - 0.00669437999013 * s2);
    /* поправка свободного воздуха 2-го порядка */
    return g0 - (3.0877e-6 - 4.4e-9 * s2) * h + 7.2e-13 * h * h;
}

/* --------------------------- GeoFrame --------------------------- */

GeoFrame::GeoFrame()
{
    m_originEcef = VecEcef::zero();
    m_Ref = Mat3::eye();
    m_Rfe = Mat3::eye();
    m_gEnu = VecEnu(0.0, 0.0, -9.80665);
    m_lat = 0.0;
    m_lon = 0.0;
}

GeoFrame::GeoFrame(Geodetic site)
{
    reset(site);
}

/* Привязка ENU-фрейма к площадке: начало, повороты и константная сила
 * тяжести (J2 + центробежная -- статическая площадка), повёрнутая в ENU. */
void GeoFrame::reset(Geodetic site)
{
    m_lat = site.lat;
    m_lon = site.lon;
    m_originEcef = EarthGeodesy::geodeticToEcef(site);
    m_Ref = EarthGeodesy::ecefToEnuDcm(site.lat, site.lon);
    m_Rfe = m_Ref.transpose();
    m_gEnu =
        VecEnu::of(m_Ref.mulv(EarthGravity::gravityEcef(m_originEcef).vec()));
}

// Матрицы m_Ref (ECEF->ENU) и m_Rfe (ENU->ECEF) -- обычные повороты
// между РАЗНЫМИ геоцентрическими/локальными системами, поэтому переход
// выполняется покомпонентно через Vec3-мост и типизируется на выходе.
VecEnu GeoFrame::ecefPosToEnu(VecEcef rEcef) const
{
    return VecEnu::of(m_Ref.mulv(rEcef.sub(m_originEcef).vec()));
}
VecEcef GeoFrame::enuPosToEcef(VecEnu rEnu) const
{
    return m_originEcef.add(VecEcef::of(m_Rfe.mulv(rEnu.vec())));
}
VecEnu GeoFrame::ecefVecToEnu(VecEcef vEcef) const
{
    return VecEnu::of(m_Ref.mulv(vEcef.vec()));
}
VecEcef GeoFrame::enuVecToEcef(VecEnu vEnu) const
{
    return VecEcef::of(m_Rfe.mulv(vEnu.vec()));
}

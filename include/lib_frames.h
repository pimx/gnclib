/*============================================================================
 * lib_frames.h -- системы координат библиотеки (этап 297: перенесены из
 * lib_linalg.h). Типизированные фреймы: VecEnu (восток-север-верх),
 * VecXBody (связанная, ось X вдоль оси ракеты), VecEcef, VecGosk, а
 * также VecXBody и VecZBody (связанная, ось Z вдоль
 * оси ракеты). На старте при вертикальной ориентации ракеты оси
 * VecZBody совпадают с E, N, U (x=E, y=N, z=U), оси VecXBody -- с
 * U, E, N (x=U, y=E, z=N). Здесь же -- тела методов сопряжения
 * Mat3/Quat с фреймами.
 *============================================================================*/
#ifndef LIB_FRAMES_H
#define LIB_FRAMES_H

#include "lib_linalg.h"

struct VecUen; /* объявления: мосты определены после структур */
struct VecNed;

struct VecEnu
{
    double e, n, u;

    VecEnu() : e(0.0), n(0.0), u(0.0)
    {
    }
    VecEnu(double e_, double n_, double u_) : e(e_), n(n_), u(u_)
    {
    }

    VecEnu add(VecEnu b) const
    {
        return VecEnu(e + b.e, n + b.n, u + b.u);
    }
    VecEnu sub(VecEnu b) const
    {
        return VecEnu(e - b.e, n - b.n, u - b.u);
    }
    VecEnu scale(double s) const
    {
        return VecEnu(e * s, n * s, u * s);
    }
    double dot(VecEnu b) const
    {
        return e * b.e + n * b.n + u * b.u;
    }
    VecEnu cross(VecEnu b) const
    {
        return VecEnu(n * b.u - u * b.n, u * b.e - e * b.u, e * b.n - n * b.e);
    }
    double norm() const
    {
        return sqrt(dot(*this));
    }
    double norm2() const
    {
        return dot(*this);
    }
    VecEnu neg() const
    {
        return VecEnu(-e, -n, -u);
    }
    VecEnu unit() const
    {
        double q = norm();
        if (q < 1e-300)
        {
            return VecEnu(0.0, 0.0, 0.0);
        }
        return scale(1.0 / q);
    }
    VecEnu clampNorm(double maxNorm) const
    {
        double q = norm();
        if (q <= maxNorm || q < 1e-300)
        {
            return *this;
        }
        return scale(maxNorm / q);
    }
    VecEnu lerp(VecEnu b, double t) const
    {
        return VecEnu(e + (b.e - e) * t, n + (b.n - n) * t, u + (b.u - u) * t);
    }
    static VecEnu zero()
    {
        return VecEnu(0.0, 0.0, 0.0);
    }

    /* Переупорядочения компонент той же площадки (определены ниже):
     * VecUen (U,E,N) и VecNed (N,E,D = -U). Точные и обратимые. */
    VecUen uen() const;
    VecNed ned() const;

    /* Мосты к безымянной тройке (порядок E, N, U). */
    Vec3 vec() const
    {
        return Vec3(e, n, u);
    }
    static VecEnu of(Vec3 v)
    {
        return VecEnu(v.x, v.y, v.z);
    }
};

/*============================================================================
 * ПЕРЕУПОРЯДОЧЕНИЯ ЛОКАЛЬНОЙ ПЛОЩАДКИ: VecUen и VecNed.
 *
 * Обе структуры описывают ТУ ЖЕ физическую площадку, что и VecEnu, --
 * меняется лишь порядок (и для NED знак вертикали) хранения компонент:
 *   VecUen: u (вверх), e (восток), n (север);  u x e = n -- правая;
 *   VecNed: n (север), e (восток), d (вниз, d = -u);  n x e = d -- правая.
 * Переходы -- точные перестановки без вращения, полностью обратимые.
 *
 * Замечание о представлении ориентации: матрица/кватернион -- это
 * запись поворота В ПАРЕ БАЗИСОВ. Стартовая ориентация аппарата с
 * продольной осью x вверх в паре (XBody, VecEnu) -- перестановочный
 * кватернион (0.5,-0.5,-0.5,-0.5), а в паре (XBody, VecUen) была бы
 * единичной: физика инвариантна, меняется представление. Сопряжение
 * ориентации в библиотеке зафиксировано парой (VecXBody, VecEnu);
 * VecUen/VecNed -- конвертационные типы для вывода, обмена и стыковки
 * с внешними системами (NED -- авиационная и JSBSim-локальная).
 *============================================================================*/
struct VecUen
{
    double u, e, n;

    VecUen() : u(0.0), e(0.0), n(0.0)
    {
    }
    VecUen(double u_, double e_, double n_) : u(u_), e(e_), n(n_)
    {
    }

    VecUen add(VecUen b) const
    {
        return VecUen(u + b.u, e + b.e, n + b.n);
    }
    VecUen sub(VecUen b) const
    {
        return VecUen(u - b.u, e - b.e, n - b.n);
    }
    VecUen scale(double s) const
    {
        return VecUen(u * s, e * s, n * s);
    }
    double dot(VecUen b) const
    {
        return u * b.u + e * b.e + n * b.n;
    }
    double norm() const
    {
        return sqrt(dot(*this));
    }
    VecUen neg() const
    {
        return VecUen(-u, -e, -n);
    }
    static VecUen zero()
    {
        return VecUen(0.0, 0.0, 0.0);
    }

    /* Обратно в порядок ENU (точно). */
    VecEnu enu() const
    {
        return VecEnu(e, n, u);
    }
};

struct VecNed
{
    double n, e, d;

    VecNed() : n(0.0), e(0.0), d(0.0)
    {
    }
    VecNed(double n_, double e_, double d_) : n(n_), e(e_), d(d_)
    {
    }

    VecNed add(VecNed b) const
    {
        return VecNed(n + b.n, e + b.e, d + b.d);
    }
    VecNed sub(VecNed b) const
    {
        return VecNed(n - b.n, e - b.e, d - b.d);
    }
    VecNed scale(double s) const
    {
        return VecNed(n * s, e * s, d * s);
    }
    double dot(VecNed b) const
    {
        return n * b.n + e * b.e + d * b.d;
    }
    double norm() const
    {
        return sqrt(dot(*this));
    }
    VecNed neg() const
    {
        return VecNed(-n, -e, -d);
    }
    static VecNed zero()
    {
        return VecNed(0.0, 0.0, 0.0);
    }

    /* Обратно в ENU: e = e, n = n, u = -d (точно). */
    VecEnu enu() const
    {
        return VecEnu(e, n, -d);
    }
};

inline VecUen VecEnu::uen() const
{
    return VecUen(u, e, n);
}
inline VecNed VecEnu::ned() const
{
    return VecNed(n, e, -u);
}

struct VecXBody
{
    double x, y, z;

    VecXBody() : x(0.0), y(0.0), z(0.0)
    {
    }
    VecXBody(double x_, double y_, double z_) : x(x_), y(y_), z(z_)
    {
    }

    VecXBody add(VecXBody b) const
    {
        return VecXBody(x + b.x, y + b.y, z + b.z);
    }
    VecXBody sub(VecXBody b) const
    {
        return VecXBody(x - b.x, y - b.y, z - b.z);
    }
    VecXBody scale(double s) const
    {
        return VecXBody(x * s, y * s, z * s);
    }
    double dot(VecXBody b) const
    {
        return x * b.x + y * b.y + z * b.z;
    }
    VecXBody cross(VecXBody b) const
    {
        return VecXBody(y * b.z - z * b.y, z * b.x - x * b.z,
                        x * b.y - y * b.x);
    }
    double norm() const
    {
        return sqrt(dot(*this));
    }
    double norm2() const
    {
        return dot(*this);
    }
    VecXBody neg() const
    {
        return VecXBody(-x, -y, -z);
    }
    VecXBody unit() const
    {
        double q = norm();
        if (q < 1e-300)
        {
            return VecXBody(0.0, 0.0, 0.0);
        }
        return scale(1.0 / q);
    }
    VecXBody clampNorm(double maxNorm) const
    {
        double q = norm();
        if (q <= maxNorm || q < 1e-300)
        {
            return *this;
        }
        return scale(maxNorm / q);
    }
    VecXBody lerp(VecXBody b, double t) const
    {
        return VecXBody(x + (b.x - x) * t, y + (b.y - y) * t,
                        z + (b.z - z) * t);
    }
    static VecXBody zero()
    {
        return VecXBody(0.0, 0.0, 0.0);
    }

    /* Мосты к безымянной тройке (порядок X, Y, Z). */
    Vec3 vec() const
    {
        return Vec3(x, y, z);
    }
    static VecXBody of(Vec3 v)
    {
        return VecXBody(v.x, v.y, v.z);
    }
};

/*============================================================================
 * ГЕОЦЕНТРИЧЕСКИЕ СИСТЕМЫ: ECEF и ГОСК (ГОСТ 20058-80).
 *
 * VecEcef -- гринвичская прямоугольная (WGS-84): начало в центре масс
 * Земли, ось x -- на пересечение Гринвичского меридиана и экватора,
 * ось z -- по оси вращения на север, ось y дополняет до правой
 * (лямбда = 90 град в.д.).
 *
 * VecGosk -- ГРИНВИЧСКАЯ СИСТЕМА ПО ГОСТ 20058-80 (ГОСК): та же точка
 * пространства в аэрокосмической конвенции ГОСТ, где вертикальной осью
 * служит y: ось x -- на пересечение Гринвичского меридиана и экватора
 * (как в ECEF), ось y -- по оси вращения на север, ось z дополняет до
 * ПРАВОЙ системы (лямбда = 90 град з.д.).
 *
 * Системы ЭКВИВАЛЕНТНЫ: связаны перестановкой осей без вращения и без
 * изменения начала, поэтому переход точен и обратим, а норма и взаимные
 * углы сохраняются:
 *     x_g = x_e,  y_g = z_e,  z_g = -y_e
 *     x_e = x_g,  y_e = -z_g, z_e = y_g
 * Обе тройки правые: x_g * y_g = z_g (проверяется в lib_check).
 *============================================================================*/
struct VecGosk;

struct VecEcef
{
    double x, y, z;

    VecEcef() : x(0.0), y(0.0), z(0.0)
    {
    }
    VecEcef(double x_, double y_, double z_) : x(x_), y(y_), z(z_)
    {
    }

    VecEcef add(VecEcef b) const
    {
        return VecEcef(x + b.x, y + b.y, z + b.z);
    }
    VecEcef sub(VecEcef b) const
    {
        return VecEcef(x - b.x, y - b.y, z - b.z);
    }
    VecEcef scale(double s) const
    {
        return VecEcef(x * s, y * s, z * s);
    }
    double dot(VecEcef b) const
    {
        return x * b.x + y * b.y + z * b.z;
    }
    VecEcef cross(VecEcef b) const
    {
        return VecEcef(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }
    double norm() const
    {
        return sqrt(dot(*this));
    }
    double norm2() const
    {
        return dot(*this);
    }
    VecEcef neg() const
    {
        return VecEcef(-x, -y, -z);
    }
    VecEcef unit() const
    {
        double q = norm();
        if (q < 1e-300)
        {
            return VecEcef(0.0, 0.0, 0.0);
        }
        return scale(1.0 / q);
    }
    static VecEcef zero()
    {
        return VecEcef(0.0, 0.0, 0.0);
    }

    /* Переход в ГОСК (определён после VecGosk). */
    VecGosk gosk() const;

    Vec3 vec() const
    {
        return Vec3(x, y, z);
    }
    static VecEcef of(Vec3 v)
    {
        return VecEcef(v.x, v.y, v.z);
    }
};

struct VecGosk
{
    double x, y, z;

    VecGosk() : x(0.0), y(0.0), z(0.0)
    {
    }
    VecGosk(double x_, double y_, double z_) : x(x_), y(y_), z(z_)
    {
    }

    VecGosk add(VecGosk b) const
    {
        return VecGosk(x + b.x, y + b.y, z + b.z);
    }
    VecGosk sub(VecGosk b) const
    {
        return VecGosk(x - b.x, y - b.y, z - b.z);
    }
    VecGosk scale(double s) const
    {
        return VecGosk(x * s, y * s, z * s);
    }
    double dot(VecGosk b) const
    {
        return x * b.x + y * b.y + z * b.z;
    }
    VecGosk cross(VecGosk b) const
    {
        return VecGosk(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }
    double norm() const
    {
        return sqrt(dot(*this));
    }
    double norm2() const
    {
        return dot(*this);
    }
    VecGosk neg() const
    {
        return VecGosk(-x, -y, -z);
    }
    VecGosk unit() const
    {
        double q = norm();
        if (q < 1e-300)
        {
            return VecGosk(0.0, 0.0, 0.0);
        }
        return scale(1.0 / q);
    }
    static VecGosk zero()
    {
        return VecGosk(0.0, 0.0, 0.0);
    }

    /* Переход в ECEF: x_e = x_g, y_e = -z_g, z_e = y_g. */
    VecEcef ecef() const
    {
        return VecEcef(x, -z, y);
    }

    Vec3 vec() const
    {
        return Vec3(x, y, z);
    }
    static VecGosk of(Vec3 v)
    {
        return VecGosk(v.x, v.y, v.z);
    }
};

/* Переход ECEF -> ГОСК: x_g = x_e, y_g = z_e, z_g = -y_e. */
inline VecGosk VecEcef::gosk() const
{
    return VecGosk(x, z, -y);
}

/* Тела методов сопряжения Mat3/Quat с фреймами (объявления -- в
 * lib_linalg.h). */
inline VecEnu Mat3::mulv(VecXBody v) const
{
    return VecEnu(m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
                  m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
                  m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z);
}
inline VecXBody Mat3::mulTv(VecEnu v) const
{
    return VecXBody(m[0][0] * v.e + m[1][0] * v.n + m[2][0] * v.u,
                    m[0][1] * v.e + m[1][1] * v.n + m[2][1] * v.u,
                    m[0][2] * v.e + m[1][2] * v.n + m[2][2] * v.u);
}
inline VecEnu Quat::rotate(VecXBody v) const
{
    return toDcm().mulv(v);
}
inline Quat Quat::deriv(VecXBody wBody) const
{
    return deriv(Vec3(wBody.x, wBody.y, wBody.z));
}

/* VecZBody -- связанная система, ось Z вдоль оси ракеты. */
struct VecZBody
{
    double x, y, z;

    VecZBody() : x(0.0), y(0.0), z(0.0)
    {
    }
    VecZBody(double x_, double y_, double z_) : x(x_), y(y_), z(z_)
    {
    }
};

/* Преобразования VecXBody <-> VecZBody: один и тот же физический
 * вектор в двух раскладках осей (на старте XBody = U,E,N; ZBody =
 * E,N,U): xb = (zb.z, zb.x, zb.y), zb = (xb.y, xb.z, xb.x). */
VecXBody xbodyFromZbody(const VecZBody& zb);
VecZBody zbodyFromXbody(const VecXBody& xb);

#endif /* LIB_FRAMES_H */

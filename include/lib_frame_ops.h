//============================================================================
// lib_frame_ops.h -- операции проекта над семействами систем координат gnc
// (V61: прежняя типизация VecEnu/VecXBody удалена, проект целиком на
// семействах именованных компонент cflib).
//
// Три группы средств:
//   1) конструкторы mk* -- присваивание ПО ИМЕНАМ полей (порядок
//      аргументов конструктора зафиксирован сигнатурой, позиционная
//      инициализация агрегатов gnc в проекте запрещена);
//   2) поэлементная алгебра add/sub/scale/neg/dot/cross/norm/unit/
//      clampNorm/lerp для CfLocalVector / CfBodyVector / CfLocalPoint
//      (точка +- вектор) -- те же формулы, что у прежних типизированных
//      векторов (порядок операций сохранён: результаты прежних прогонов
//      воспроизводимы);
//   3) АЛГЕБРАИЧЕСКИЙ мост кватерниона прежней пары (quatOldOf /
//      rotOfQuatOld) -- граница с ядром mrk6 и историческими записями;
//      от lib_linalg НЕ зависит. Матричная математика ESKF (V62:
//      единственный потребитель lib_linalg) -- src/eskf_math.h.
//
// Род величины связанной системы (вектор / угловая скорость / момент)
// задаётся типами gnc; переводы рода -- явные (omgVec/omgOf, momVec).
//============================================================================
#ifndef LIB_FRAME_OPS_H
#define LIB_FRAME_OPS_H

#include "lib_frames.h"
#include "lib_rotation.h"

#include <math.h>

//---------------------------------------------------------------------------
// Конструкторы по именам полей
//---------------------------------------------------------------------------
inline gnc::CfLocalVector mkLoc(double e, double n, double u)
{
    gnc::CfLocalVector r { };
    r.e = e;
    r.n = n;
    r.u = u;
    return r;
}

inline gnc::CfLocalPoint mkLop(double e, double n, double u)
{
    gnc::CfLocalPoint r { };
    r.e = e;
    r.n = n;
    r.u = u;
    return r;
}

inline gnc::CfBodyVector mkBody(double bu, double be, double bn)
{
    gnc::CfBodyVector r { };
    r.bu = bu;
    r.be = be;
    r.bn = bn;
    return r;
}

inline gnc::CfBodyOmega mkOmg(double bu, double be, double bn)
{
    gnc::CfBodyOmega r { };
    r.bu = bu;
    r.be = be;
    r.bn = bn;
    return r;
}

inline gnc::CfBodyMoment mkMom(double bu, double be, double bn)
{
    gnc::CfBodyMoment r { };
    r.bu = bu;
    r.be = be;
    r.bn = bn;
    return r;
}

inline gnc::CfEcefPoint mkEcp(double x, double y, double z)
{
    gnc::CfEcefPoint r { };
    r.x = x;
    r.y = y;
    r.z = z;
    return r;
}

inline gnc::CfEcefVector mkEcf(double x, double y, double z)
{
    gnc::CfEcefVector r { };
    r.x = x;
    r.y = y;
    r.z = z;
    return r;
}

inline gnc::CfGeodetic mkGdt(double lat, double lon, double alt)
{
    gnc::CfGeodetic g { };
    g.lat = lat;
    g.lon = lon;
    g.alt = alt;
    return g;
}

//---------------------------------------------------------------------------
// Переводы рода величины (явные) и точка <-> вектор
//---------------------------------------------------------------------------
inline gnc::CfBodyVector omgVec(const gnc::CfBodyOmega& w)
{
    return mkBody(w.bu, w.be, w.bn);
}

inline gnc::CfBodyOmega omgOf(const gnc::CfBodyVector& v)
{
    return mkOmg(v.bu, v.be, v.bn);
}

inline gnc::CfBodyVector momVec(const gnc::CfBodyMoment& m)
{
    return mkBody(m.bu, m.be, m.bn);
}

inline gnc::CfBodyMoment momOf(const gnc::CfBodyVector& v)
{
    return mkMom(v.bu, v.be, v.bn);
}

inline gnc::CfLocalVector vecOf(const gnc::CfLocalPoint& p)
{
    return mkLoc(p.e, p.n, p.u);
}

inline gnc::CfLocalPoint lopOf(const gnc::CfLocalVector& v)
{
    return mkLop(v.e, v.n, v.u);
}

//---------------------------------------------------------------------------
// Алгебра площадки: CfLocalVector (+ точка +- вектор)
//---------------------------------------------------------------------------
inline gnc::CfLocalVector add(const gnc::CfLocalVector& a,
                              const gnc::CfLocalVector& b)
{
    return mkLoc(a.e + b.e, a.n + b.n, a.u + b.u);
}

inline gnc::CfLocalVector sub(const gnc::CfLocalVector& a,
                              const gnc::CfLocalVector& b)
{
    return mkLoc(a.e - b.e, a.n - b.n, a.u - b.u);
}

inline gnc::CfLocalVector scale(const gnc::CfLocalVector& a, double s)
{
    return mkLoc(a.e * s, a.n * s, a.u * s);
}

inline gnc::CfLocalVector neg(const gnc::CfLocalVector& a)
{
    return mkLoc(-a.e, -a.n, -a.u);
}

inline double dot(const gnc::CfLocalVector& a, const gnc::CfLocalVector& b)
{
    return a.e * b.e + a.n * b.n + a.u * b.u;
}

inline gnc::CfLocalVector cross(const gnc::CfLocalVector& a,
                                const gnc::CfLocalVector& b)
{
    return mkLoc(a.n * b.u - a.u * b.n, a.u * b.e - a.e * b.u,
                 a.e * b.n - a.n * b.e);
}

inline double norm2(const gnc::CfLocalVector& a)
{
    return dot(a, a);
}

inline double norm(const gnc::CfLocalVector& a)
{
    return sqrt(dot(a, a));
}

inline gnc::CfLocalVector unit(const gnc::CfLocalVector& a)
{
    double q = norm(a);
    if (q < 1e-300)
    {
        return mkLoc(0.0, 0.0, 0.0);
    }
    return scale(a, 1.0 / q);
}

inline gnc::CfLocalVector clampNorm(const gnc::CfLocalVector& a,
                                    double maxNorm)
{
    double q = norm(a);
    if (q <= maxNorm || q < 1e-300)
    {
        return a;
    }
    return scale(a, maxNorm / q);
}

/* Линейная интерполяция a -> b, t в [0,1] (не клампится). */
inline gnc::CfLocalVector lerp(const gnc::CfLocalVector& a,
                               const gnc::CfLocalVector& b, double t)
{
    return add(a, scale(sub(b, a), t));
}

inline gnc::CfLocalPoint add(const gnc::CfLocalPoint& p,
                             const gnc::CfLocalVector& v)
{
    return mkLop(p.e + v.e, p.n + v.n, p.u + v.u);
}

inline gnc::CfLocalVector sub(const gnc::CfLocalPoint& a,
                              const gnc::CfLocalPoint& b)
{
    return mkLoc(a.e - b.e, a.n - b.n, a.u - b.u);
}

//---------------------------------------------------------------------------
// Алгебра связанной системы: CfBodyVector (+ угловая скорость)
//---------------------------------------------------------------------------
inline gnc::CfBodyVector add(const gnc::CfBodyVector& a,
                             const gnc::CfBodyVector& b)
{
    return mkBody(a.bu + b.bu, a.be + b.be, a.bn + b.bn);
}

inline gnc::CfBodyVector sub(const gnc::CfBodyVector& a,
                             const gnc::CfBodyVector& b)
{
    return mkBody(a.bu - b.bu, a.be - b.be, a.bn - b.bn);
}

inline gnc::CfBodyVector scale(const gnc::CfBodyVector& a, double s)
{
    return mkBody(a.bu * s, a.be * s, a.bn * s);
}

inline gnc::CfBodyVector neg(const gnc::CfBodyVector& a)
{
    return mkBody(-a.bu, -a.be, -a.bn);
}

inline double dot(const gnc::CfBodyVector& a, const gnc::CfBodyVector& b)
{
    return a.bu * b.bu + a.be * b.be + a.bn * b.bn;
}

/* Векторное произведение по циклике компонент (bu, be, bn) -- те же
 * формулы, что у прежнего типизированного вектора (x, y, z). */
inline gnc::CfBodyVector cross(const gnc::CfBodyVector& a,
                               const gnc::CfBodyVector& b)
{
    return mkBody(a.be * b.bn - a.bn * b.be, a.bn * b.bu - a.bu * b.bn,
                  a.bu * b.be - a.be * b.bu);
}

inline double norm2(const gnc::CfBodyVector& a)
{
    return dot(a, a);
}

inline double norm(const gnc::CfBodyVector& a)
{
    return sqrt(dot(a, a));
}

inline gnc::CfBodyVector unit(const gnc::CfBodyVector& a)
{
    double q = norm(a);
    if (q < 1e-300)
    {
        return mkBody(0.0, 0.0, 0.0);
    }
    return scale(a, 1.0 / q);
}

inline gnc::CfBodyVector clampNorm(const gnc::CfBodyVector& a, double maxNorm)
{
    double q = norm(a);
    if (q <= maxNorm || q < 1e-300)
    {
        return a;
    }
    return scale(a, maxNorm / q);
}

inline gnc::CfBodyOmega add(const gnc::CfBodyOmega& a,
                            const gnc::CfBodyOmega& b)
{
    return mkOmg(a.bu + b.bu, a.be + b.be, a.bn + b.bn);
}

inline gnc::CfBodyOmega scale(const gnc::CfBodyOmega& a, double s)
{
    return mkOmg(a.bu * s, a.be * s, a.bn * s);
}

/* Кватернион ПРЕЖНЕЙ ПАРЫ (Гамильтон, скаляр первым; на вертикальной
 * стойке (0.5, -0.5, -0.5, -0.5)) <-> CfRotation -- АЛГЕБРАИЧЕСКИЙ
 * мост без построения DCM и без Шепперда: знак детерминирован и
 * непрерывен вместе с кватернионом CfRotation.
 *
 * Вывод: q_old = q_L * qPark, где q_L -- кватернион CfRotation с
 * векторной частью, переименованной в оси (x, y, z) = (e, n, u)
 * прежней записи (q_L = (Qw, Qe, Qn, Qu)), qPark =
 * (0.5, -0.5, -0.5, -0.5) -- кватернион перестановки стойки
 * (тождество: quatOldOf(CfRotation()) = qPark; согласие с парой Quat
 * проверяется приёмкой lib_check). */
inline void quatOldOf(const gnc::CfRotation& r, double* w, double* x,
                      double* y, double* z)
{
    double aw = r.Qw();
    double ax = r.Qe();
    double ay = r.Qn();
    double az = r.Qu();
    /* произведение Гамильтона q_L * qPark, qPark = (0.5, -.5, -.5, -.5) */
    *w = 0.5 * (aw + ax + ay + az);
    *x = 0.5 * (-aw + ax - ay + az);
    *y = 0.5 * (-aw + ax + ay - az);
    *z = 0.5 * (-aw - ax + ay + az);
}

/* Обратный мост: CfRotation из кватерниона прежней пары.
 * q_L = q_old * qPark^-1, qPark^-1 = (0.5, 0.5, 0.5, 0.5);
 * далее SetQuat(w, an, ae, au) = (w_L, y_L, x_L, z_L). */
inline gnc::CfRotation rotOfQuatOld(double w, double x, double y, double z)
{
    double lw = 0.5 * (w - x - y - z);
    double lx = 0.5 * (x + w - z + y);
    double ly = 0.5 * (y + z + w - x);
    double lz = 0.5 * (z - y + x + w);
    gnc::CfRotation r;
    r.SetQuat(lw, ly, lx, lz);
    return r;
}

#endif // LIB_FRAME_OPS_H

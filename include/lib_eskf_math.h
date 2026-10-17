// gnclib 1.1.0: navigation ESKF (with the Earth rate, IMU-lever noise and
// gate recovery).
//============================================================================
// lib_eskf_math.h -- матричная математика фильтра навигации.
//
// ЕДИНСТВЕННОЕ место библиотеки, где используется безымянная линейная
// алгебра lib_linalg (Vec3/Mat3): якобианы пространства ошибок ESKF
// (skew, fromAxisAngle, произведения) записаны в матричной форме.
// Остальная библиотека -- на семействах gnc; сопряжение -- мосты
// этого заголовка (раскладка матрицы прежней пары: v_loc = R * v_body,
// строки -- e, n, u; столбцы -- bu, be, bn). Включается ТОЛЬКО из
// lib_eskf.h.
//============================================================================
#ifndef GNCLIB_LIB_ESKF_MATH_H
#define GNCLIB_LIB_ESKF_MATH_H

#include "lib_linalg.h"
#include "lib_frame_ops.h"

namespace gnc
{

// the global frame operations of lib_frame_ops.h join the overload sets of
// the gnc namespace (otherwise gnc::add(Vec3, ...) etc. hide them here)
using ::add;
using ::sub;
using ::scale;
using ::cross;
using ::dot;
using ::norm;
using ::norm2;
using ::omgVec;
using ::omgOf;
using ::momVec;
using ::momOf;
using ::mkBody;
using ::mkLoc;
using ::mkLop;
using ::mkOmg;
using ::mkMom;
using ::clampNorm;


//---------------------------------------------------------------------------
// Типизированные переходы по матрице прежней пары (формулы дословно
// прежних Mat3::mulv/mulTv c типизированными векторами)
//---------------------------------------------------------------------------
inline gnc::CfLocalVector mulv(const Mat3& R, const gnc::CfBodyVector& v)
{
    return mkLoc(R.m[0][0] * v.bu + R.m[0][1] * v.be + R.m[0][2] * v.bn,
                 R.m[1][0] * v.bu + R.m[1][1] * v.be + R.m[1][2] * v.bn,
                 R.m[2][0] * v.bu + R.m[2][1] * v.be + R.m[2][2] * v.bn);
}

//---------------------------------------------------------------------------
// Съём компонент в безымянный Vec3 (оси Родрига, якобианы)
//---------------------------------------------------------------------------
inline Vec3 vec3Of(const gnc::CfBodyVector& v)
{
    return Vec3(v.bu, v.be, v.bn);
}

// Вектор ESKF из вектора площадки gnc.
inline Vec3 vec3Of(const gnc::CfLocalVector& v)
{
    return Vec3(v.e, v.n, v.u);
}

//---------------------------------------------------------------------------
// Мосты Mat3 <-> gnc::CfRotation (точный перенос элементов DCM;
// Шепперд -- внутри SetDcm)
//---------------------------------------------------------------------------
inline gnc::CfRotation rotOfMat3(const Mat3& R)
{
    gnc::DCM d { };
    d.e_bu = R.m[0][0];
    d.e_be = R.m[0][1];
    d.e_bn = R.m[0][2];
    d.n_bu = R.m[1][0];
    d.n_be = R.m[1][1];
    d.n_bn = R.m[1][2];
    d.u_bu = R.m[2][0];
    d.u_be = R.m[2][1];
    d.u_bn = R.m[2][2];
    gnc::CfRotation r;
    r.SetDcm(d);
    return r;
}

// Матрица поворота 3x3 из gnc::CfRotation.
inline Mat3 mat3OfRot(const gnc::CfRotation& r)
{
    const gnc::DCM& d = r.Mat();
    Mat3 R;
    R.m[0][0] = d.e_bu;
    R.m[0][1] = d.e_be;
    R.m[0][2] = d.e_bn;
    R.m[1][0] = d.n_bu;
    R.m[1][1] = d.n_be;
    R.m[1][2] = d.n_bn;
    R.m[2][0] = d.u_bu;
    R.m[2][1] = d.u_be;
    R.m[2][2] = d.u_bn;
    return R;
}

//---------------------------------------------------------------------------
// Шимы Vec3: свободные функции в стиле выражений проекта
//---------------------------------------------------------------------------
inline Vec3 add(const Vec3& a, const Vec3& b)
{
    return a.add(b);
}

// Разность векторов.
inline Vec3 sub(const Vec3& a, const Vec3& b)
{
    return a.sub(b);
}

// Умножение вектора на скаляр.
inline Vec3 scale(const Vec3& a, double s)
{
    return a.scale(s);
}

// Скалярное произведение.
inline double dot(const Vec3& a, const Vec3& b)
{
    return a.dot(b);
}

// Векторное произведение.
inline Vec3 cross(const Vec3& a, const Vec3& b)
{
    return a.cross(b);
}

// Норма вектора.
inline double norm(const Vec3& a)
{
    return a.norm();
}

// Произведение матрицы на вектор R*v.
inline Vec3 mulv(const Mat3& R, const Vec3& v)
{
    return R.mulv(v);
}

}   // namespace gnc

#endif // GNCLIB_LIB_ESKF_MATH_H

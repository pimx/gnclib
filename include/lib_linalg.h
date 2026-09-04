/* lib_linalg.h -- минимальная 3D линейная алгебра для GNC.
 *
 * Vec3 / Mat3 / Quat; все операции -- методами (без перегрузки операторов).
 * Методы, возвращающие другой тип (Mat3::toQuat, Quat::toDcm), определены
 * вне классов после объявления всех трёх типов.
 * Без STL, без исключений, без динамических аллокаций.
 */
#ifndef LIB_LINALG_H
#define LIB_LINALG_H

#include <math.h>

struct Quat; /* опережающее объявление: используется в Mat3::toQuat */

/* ----------------------------- Vec3 ----------------------------- */
struct Vec3
{
    double x, y, z;

    Vec3() : x(0.0), y(0.0), z(0.0)
    {
    }
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_)
    {
    }

    Vec3 add(Vec3 b) const
    {
        return Vec3(x + b.x, y + b.y, z + b.z);
    }
    Vec3 sub(Vec3 b) const
    {
        return Vec3(x - b.x, y - b.y, z - b.z);
    }
    Vec3 scale(double s) const
    {
        return Vec3(x * s, y * s, z * s);
    }
    double dot(Vec3 b) const
    {
        return x * b.x + y * b.y + z * b.z;
    }
    Vec3 cross(Vec3 b) const
    {
        return Vec3(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }
    double norm() const
    {
        return sqrt(dot(*this));
    }

    /* Единичный вектор; нулевой вектор возвращается нулём. */
    Vec3 unit() const
    {
        double n = norm();
        if (n < 1e-300)
        {
            return Vec3(0, 0, 0);
        }
        return scale(1.0 / n);
    }

    /* Ограничение нормы сверху (полезно для клампов команд). */
    Vec3 clampNorm(double maxNorm) const
    {
        double n = norm();
        if (n <= maxNorm || n < 1e-300)
        {
            return *this;
        }
        return scale(maxNorm / n);
    }
};

/* ----------------------------- Mat3 ----------------------------- */
/* Построчное хранение 3x3: m[строка][столбец].
 * ТИПИЗАЦИЯ СИСТЕМ КООРДИНАТ -- семейства gnc (lib_frames.h,
 * lib_rotation.h); с этапа V61 прежние VecEnu/VecXBody удалены,
 * сопряжение типизированных величин -- gnc::CfRotation (B2L/L2B).
 * Vec3/Mat3/Quat -- безымянная линейная алгебра для внутренней
 * математики (фильтры, ковариации, интеграторы). */
struct Mat3
{
    double m[3][3];

    /* Умножение на вектор: r = M * v. */
    Vec3 mulv(Vec3 v) const
    {
        return Vec3(m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
                    m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
                    m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z);
    }

    /* Умножение матриц: R = (*this) * B. */
    Mat3 mul(Mat3 B) const
    {
        Mat3 R = Mat3::zero();
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                double s = 0.0;
                for (int k = 0; k < 3; k++)
                {
                    s += m[i][k] * B.m[k][j];
                }
                R.m[i][j] = s;
            }
        }
        return R;
    }
    Mat3 scale(double s) const
    {
        Mat3 R;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                R.m[i][j] = m[i][j] * s;
            }
        }
        return R;
    }
    Mat3 add(Mat3 B) const
    {
        Mat3 R;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                R.m[i][j] = m[i][j] + B.m[i][j];
            }
        }
        return R;
    }
    Mat3 sub(Mat3 B) const
    {
        Mat3 R;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                R.m[i][j] = m[i][j] - B.m[i][j];
            }
        }
        return R;
    }

    double trace() const
    {
        return m[0][0] + m[1][1] + m[2][2];
    }

    Vec3 row(int i) const
    {
        return Vec3(m[i][0], m[i][1], m[i][2]);
    }

    /* Переортонормирование DCM методом Грама-Шмидта (по строкам);
     * ограничивает накопление численного дрейфа. */
    Mat3 orthonormalize() const
    {
        Vec3 r0 = row(0);
        Vec3 r1 = row(1);
        r0 = r0.unit();
        r1 = r1.sub(r0.scale(r0.dot(r1)));
        r1 = r1.unit();
        Vec3 r2 = r0.cross(r1);
        Mat3 R;
        R.m[0][0] = r0.x;
        R.m[0][1] = r0.y;
        R.m[0][2] = r0.z;
        R.m[1][0] = r1.x;
        R.m[1][1] = r1.y;
        R.m[1][2] = r1.z;
        R.m[2][0] = r2.x;
        R.m[2][1] = r2.y;
        R.m[2][2] = r2.z;
        return R;
    }

    static Mat3 zero()
    {
        Mat3 R;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                R.m[i][j] = 0.0;
            }
        }
        return R;
    }
    static Mat3 eye()
    {
        Mat3 R = Mat3::zero();
        R.m[0][0] = R.m[1][1] = R.m[2][2] = 1.0;
        return R;
    }

    /* Кососимметричная матрица [w]x: skew(w).mulv(v) = w x v. */
    static Mat3 skew(Vec3 w)
    {
        Mat3 R = Mat3::zero();
        R.m[0][1] = -w.z;
        R.m[0][2] = w.y;
        R.m[1][0] = w.z;
        R.m[1][2] = -w.x;
        R.m[2][0] = -w.y;
        R.m[2][1] = w.x;
        return R;
    }

    /* Матрица поворота по оси (единичной) и углу: формула Родрига
     * R = I + sin(a) K + (1-cos(a)) K^2. */
    static Mat3 fromAxisAngle(Vec3 axisUnit, double angleRad)
    {
        Mat3 K = skew(axisUnit);
        Mat3 K2 = K.mul(K);
        return Mat3::eye()
            .add(K.scale(sin(angleRad)))
            .add(K2.scale(1.0 - cos(angleRad)));
    }

    Quat toQuat() const; /* определён вне класса (возвращает Quat) */
};

/* ----------------------------- Quat ----------------------------- */
/* Соглашение Гамильтона, скаляр первым: q = (w, x, y, z). */
struct Quat
{
    double w, x, y, z;

    Quat() : w(1.0), x(0.0), y(0.0), z(0.0)
    {
    }
    Quat(double w_, double x_, double y_, double z_)
        : w(w_), x(x_), y(y_), z(z_)
    {
    }

    double norm() const
    {
        return sqrt(w * w + x * x + y * y + z * z);
    }

    Quat normalize() const
    {
        double n = norm();
        if (n < 1e-300)
        {
            return Quat::id();
        }
        return Quat(w / n, x / n, y / n, z / n);
    }

    /* Произведение Гамильтона: q = (*this) * b. */
    Quat mul(Quat b) const
    {
        return Quat(w * b.w - x * b.x - y * b.y - z * b.z,
                    w * b.x + x * b.w + y * b.z - z * b.y,
                    w * b.y - x * b.z + y * b.w + z * b.x,
                    w * b.z + x * b.y - y * b.x + z * b.w);
    }

    /* Поворот вектора кватернионом (без построения DCM):
     * v' = v + 2 qv x (qv x v + w v). */
    Vec3 rotate(Vec3 v) const
    {
        Vec3 qv(x, y, z);
        Vec3 t = qv.cross(v).scale(2.0);
        return v.add(t.scale(w)).add(qv.cross(t));
    }

    /* Кватернион по оси (единичной) и углу. */
    static Quat fromAxisAngle(Vec3 axisUnit, double angleRad)
    {
        double h = 0.5 * angleRad;
        double s = sin(h);
        return Quat(cos(h), axisUnit.x * s, axisUnit.y * s, axisUnit.z * s);
    }

    static Quat id()
    {
        return Quat(1, 0, 0, 0);
    }

    Mat3 toDcm() const; /* определён вне класса (возвращает Mat3) */
};

/* --- методы, возвращающие другой тип (нужны оба полных типа) --- */

/* Кватернион (body->ref) в DCM: возвращает R такую, что v_ref = R * v_body. */
inline Mat3 Quat::toDcm() const
{
    Quat q = normalize();
    double qw = q.w, qx = q.x, qy = q.y, qz = q.z;
    Mat3 R;
    R.m[0][0] = 1 - 2 * (qy * qy + qz * qz);
    R.m[0][1] = 2 * (qx * qy - qw * qz);
    R.m[0][2] = 2 * (qx * qz + qw * qy);
    R.m[1][0] = 2 * (qx * qy + qw * qz);
    R.m[1][1] = 1 - 2 * (qx * qx + qz * qz);
    R.m[1][2] = 2 * (qy * qz - qw * qx);
    R.m[2][0] = 2 * (qx * qz - qw * qy);
    R.m[2][1] = 2 * (qy * qz + qw * qx);
    R.m[2][2] = 1 - 2 * (qx * qx + qy * qy);
    return R;
}

/* DCM (v_ref = R v_body) в кватернион, метод Шепперда (численно устойчив). */
inline Quat Mat3::toQuat() const
{
    double t = trace();
    Quat q;
    if (t > 0.0)
    {
        double s = sqrt(t + 1.0) * 2.0;
        q.w = 0.25 * s;
        q.x = (m[2][1] - m[1][2]) / s;
        q.y = (m[0][2] - m[2][0]) / s;
        q.z = (m[1][0] - m[0][1]) / s;
    }
    else if (m[0][0] > m[1][1] && m[0][0] > m[2][2])
    {
        double s = sqrt(1.0 + m[0][0] - m[1][1] - m[2][2]) * 2.0;
        q.w = (m[2][1] - m[1][2]) / s;
        q.x = 0.25 * s;
        q.y = (m[0][1] + m[1][0]) / s;
        q.z = (m[0][2] + m[2][0]) / s;
    }
    else if (m[1][1] > m[2][2])
    {
        double s = sqrt(1.0 + m[1][1] - m[0][0] - m[2][2]) * 2.0;
        q.w = (m[0][2] - m[2][0]) / s;
        q.x = (m[0][1] + m[1][0]) / s;
        q.y = 0.25 * s;
        q.z = (m[1][2] + m[2][1]) / s;
    }
    else
    {
        double s = sqrt(1.0 + m[2][2] - m[0][0] - m[1][1]) * 2.0;
        q.w = (m[1][0] - m[0][1]) / s;
        q.x = (m[0][2] + m[2][0]) / s;
        q.y = (m[1][2] + m[2][1]) / s;
        q.z = 0.25 * s;
    }
    return q.normalize();
}

#endif /* LIB_LINALG_H */

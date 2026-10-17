//============================================================================
// lib_util.h -- мелкие функции общего назначения (1.6.0): ограничение по
// модулю и диапазоном, перевод углов.
//
// В отличие от макроса clampd (lib_func.h) -- функции: аргумент
// вычисляется один раз. Подключается явно (не через lib.h): имена Limit,
// Clamp, Deg, ToDeg короткие, а потребители lib.h могут держать свои.
// C++17, без STL, без кучи, без исключений.
//============================================================================
#ifndef GNCLIB_LIB_UTIL_H
#define GNCLIB_LIB_UTIL_H

#include "lib_frames.h"

namespace gnc
{

// Ограничение по модулю: v в [-limit, limit].
//   v -- значение; limit -- предел (>= 0)
inline double Limit(double v, double limit)
{
    return v > limit ? limit : (v < -limit ? -limit : v);
}

// Ограничение диапазоном: v в [lo, hi].
//   v -- значение; lo, hi -- границы (lo <= hi)
inline double Clamp(double v, double lo, double hi)
{
    return v < lo ? lo : (v > hi ? hi : v);
}

// Градусы в радианы, выражение d * (PI / 180). Внимание: d * PI / 180
// округляется иначе (например, для 12 град отличие в 1 ЕМР).
//   d -- угол, град
constexpr double Deg(double d)
{
    return d * (PI / 180.0);
}

// Радианы в градусы.
//   r -- угол, рад
constexpr double ToDeg(double r)
{
    return r * (180.0 / PI);
}

} // namespace gnc

#endif // GNCLIB_LIB_UTIL_H

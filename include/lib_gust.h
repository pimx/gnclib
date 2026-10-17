//============================================================================
// lib_gust.h -- порывы ветра как функция времени (1.6.0).
//
// ПОЧЕМУ ФУНКЦИЯ ВРЕМЕНИ, А НЕ ФИЛЬТР. Модель интегрируется методом
// Рунге -- Кутты: производная запрашивается несколько раз внутри шага и не
// по порядку времени. Порыв как состояние фильтра (EnvWind, lib_wind.h)
// ответил бы на один и тот же момент по-разному. Здесь порыв -- ЧИСТАЯ
// функция времени: сумма GUST_TONES синусоид на каждую составляющую,
// периоды от tau до 20 tau, амплитуды ~ sqrt(период) (спектр, спадающий
// к высоким частотам), фазы -- хеш зерна, составляющей и номера тона.
// Нормировка: среднеквадратичное значение составляющей -- sigma.
// Горизонтальные составляющие независимы; вертикальная -- доля vertical
// горизонтальной.
// C++17, без STL, без кучи, без исключений.
//============================================================================
#ifndef GNCLIB_LIB_GUST_H
#define GNCLIB_LIB_GUST_H

#include <math.h>

#include "lib_frames.h"

namespace gnc
{

inline constexpr int GUST_TONES = 8; // синусоид на составляющую

// Все поля задаёт потребитель; при sigma = 0 порывов нет.
struct GustTones
{
    double sigma = 0.0;    // СКО горизонтального порыва, м/с; 0 -- порывов нет
    double tau = 1.0;      // время корреляции, с
    int seed = 0;          // фазы: одно зерно -- один и тот же воздух
    double vertical = 1.0; // вертикальный порыв в долях горизонтального
};

namespace detail
{

// Фаза тона по зерну, составляющей и номеру тона, рад.
inline double GustPhase(int seed, int component, int tone)
{
    unsigned int x = (unsigned int)(seed * 2654435761u) ^ (unsigned int)(component * 40503u)
                     ^ (unsigned int)(tone * 2246822519u);
    x = x ^ (x >> 13);
    x = x * 1274126177u;
    x = x ^ (x >> 16);
    return 2.0 * PI * ((double)(x % 100000u) / 100000.0);
}

// Одна составляющая единичного СКО в момент time.
inline double GustComponent(const GustTones& set, double time, int component)
{
    double sum = 0.0;
    double norm = 0.0;
    for (int k = 0; k < GUST_TONES; ++k)
    {
        double period = set.tau * (1.0 + 19.0 * (double)k / (double)(GUST_TONES - 1));
        double amp = sqrt(period / set.tau);
        double w = 2.0 * PI / period;
        sum = sum + amp * sin(w * time + GustPhase(set.seed, component, k));
        norm = norm + 0.5 * amp * amp;
    }
    return norm > 1.0e-12 ? sum / sqrt(norm) : 0.0;
}

} // namespace detail

// Порыв в момент времени.
//   set -- настройка; time -- время, с
// Возвращает скорость порыва в осях площадки, м/с (n, e -- горизонталь,
// u -- вертикаль); при sigma <= 0 или tau <= 0 -- ноль.
inline CfLocalVector GustTonesAt(const GustTones& set, double time)
{
    CfLocalVector g { };
    if (set.sigma <= 0.0 || set.tau <= 0.0)
    {
        return g;
    }
    g.n = set.sigma * detail::GustComponent(set, time, 0);
    g.e = set.sigma * detail::GustComponent(set, time, 1);
    g.u = set.sigma * set.vertical * detail::GustComponent(set, time, 2);
    return g;
}

} // namespace gnc

#endif // GNCLIB_LIB_GUST_H

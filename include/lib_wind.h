/* lib_wind.h -- модели ветра.
 *
 * class WindModel -- ветер в ENU: постоянная составляющая + порывы.
 * Порывы -- процесс Гаусса-Маркова 1-го порядка по каждой оси:
 *   dw/dt = -w/tau + sigma sqrt(2/tau) n(t)
 * (упрощённый аналог модели Драйдена, достаточный для малых скоростей;
 * вертикальные порывы ослаблены множителем 0.3).
 *
 * class WindProfile -- статические профили среднего ветра по высоте:
 * степенной и логарифмический (приземный слой).
 *
 * Библиотека не зависит от структур обмена проекта: вход/выход -- VecEnu и
 * секунды double. Без STL, без исключений, без динамических аллокаций.
 */
#ifndef LIB_WIND_H
#define LIB_WIND_H

#include "lib_frames.h"
#include "lib_rng.h"

class WindModel
{
public:
    WindModel();

    /* Настройка: средний ветер в ENU, СКО порывов [м/с], постоянная
     * времени порывов [с], зерно ГСЧ, стартовое время [с]. */
    void setup(VecEnu meanEnu, double gustSigma, double gustTau,
               unsigned long long seed, double t0);

    /* Продвижение модели порывов к моменту t [с]; возвращает текущую
     * скорость ветра в ENU. Точная дискретизация ОУ-процесса:
     * a = exp(-dt/tau), s = sigma sqrt(1 - a^2). */
    VecEnu update(double t);

    /* Текущая скорость ветра в ENU без продвижения, м/с. */
    VecEnu velocityEnu() const
    {
        return m_mean.add(m_gust);
    }

private:
    VecEnu m_mean;
    VecEnu m_gust;
    double m_sigma;
    double m_tau;
    Rng m_rng;
    double m_time;
};

class WindProfile
{
public:
    /* Степенной профиль: V(h) = vRef (h/hRef)^alpha; alpha ~0.14 для
     * открытой местности. h и hRef клампятся снизу 0.1 м. */
    static double powerLaw(double vRef, double hRef, double h, double alpha);

    /* Логарифмический профиль приземного слоя:
     * V(h) = vRef ln(h/z0)/ln(hRef/z0); z0 -- шероховатость [м].
     * Ниже z0 возврат 0. */
    static double logLaw(double vRef, double hRef, double h, double z0);
};

#endif /* LIB_WIND_H */

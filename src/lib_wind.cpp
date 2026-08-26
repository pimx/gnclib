#define _USE_MATH_DEFINES

/* lib_wind.cpp -- ветер: Гаусс-Марков порывы + профили по высоте. */
#include "lib_wind.h"

#include <math.h>

WindModel::WindModel()
    : m_mean(VecEnu::zero()),
      m_gust(VecEnu::zero()),
      m_sigma(0.0),
      m_tau(1.0),
      m_rng(123456789ULL),
      m_time(0.0)
{
}

/* Настройка модели: параметры порывов и зерно. */
void WindModel::setup(VecEnu meanEnu, double gustSigma, double gustTau,
                      unsigned long long seed, double t0)
{
    m_mean = meanEnu;
    m_sigma = gustSigma;
    m_tau = (gustTau > 1e-3) ? gustTau : 1e-3;
    m_gust = VecEnu::zero();
    m_rng.reseed(seed);
    m_time = t0;
}

/* Продвижение порывов на шаг: точная дискретизация процесса
 * Орнштейна-Уленбека; вертикальные порывы слабее (x0.3). */
VecEnu WindModel::update(double t)
{
    double dt = t - m_time;
    m_time = t;
    if (dt <= 0.0)
    {
        return velocityEnu();
    }
    double a = exp(-dt / m_tau);
    double s = m_sigma * sqrt(1.0 - a * a);
    m_gust.e = a * m_gust.e + s * m_rng.gauss();
    m_gust.n = a * m_gust.n + s * m_rng.gauss();
    m_gust.u = a * m_gust.u + s * (0.3 * m_rng.gauss());
    return velocityEnu();
}

/* Степенной профиль среднего ветра по высоте. */
double WindProfile::powerLaw(double vRef, double hRef, double h, double alpha)
{
    if (h < 0.1)
    {
        h = 0.1;
    }
    if (hRef < 0.1)
    {
        hRef = 0.1;
    }
    return vRef * pow(h / hRef, alpha);
}

/* Логарифмический профиль приземного слоя. */
double WindProfile::logLaw(double vRef, double hRef, double h, double z0)
{
    if (z0 < 1e-4)
    {
        z0 = 1e-4;
    }
    if (h <= z0)
    {
        return 0.0;
    }
    if (hRef <= z0)
    {
        hRef = z0 * 2.0;
    }
    return vRef * log(h / z0) / log(hRef / z0);
}

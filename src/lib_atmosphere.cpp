#define _USE_MATH_DEFINES

/* lib_atmosphere.cpp -- послойная модель US Standard Atmosphere 1976.
 * Высота трактуется как геометрическая (ошибка <0.5% ниже 50 км), что
 * достаточно для сопротивления и противодавления двигателя. */
#include "lib_atmosphere.h"

#include <math.h>

const double Atmosphere::P_SL = 101325.0;

static const double G0 = 9.80665; /* стандартная гравитация, м/с^2 */
static const double R_AIR =
    287.0528;                        /* газовая постоянная воздуха, Дж/(кг K) */
static const double R_E = 6356766.0; /* радиус для геопотенциальной высоты, м */

/* Конструктор: таблица слоёв (базовая высота, базовая температура,
 * градиент, базовое давление) модели US-76. */
Atmosphere::Atmosphere(double dTemp) : m_dTemp(dTemp)
{
    const Layer tab[N_LAYERS] = {
        {0.0, 288.15, -0.0065, 101325.0},
        {11000.0, 216.65, 0.0, 22632.06},
        {20000.0, 216.65, 0.001, 5474.889},
        {32000.0, 228.65, 0.0028, 868.0187},
        {47000.0, 270.65, 0.0, 110.9063},
        {51000.0, 270.65, -0.0028, 66.93887},
        {71000.0, 214.65, -0.002, 3.956420},
    };
    for (int i = 0; i < N_LAYERS; i++)
    {
        m_layers[i] = tab[i];
    }
}

/* Слой, чья база на уровне h или ниже (последний слой выше всех баз). */
int Atmosphere::layerIndex(double h) const
{
    int i = N_LAYERS - 1;
    for (int k = 0; k < N_LAYERS - 1; k++)
    {
        if (h < m_layers[k + 1].h_b)
        {
            i = k;
            break;
        }
    }
    return i;
}

/* Температура ISA: линейный градиент в слое; константа выше модели. */
double Atmosphere::temperatureIsa(double h) const
{
    if (h <= 0.0)
    {
        h = 0.0;
    }
    if (h >= 86000.0)
    {
        return 186.0;
    }
    int i = layerIndex(h);
    double T = m_layers[i].T_b + m_layers[i].L_b * (h - m_layers[i].h_b);
    if (T < 1e-3)
    {
        T = 1e-3;
    }
    return T;
}

/* Температура с учётом сдвига от ISA. */
double Atmosphere::temperature(double h) const
{
    double T = temperatureIsa(h) + m_dTemp;
    if (T < 1e-3)
    {
        T = 1e-3;
    }
    return T;
}

/* Давление: гидростатическое интегрирование в слое -- степенной закон при
 * ненулевом градиенте, экспонента для изотермического слоя. По соглашению
 * нестандартного дня давление НЕ сдвигается (профиль ISA). */
double Atmosphere::pressure(double h) const
{
    if (h <= 0.0)
    {
        h = 0.0;
    }
    if (h >= 86000.0)
    {
        return 0.0;
    }
    int i = layerIndex(h);
    double hb = m_layers[i].h_b, Tb = m_layers[i].T_b, Lb = m_layers[i].L_b,
           Pb = m_layers[i].P_b;
    double dh = h - hb;
    if (fabs(Lb) > 1e-12)
    {
        double T = Tb + Lb * dh;
        if (T < 1e-3)
        {
            T = 1e-3;
        }
        return Pb * pow(T / Tb, -G0 / (R_AIR * Lb));
    }
    return Pb * exp(-G0 * dh / (R_AIR * Tb));
}

/* Плотность из уравнения состояния: rho = P/(R T); T учитывает dTemp. */
double Atmosphere::density(double h) const
{
    if (h <= 0.0)
    {
        h = 0.0;
    }
    if (h >= 86000.0)
    {
        return 0.0;
    }
    double T = temperature(h);
    double P = pressure(h);
    return P / (R_AIR * T);
}

/* Отношение давления к уровню моря, кламп [0,1]. */
double Atmosphere::pressureRatio(double h) const
{
    double pr = pressure(h) / P_SL;
    if (pr > 1.0)
    {
        pr = 1.0;
    }
    if (pr < 0.0)
    {
        pr = 0.0;
    }
    return pr;
}

/* Скорость звука для gamma = 1.4 (идеальный двухатомный воздух). */
double Atmosphere::soundSpeed(double h) const
{
    return sqrt(1.4 * R_AIR * temperature(h));
}

/* Число Маха; 0 при вырождении скорости звука. */
double Atmosphere::mach(double h, double V) const
{
    double a = soundSpeed(h);
    return (a > 1e-6) ? V / a : 0.0;
}

/* Скоростной напор q = 1/2 rho V^2. */
double Atmosphere::dynamicPressure(double h, double V) const
{
    return 0.5 * density(h) * V * V;
}

/* Динамическая вязкость по формуле Сазерленда. */
double Atmosphere::dynamicViscosity(double h) const
{
    double T = temperature(h);
    return 1.458e-6 * T * sqrt(T) / (T + 110.4);
}

/* Кинематическая вязкость nu = mu/rho; 0 при нулевой плотности. */
double Atmosphere::kinematicViscosity(double h) const
{
    double rho = density(h);
    if (rho <= 0.0)
    {
        return 0.0;
    }
    return dynamicViscosity(h) / rho;
}

/* Число Рейнольдса Re = V L / nu; 0 при вырождении вязкости. */
double Atmosphere::reynolds(double h, double V, double L) const
{
    double nu = kinematicViscosity(h);
    if (nu <= 0.0)
    {
        return 0.0;
    }
    return V * L / nu;
}

/* Геометрическая -> геопотенциальная высота. */
double Atmosphere::geometricToGeopotential(double h)
{
    return R_E * h / (R_E + h);
}

/* Геопотенциальная -> геометрическая высота. */
double Atmosphere::geopotentialToGeometric(double H)
{
    return R_E * H / (R_E - H);
}

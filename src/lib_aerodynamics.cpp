#define _USE_MATH_DEFINES

/* lib_aerodynamics.cpp -- модели аэродинамических сил (сопротивление,
 * полосовая модель цилиндра, рули, экранный эффект). */
#include "lib_aerodynamics.h"

#include <math.h>

Aerodynamics::Aerodynamics(const Atmosphere& atmo) : m_atmo(atmo)
{
}

/* Сопротивление с постоянным Cd против воздушной скорости. Ноль при
 * вырождении потока или плотности. */
// Скалярное ядро: модуль силы сопротивления с постоянным Cd.
// Возврат 0 при пренебрежимой скорости или нулевой плотности.
double Aerodynamics::dragMagnitude(double h, double V, double Cd,
                                   double aRef) const
{
    double rho = m_atmo.density(h);
    if (V < 1e-6 || rho <= 0.0)
    {
        return 0.0;
    }
    return 0.5 * rho * V * V * Cd * aRef;
}

// Перегрузки по системе координат: сила направлена ПРОТИВ воздушной
// скорости, поэтому результат остаётся в системе аргумента.
VecXBody Aerodynamics::dragConstantCd(double h, VecXBody vAir, double Cd,
                                     double aRef) const
{
    double V = vAir.norm();
    double Fmag = dragMagnitude(h, V, Cd, aRef);
    if (Fmag <= 0.0)
    {
        return VecXBody::zero();
    }
    return vAir.scale(-Fmag / V);
}

VecEnu Aerodynamics::dragConstantCd(double h, VecEnu vAir, double Cd,
                                    double aRef) const
{
    double V = vAir.norm();
    double Fmag = dragMagnitude(h, V, Cd, aRef);
    if (Fmag <= 0.0)
    {
        return VecEnu::zero();
    }
    return vAir.scale(-Fmag / V);
}

/* Посадочный Cx(M): трансзвуковой горб на сверхзвуковом плато. */
double Aerodynamics::cxMachLanding(double M)
{
    /* Табличное задание Cx(M) лобового сопротивления (тупой торец,
     * посадочная конфигурация): узлы -- линейная интерполяция, за
     * пределами таблицы -- ближайший узел. Таблица получена
     * дискретизацией прежней аналитической аппроксимации и подлежит
     * замене продувочными данными без изменения кода. */
    static const double tabM[] = {0.0,  0.3,  0.6,  0.8,  0.9,  1.0,
                                  1.1,  1.2,  1.4,  1.6,  2.0,  2.5,
                                  3.0,  4.0,  5.0,  6.0,  10.0};
    static const double tabCx[] = {1.5212, 1.5102, 1.5506, 1.6460, 1.6870,
                                   1.6964, 1.6653, 1.6029, 1.4659, 1.3886,
                                   1.3332, 1.2940, 1.2715, 1.2546, 1.2509,
                                   1.2502, 1.2500};
    const int n = (int)(sizeof(tabM) / sizeof(tabM[0]));
    if (M <= tabM[0])
    {
        return tabCx[0];
    }
    if (M >= tabM[n - 1])
    {
        return tabCx[n - 1];
    }
    int i = 0;
    while (i < n - 2 && M > tabM[i + 1])
    {
        ++i;
    }
    double f = (M - tabM[i]) / (tabM[i + 1] - tabM[i]);
    return tabCx[i] + (tabCx[i + 1] - tabCx[i]) * f;
}

/* Маховская поправка осевого Cx: единица ниже M = 0.3, выше -- рост по
 * закону cxMachLanding, нормированному к дозвуковому значению. */
double Aerodynamics::machFactorAxial(double M)
{
    if (M <= 0.3)
    {
        return 1.0;
    }
    double f = cxMachLanding(M) / cxMachLanding(0.3);
    return (f > 1.0) ? f : 1.0;
}

/* Экранирование ретро-струёй: линейное снижение с долей тяги, пол 5%
 * сохраняет корректность интегрирования. */
double Aerodynamics::cxRetroPlume(double Cx, double k, double thrFrac)
{
    if (k <= 0.0 || thrFrac <= 0.0)
    {
        return Cx;
    }
    if (thrFrac > 1.0)
    {
        thrFrac = 1.0;
    }
    double f = 1.0 - k * thrFrac;
    if (f < 0.05)
    {
        f = 0.05;
    }
    return Cx * f;
}

/* Поправка сжимаемости Прандтля-Глауэрта: 1/sqrt(1 - M^2), кламп M 0.95. */
double Aerodynamics::prandtlGlauert(double M)
{
    if (M < 0.0)
    {
        M = 0.0;
    }
    if (M > 0.95)
    {
        M = 0.95;
    }
    return 1.0 / sqrt(1.0 - M * M);
}

/* Экранный эффект по Чизмену-Беннетту: T_ige/T_oge = 1/(1-(R/4h)^2). */
double Aerodynamics::groundEffectFactor(double diskRadius, double h,
                                        double maxFactor)
{
    if (diskRadius <= 0.0)
    {
        return 1.0;
    }
    if (maxFactor < 1.0)
    {
        maxFactor = 1.0;
    }
    double quarter = 0.25 * diskRadius;
    if (h <= quarter)
    {
        return maxFactor;
    }
    double r = quarter / h;
    double f = 1.0 / (1.0 - r * r);
    if (f > maxFactor)
    {
        f = maxFactor;
    }
    if (f < 1.0)
    {
        f = 1.0;
    }
    return f;
}

/* Посадочная сила сопротивления: Cx(M) со снижением струёй, против
 * воздушной скорости. */
// Скалярное ядро посадочного сопротивления: Cx(M) со снижением струёй.
double Aerodynamics::dragLandingMagnitude(double h, double V, double aRef,
                                          double retroPlumeK,
                                          double thrFrac) const
{
    double rho = m_atmo.density(h);
    if (V < 1e-6 || rho <= 0.0)
    {
        return 0.0;
    }
    double M = m_atmo.mach(h, V);
    double Cx = cxRetroPlume(cxMachLanding(M), retroPlumeK, thrFrac);
    return 0.5 * rho * V * V * Cx * aRef;
}

VecXBody Aerodynamics::dragLanding(double h, VecXBody vAir, double aRef,
                                  double retroPlumeK, double thrFrac) const
{
    double V = vAir.norm();
    double Fmag = dragLandingMagnitude(h, V, aRef, retroPlumeK, thrFrac);
    if (Fmag <= 0.0)
    {
        return VecXBody::zero();
    }
    return vAir.scale(-Fmag / V);
}

VecEnu Aerodynamics::dragLanding(double h, VecEnu vAir, double aRef,
                                 double retroPlumeK, double thrFrac) const
{
    double V = vAir.norm();
    double Fmag = dragLandingMagnitude(h, V, aRef, retroPlumeK, thrFrac);
    if (Fmag <= 0.0)
    {
        return VecEnu::zero();
    }
    return vAir.scale(-Fmag / V);
}

/* Скалярное посадочное замедление для точечных прогнозов. */
double Aerodynamics::dragDecelLanding(double h, double V, double aRef, double m,
                                      double retroPlumeK, double thrFrac) const
{
    if (m < 1.0)
    {
        m = 1.0;
    }
    double rho = m_atmo.density(h);
    if (V < 1e-6 || rho <= 0.0)
    {
        return 0.0;
    }
    double Cx =
        cxRetroPlume(cxMachLanding(m_atmo.mach(h, V)), retroPlumeK, thrFrac);
    return 0.5 * rho * V * V * Cx * aRef / m;
}

/* Потенциал пары рулей на полном отклонении: 2 q Cl_delta A_fin delta_max. */
double Aerodynamics::finForceMax(double h, double V, double finArea,
                                 double finClDelta, double finMaxDeg) const
{
    double q = m_atmo.dynamicPressure(h, V);
    double dmax = finMaxDeg * M_PI / 180.0;
    return 2.0 * q * finClDelta * finArea * dmax;
}

/* Полосовая модель цилиндра: поперечное обтекание по полосам с учётом
 * вращения + осевое сопротивление по миделю (маховская поправка и
 * экранирование ретро-струёй на спуске). */
void Aerodynamics::cylinderStrip(VecXBody vRelBody, VecXBody omegaBody,
                                 double rho, double mach, double thrFrac,
                                 double length, double radius,
                                 double cmFromTail, double cdLateral,
                                 double cdAxial, double retroPlumeK,
                                 int nStrips, VecXBody* forceBody,
                                 VecXBody* momentBody) const
{
    *forceBody = VecXBody::zero();
    *momentBody = VecXBody::zero();
    if (rho <= 0.0 || length <= 0.0 || radius <= 0.0)
    {
        return;
    }
    if (nStrips < 1)
    {
        nStrips = 1;
    }
    if (nStrips > 32)
    {
        nStrips = 32;
    }

    double dx = length / (double)nStrips;
    double dia = 2.0 * radius;

    /* Поперечное обтекание по полосам. Ось тела -- X (VecXBody: x вдоль
     * ракеты к носу, y и z поперечные); полосы расположены вдоль x,
     * поперечная составляющая местной скорости -- (0, y, z). */
    for (int i = 0; i < nStrips; ++i)
    {
        double xi = (i + 0.5) * dx - cmFromTail; /* полоса относительно ЦМ */
        VecXBody r(xi, 0.0, 0.0);
        VecXBody vLoc = vRelBody.add(omegaBody.cross(r));
        VecXBody vLat(0.0, vLoc.y, vLoc.z); /* поперечная составляющая */
        double vm = vLat.norm();
        if (vm < 1e-9)
        {
            continue;
        }
        VecXBody dF = vLat.scale(-0.5 * rho * cdLateral * dia * dx * vm);
        *forceBody = forceBody->add(dF);
        *momentBody = momentBody->add(r.cross(dF));
    }

    /* Осевое сопротивление (площадь миделя): маховская поправка и
     * экранирование ретро-струёй (струя навстречу потоку -- спуск, то есть
     * воздушная скорость вдоль оси отрицательна: поток набегает на хвост). */
    double aAx = M_PI * radius * radius;
    double vx = vRelBody.x;
    double cdAx = cdAxial * machFactorAxial(mach);
    if (vx < 0.0) /* спуск: струя против потока */
    {
        cdAx = cxRetroPlume(cdAx, retroPlumeK, thrFrac);
    }
    forceBody->x += -0.5 * rho * cdAx * aAx * fabs(vx) * vx;
}

//============================================================================
// lib_slosh_params.cpp -- параметры аналога плескания (1.7.0): строки
// цилиндра с плоским днищем, таблицы по заполнению, перевод строки в
// параметры бака, оценки демпфирования. См. lib_slosh.h.
//============================================================================
#include "lib_slosh.h"

#include <math.h>

#include "lib_table.h"
#include "lib_util.h"

namespace gnc
{

namespace
{

const int NSERIES = 100;        // членов ряда импульсной инерции

// первые нули J1'(x)
const double XI_TAB[10] =
{
    1.8411837813, 5.3314427735, 8.5363163664, 11.7060049026, 14.8635886339,
    18.0155278626, 21.1643698588, 24.3113268572, 27.4570505710, 30.6019229726
};

double Sq(double x)
{
    return x * x;
}

double Lerp(double a, double b, double t)
{
    return a + (b - a) * t;
}

} // namespace

//----------------------------------------------------------------------------
// n-й нуль J1'(x). При n > 10 -- разложение Макмагона
//   x = b - (mu + 3) / (8 b) - 4 (7 mu^2 + 82 mu - 9) / (3 (8 b)^3),
//   b = (n - 1/4) pi, mu = 4;
// ошибка меньше 5e-8.
//   n -- номер нуля, 1, 2, ...
//----------------------------------------------------------------------------
double SloshBesselZero(int n)
{
    if (n < 1)
    {
        return 0.0;
    }
    if (n <= 10)
    {
        return XI_TAB[n - 1];
    }
    double b = (n - 0.25) * PI;
    double b8 = 8.0 * b;
    return b - 7.0 / b8 - 4.0 * 431.0 / (3.0 * b8 * b8 * b8);
}

//----------------------------------------------------------------------------
// Приведённый аналог цилиндра с плоским днищем (см. lib_slosh.h).
// Алгоритм: тона n = 1..NSERIES решения для потенциального течения;
// первые nModes -- явные, остальные присоединяются к неподвижной массе.
// Импульсная инерция жидкости относительно её ЦМ (потенциал
// Стокса-Жуковского поворота по тангажу с p = 0 на свободной поверхности,
// члены стенки и днища ряда Бесселя-Дини в форме без переполнения -- через
// th и 1/ch):
//   I/rho = pi R (R h^3/12 + sum g (2T/k^3 - h^2 T/(4k) - h/(2k^2)
//           - 3h/(2 C k^2))) + pi (-h R^4/8 + sum g (2T/k - h/(2C)) R^3/xi^2),
//   g = 2R/(xi^2 - 1), k = xi/R, T = th(k h), C = ch(k h).
// I0 относительно неподвижной массы = I + sum(присоединённые m_n l_n^2)
// - m0 l0^2; плечо l0 неподвижной массы -- из ЦМ жидкости. Ряд сверен с
// независимым решением методом конечных разностей (HISTORY, 1.7.0).
//   radius, depth -- R, h, м; nModes -- 1..3; row -- выход
//----------------------------------------------------------------------------
bool SloshCylinderRow(double radius, double depth, int nModes, SloshRow* row)
{
    if (row == 0 || radius <= 0.0 || depth <= 0.0 || nModes < 1)
    {
        return false;
    }
    int nm = (nModes > SLOSH_NMODE) ? SLOSH_NMODE : nModes;
    SloshRow r;
    double R = radius;
    double h = depth;
    r.zSurf = h;
    r.zLiq = 0.5 * h;
    r.iRollr = 0.5 * R * R;
    double sumWall = 0.0;
    double sumBot = 0.0;
    double sumM = 0.0;
    double sumML = 0.0;
    double sumLump = 0.0;
    for (int n = 1; n <= NSERIES; ++n)
    {
        double xi = SloshBesselZero(n);
        double k = xi / R;
        double kh = k * h;
        double t = tanh(kh);
        double iC = (kh > 700.0) ? 0.0 : 1.0 / cosh(kh);
        double g = 2.0 * R / (xi * xi - 1.0);
        sumWall += g * (2.0 * t / (k * k * k) - h * h * t / (4.0 * k) - h / (2.0 * k * k)
                        - 1.5 * h * iC / (k * k));
        sumBot += g * (2.0 * t / k - 0.5 * h * iC) * R * R * R / (xi * xi);
        double mr = 2.0 * R * t / (xi * (xi * xi - 1.0) * h);
        double ln = 0.5 * h - (2.0 * R / xi) * tanh(0.5 * kh);
        if (n <= nm)
        {
            r.mr[n - 1] = mr;
            r.z[n - 1] = 0.5 * h + ln;
            r.w2a[n - 1] = k * t;
            r.eta[n - 1] = 2.0 * xi * t / (xi * xi - 1.0);
            sumM += mr;
            sumML += mr * ln;
        }
        else
        {
            sumLump += mr * ln * ln;
        }
    }
    double iImp = (R * (R * h * h * h / 12.0 + sumWall) - h * R * R * R * R / 8.0 + sumBot)
                  / (R * R * h);
    r.m0r = 1.0 - sumM;
    double l0 = (r.m0r > 1.0e-12) ? -sumML / r.m0r : 0.0;
    r.z0 = 0.5 * h + l0;
    r.i0r = iImp + sumLump - r.m0r * l0 * l0;
    *row = r;
    return true;
}

//----------------------------------------------------------------------------
// Таблица цилиндра: q = qa + (qb - qa) s^2, s = 0..1 -- гуще при малом
// заполнении, где th(xi h / R) меняется быстрее всего; qa = 0.005,
// qb = 0.995.
//   radius, length -- м; nModes -- 1..3; nRows -- 2..SLOSH_NROW; tab -- выход
//----------------------------------------------------------------------------
bool SloshTableCylinder(double radius, double length, int nModes, int nRows, SloshTable* tab)
{
    if (tab == 0 || length <= 0.0 || nRows < 2 || nRows > SLOSH_NROW)
    {
        return false;
    }
    const double qa = 0.005;
    const double qb = 0.995;
    tab->n = 0;
    for (int i = 0; i < nRows; ++i)
    {
        double s = (double)i / (double)(nRows - 1);
        double q = qa + (qb - qa) * s * s;
        if (!SloshCylinderRow(radius, q * length, nModes, &tab->row[i]))
        {
            tab->n = 0;
            return false;
        }
        tab->q[i] = q;
        tab->n = i + 1;
    }
    return true;
}

//----------------------------------------------------------------------------
// Строка при q: интервал -- TableLocate (края удерживаются), каждое поле --
// линейно.
//   tab -- таблица; q -- доля заполнения; row -- выход
//----------------------------------------------------------------------------
void SloshTableRow(const SloshTable& tab, double q, SloshRow* row)
{
    if (tab.n <= 0)
    {
        SloshRow z;
        *row = z;
        return;
    }
    if (tab.n == 1)
    {
        *row = tab.row[0];
        return;
    }
    double t = 0.0;
    int i = TableLocate(tab.q, tab.n, q, &t);
    const SloshRow& a = tab.row[i];
    const SloshRow& b = tab.row[i + 1];
    SloshRow r;
    r.zSurf = Lerp(a.zSurf, b.zSurf, t);
    r.zLiq = Lerp(a.zLiq, b.zLiq, t);
    r.m0r = Lerp(a.m0r, b.m0r, t);
    r.z0 = Lerp(a.z0, b.z0, t);
    r.i0r = Lerp(a.i0r, b.i0r, t);
    r.iRollr = Lerp(a.iRollr, b.iRollr, t);
    for (int n = 0; n < SLOSH_NMODE; ++n)
    {
        r.mr[n] = Lerp(a.mr[n], b.mr[n], t);
        r.z[n] = Lerp(a.z[n], b.z[n], t);
        r.w2a[n] = Lerp(a.w2a[n], b.w2a[n], t);
        r.eta[n] = Lerp(a.eta[n], b.eta[n], t);
    }
    *row = r;
}

//----------------------------------------------------------------------------
// Строка -> параметры бака (см. lib_slosh.h). Присоединение тонов
// n >= nModes: неподвижная масса и точечные массы объединяются по теореме
// о параллельных осях относительно общего ЦМ -- ЦМ жидкости сохраняется.
//   tank -- описание; row -- аналог; mL -- кг; a -- м/с^2 (со знаком);
//   aMin, bondMin -- пределы; p -- выход
//----------------------------------------------------------------------------
void SloshRowToParams(const SloshTank& tank, const SloshRow& row, double mL, double a,
                      double aMin, double bondMin, SloshTankParams* p)
{
    SloshTankParams r;
    int nm = (int)Clamp(tank.nModes, 1, SLOSH_NMODE);
    double m = (mL > 0.0) ? mL : 0.0;
    double sgn = (a >= 0.0) ? 1.0 : -1.0;
    double aAbs = fabs(a);
    r.valid = (m > 0.0) ? 1 : 0;
    if (aAbs < aMin)
    {
        aAbs = aMin;
        r.valid = 0;
    }
    r.a = sgn * aAbs;
    r.sgn = sgn;
    double mFull = tank.massFull;
    if (mFull <= 0.0)
    {
        mFull = tank.rho * PI * Sq(tank.radius) * tank.length;
    }
    r.mL = m;
    r.q = (mFull > 0.0) ? m / mFull : 0.0;
    double endArm = (sgn > 0.0) ? tank.armAft : tank.armAft + tank.length;
    double mSum = row.m0r;
    double mz = row.m0r * row.z0;
    for (int n = nm; n < SLOSH_NMODE; ++n)
    {
        mSum += row.mr[n];
        mz += row.mr[n] * row.z[n];
    }
    double z0 = (mSum > 1.0e-12) ? mz / mSum : row.zLiq;
    double i0 = row.i0r + row.m0r * Sq(row.z0 - z0);
    for (int n = nm; n < SLOSH_NMODE; ++n)
    {
        i0 += row.mr[n] * Sq(row.z[n] - z0);
    }
    r.m0 = mSum * m;
    r.arm0 = endArm + sgn * z0;
    r.i0 = i0 * m;
    r.iRoll = tank.kRoll * row.iRollr * m;
    r.armSurf = endArm + sgn * row.zSurf;
    r.armLiq = endArm + sgn * row.zLiq;
    double zw = (tank.nu > 0.0) ? SloshZetaWall(tank.nu, aAbs, tank.radius, row.zSurf) : 0.0;
    r.nModes = nm;
    for (int n = 0; n < nm; ++n)
    {
        SloshModeParams& md = r.mode[n];
        md.m = row.mr[n] * m;
        md.w = sqrt(aAbs * row.w2a[n]);
        md.k = md.m * md.w * md.w;
        md.arm = endArm + sgn * row.z[n];
        md.lPend = (row.w2a[n] > 0.0) ? 1.0 / row.w2a[n] : 0.0;
        md.eta = row.eta[n];
        md.zeta0 = (tank.zetaFix[n] > 0.0) ? tank.zetaFix[n] : zw;
    }
    if (tank.sigma > 0.0)
    {
        r.bond = tank.rho * aAbs * Sq(tank.radius) / tank.sigma;
        if (r.bond < bondMin)
        {
            r.valid = 0;
        }
    }
    *p = r;
}

//----------------------------------------------------------------------------
// Демпфирование гладкой стенкой (см. lib_slosh.h); h/R не меньше 0.05.
//   nu -- м^2/с; a -- м/с^2; radius, depth -- м
//----------------------------------------------------------------------------
double SloshZetaWall(double nu, double a, double radius, double depth)
{
    if (nu <= 0.0 || a <= 0.0 || radius <= 0.0)
    {
        return 0.0;
    }
    double hr = depth / radius;
    hr = (hr < 0.05) ? 0.05 : hr;
    double x = 1.84 * hr;
    double br = 1.0 + 0.318 / sinh(x) * (1.0 + (1.0 - hr) / cosh(x));
    return 0.79 * sqrt(nu) * pow(a, -0.25) * pow(radius, -0.75) * br;
}

//----------------------------------------------------------------------------
// Демпфирование кольцевой перегородкой по Майлсу (см. lib_slosh.h).
//   depth -- заглубление пластины, м; width -- ширина кольца, м;
//   radius -- м; wave -- высота волны у стенки, м
//----------------------------------------------------------------------------
double SloshZetaBaffle(double depth, double width, double radius, double wave)
{
    if (depth <= 0.0 || width <= 0.0 || radius <= 0.0 || wave <= 0.0)
    {
        return 0.0;
    }
    double w = Clamp(width / radius, 0.0, 1.0);
    double beta = 1.0 - (1.0 - w) * (1.0 - w);
    return 2.83 * exp(-4.6 * depth / radius) * beta * sqrt(beta) * sqrt(wave / radius);
}

} // namespace gnc

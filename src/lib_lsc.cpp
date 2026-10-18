//============================================================================
// lib_lsc.cpp -- наименьшие квадраты с границами переменных и линейными
// ограничениями на пары переменных (1.7.0), см. lib_lsc.h.
//
// Группа -- одиночная переменная либо пара. Ограничения группы: границы её
// переменных и (у пары) ограничения пары; номера: 2 j -- нижняя граница
// переменной j (-x_j <= -lb_j), 2 j + 1 -- верхняя (x_j <= ub_j), 2 n + c --
// ограничение пары c. Рабочий набор группы задаёт её параметризацию:
// пусто -- переменные свободны; одно ограничение пары с нормалью a --
// x = a h + t s, t = (-a2, a1), s -- одна свободная переменная; два (либо
// одно у одиночной) -- группа закреплена.
//============================================================================
#include "lib_lsc.h"

#include <math.h>
#include <string.h>

namespace gnc
{

Lsc::Lsc()
{
    m_n = 0;
    m_ng = 0;
    m_nc = 0;
    memset(m_wn, 0, sizeof(m_wn));
}

bool Lsc::setup(int n, const double* lb, const double* ub)
{
    m_n = 0;
    m_nc = 0;
    if (n < 1 || n > LSC_MAX_COLS)
    {
        return false;
    }
    m_n = n;
    m_ng = n;
    for (int j = 0; j < n; j++)
    {
        m_lb[j] = lb[j];
        m_ub[j] = ub[j] > lb[j] ? ub[j] : lb[j];
        m_fix[j] = (m_ub[j] - m_lb[j] <= 1e-15 * (1.0 + fabs(m_lb[j]))) ? 1 : 0;
        m_grp[j] = j;
        m_gv[j][0] = j;
        m_gv[j][1] = -1;
        m_wn[j] = 0;
        m_first[j] = -1;
        m_last[j] = -1;
    }
    return true;
}

int Lsc::pair(int j, int k)
{
    if (j < 0 || k < 0 || j >= m_n || k >= m_n || j == k)
    {
        return -1;
    }
    if (m_grp[j] != j || m_grp[k] != k || m_gv[j][1] >= 0 || m_gv[k][1] >= 0)
    {
        return -1;
    }
    m_gv[j][1] = k;
    m_gv[k][0] = -1;
    m_grp[k] = j;
    return j;
}

bool Lsc::cut(int p, double a1, double a2, double h)
{
    if (p < 0 || p >= m_n || m_gv[p][0] != p || m_gv[p][1] < 0 || m_nc >= LSC_MAX_CUTS)
    {
        return false;
    }
    double s = sqrt(a1 * a1 + a2 * a2);
    if (s <= 1e-300)
    {
        return false;
    }
    int c = m_nc;
    m_ca[c][0] = a1 / s;
    m_ca[c][1] = a2 / s;
    m_ch[c] = h / s;
    m_cg[c] = p;
    m_next[c] = -1;
    if (m_last[p] >= 0)
    {
        m_next[m_last[p]] = c;
    }
    else
    {
        m_first[p] = c;
    }
    m_last[p] = c;
    m_nc++;
    return true;
}

int Lsc::activeCount() const
{
    int c = 0;
    for (int g = 0; g < m_n; g++)
    {
        if (m_gv[g][0] < 0)
        {
            continue;
        }
        for (int s = 0; s < m_wn[g]; s++)
        {
            int id = m_ws[g][s];
            c += (id < 2 * m_n && m_fix[id / 2]) ? 0 : 1;
        }
    }
    return c;
}

void Lsc::row(int id, double* a0, double* a1, double* h) const
{
    *a0 = 0.0;
    *a1 = 0.0;
    if (id >= 2 * m_n)
    {
        int c = id - 2 * m_n;
        *a0 = m_ca[c][0];
        *a1 = m_ca[c][1];
        *h = m_ch[c];
        return;
    }
    int j = id / 2;
    int g = m_grp[j];
    double* a = m_gv[g][0] == j ? a0 : a1;
    if ((id & 1) == 0)
    {
        *a = -1.0;
        *h = -m_lb[j];
    }
    else
    {
        *a = 1.0;
        *h = m_ub[j];
    }
}

bool Lsc::independent(int g, int id) const
{
    int size = m_gv[g][1] >= 0 ? 2 : 1;
    if (m_wn[g] >= size)
    {
        return false;
    }
    if (m_wn[g] == 0)
    {
        return true;
    }
    for (int s = 0; s < m_wn[g]; s++)
    {
        if (m_ws[g][s] == id)
        {
            return false;
        }
    }
    double p0 = 0.0;
    double p1 = 0.0;
    double q0 = 0.0;
    double q1 = 0.0;
    double h = 0.0;
    row(m_ws[g][0], &p0, &p1, &h);
    row(id, &q0, &q1, &h);
    return fabs(p0 * q1 - p1 * q0) > 1e-8;
}

void Lsc::startSet(const double* x)
{
    for (int g = 0; g < m_n; g++)
    {
        m_wn[g] = 0;
        if (m_gv[g][0] < 0)
        {
            continue;
        }
        for (int s = 0; s < 2; s++)
        {
            int j = m_gv[g][s];
            if (j >= 0 && m_fix[j])
            {
                m_ws[g][m_wn[g]++] = 2 * j;
            }
        }
        for (int s = 0; s < 2; s++)
        {
            int j = m_gv[g][s];
            if (j < 0 || m_fix[j])
            {
                continue;
            }
            int id = x[j] <= m_lb[j] ? 2 * j : (x[j] >= m_ub[j] ? 2 * j + 1 : -1);
            if (id >= 0 && independent(g, id))
            {
                m_ws[g][m_wn[g]++] = id;
            }
        }
        if (m_gv[g][1] < 0)
        {
            continue;
        }
        double x0 = x[m_gv[g][0]];
        double x1 = x[m_gv[g][1]];
        for (int c = m_first[g]; c >= 0; c = m_next[c])
        {
            double slack = m_ch[c] - m_ca[c][0] * x0 - m_ca[c][1] * x1;
            int id = 2 * m_n + c;
            if (slack <= 1e-12 * (1.0 + fabs(m_ch[c])) && independent(g, id))
            {
                m_ws[g][m_wn[g]++] = id;
            }
        }
    }
}

// Отражение Хаусхолдера столбца col (строки k..m-1): R_kk в col[k],
// вектор отражения в u, возврат 2 / (u^T u) либо 0 (столбец нулевой).
//   col -- столбец; u -- вектор отражения; k -- номер шага; m -- строк
static double Householder(double* col, double* u, int k, int m)
{
    double s = 0.0;
    for (int i = k; i < m; i++)
    {
        s += col[i] * col[i];
    }
    s = sqrt(s);
    if (s <= 1e-300)
    {
        return 0.0;
    }
    double alpha = col[k] > 0.0 ? -s : s;
    double uu = 0.0;
    for (int i = k; i < m; i++)
    {
        u[i] = col[i];
    }
    u[k] -= alpha;
    for (int i = k; i < m; i++)
    {
        uu += u[i] * u[i];
    }
    col[k] = alpha;
    return 2.0 / uu;
}

// Применение отражения (I - beta u u^T) к вектору v (строки k..m-1).
//   v -- вектор; u -- вектор отражения; beta -- множитель; k, m -- строки
static void Reflect(double* v, const double* u, double beta, int k, int m)
{
    double t = 0.0;
    for (int i = k; i < m; i++)
    {
        t += u[i] * v[i];
    }
    t *= beta;
    for (int i = k; i < m; i++)
    {
        v[i] -= t * u[i];
    }
}

bool Lsc::subproblem(const double* A, const double* b, int m, const double* x)
{
    int n = m_n;
    for (int i = 0; i < m; i++)
    {
        m_r[i] = b[i];
    }
    int nf = 0;
    for (int g = 0; g < n; g++)
    {
        int j = m_gv[g][0];
        if (j < 0)
        {
            continue;
        }
        int k = m_gv[g][1];
        int size = k >= 0 ? 2 : 1;
        double base0 = 0.0;
        double base1 = 0.0;
        int wn = m_wn[g];
        if (wn == size)
        {
            double p0 = 0.0;
            double p1 = 0.0;
            double h = 0.0;
            row(m_ws[g][0], &p0, &p1, &h);
            if (size == 1)
            {
                base0 = h / p0;
            }
            else
            {
                double q0 = 0.0;
                double q1 = 0.0;
                double hq = 0.0;
                row(m_ws[g][1], &q0, &q1, &hq);
                double det = p0 * q1 - p1 * q0;
                double sp = h - p0 * x[j] - p1 * x[k];
                double sq = hq - q0 * x[j] - q1 * x[k];
                base0 = x[j] + (sp * q1 - p1 * sq) / det;
                base1 = x[k] + (p0 * sq - sp * q0) / det;
            }
        }
        else if (wn == 1)
        {
            double p0 = 0.0;
            double p1 = 0.0;
            double h = 0.0;
            row(m_ws[g][0], &p0, &p1, &h);
            double s = p0 * p0 + p1 * p1;
            double sp = (h - p0 * x[j] - p1 * x[k]) / s;
            base0 = x[j] + p0 * sp;
            base1 = x[k] + p1 * sp;
            double sn = sqrt(s);
            m_fg[nf] = g;
            m_ft[nf][0] = -p1 / sn;
            m_ft[nf][1] = p0 / sn;
            nf++;
        }
        else
        {
            m_fg[nf] = g;
            m_ft[nf][0] = 1.0;
            m_ft[nf][1] = 0.0;
            nf++;
            if (size == 2)
            {
                m_fg[nf] = g;
                m_ft[nf][0] = 0.0;
                m_ft[nf][1] = 1.0;
                nf++;
            }
        }
        m_xs[j] = base0;
        if (k >= 0)
        {
            m_xs[k] = base1;
        }
        for (int i = 0; i < m; i++)
        {
            double v = A[i * n + j] * base0;
            v += k >= 0 ? A[i * n + k] * base1 : 0.0;
            m_r[i] -= v;
        }
    }
    for (int f = 0; f < nf; f++)
    {
        int j = m_gv[m_fg[f]][0];
        int k = m_gv[m_fg[f]][1];
        for (int i = 0; i < m; i++)
        {
            double v = A[i * n + j] * m_ft[f][0];
            v += k >= 0 ? A[i * n + k] * m_ft[f][1] : 0.0;
            m_W[f * m + i] = v;
        }
    }
    for (int f = 0; f < nf; f++)
    {
        double beta = Householder(&m_W[f * m], m_u, f, m);
        if (beta == 0.0)
        {
            return false;
        }
        for (int q = f + 1; q < nf; q++)
        {
            Reflect(&m_W[q * m], m_u, beta, f, m);
        }
        Reflect(m_r, m_u, beta, f, m);
    }
    for (int f = nf - 1; f >= 0; f--)
    {
        double s = m_r[f];
        for (int q = f + 1; q < nf; q++)
        {
            s -= m_W[q * m + f] * m_z[q];
        }
        m_z[f] = s / m_W[f * m + f];
    }
    for (int f = 0; f < nf; f++)
    {
        int j = m_gv[m_fg[f]][0];
        int k = m_gv[m_fg[f]][1];
        m_xs[j] += m_ft[f][0] * m_z[f];
        if (k >= 0)
        {
            m_xs[k] += m_ft[f][1] * m_z[f];
        }
    }
    return true;
}

int Lsc::worstMultiplier(const double* A, const double* b, int m, const double* x,
                         double tol) const
{
    int n = m_n;
    double r[LSC_MAX_ROWS];
    for (int i = 0; i < m; i++)
    {
        double s = -b[i];
        for (int k = 0; k < n; k++)
        {
            s += A[i * n + k] * x[k];
        }
        r[i] = s;
    }
    int worst = -1;
    double wl = -tol;
    for (int g = 0; g < n; g++)
    {
        int j = m_gv[g][0];
        if (j < 0 || m_wn[g] == 0)
        {
            continue;
        }
        int k = m_gv[g][1];
        double g0 = 0.0;
        double g1 = 0.0;
        for (int i = 0; i < m; i++)
        {
            g0 += A[i * n + j] * r[i];
            g1 += k >= 0 ? A[i * n + k] * r[i] : 0.0;
        }
        double lam[2] = { 0.0, 0.0 };
        double p0 = 0.0;
        double p1 = 0.0;
        double h = 0.0;
        double cond = 1.0;
        row(m_ws[g][0], &p0, &p1, &h);
        if (m_wn[g] == 1)
        {
            lam[0] = -(g0 * p0 + g1 * p1) / (p0 * p0 + p1 * p1);
        }
        else
        {
            double q0 = 0.0;
            double q1 = 0.0;
            row(m_ws[g][1], &q0, &q1, &h);
            double det = p0 * q1 - q0 * p1;
            lam[0] = (-g0 * q1 + g1 * q0) / det;
            lam[1] = (-p0 * g1 + p1 * g0) / det;
            cond = 1.0 / fabs(det);
        }
        for (int s = 0; s < m_wn[g]; s++)
        {
            int id = m_ws[g][s];
            if (id < 2 * n && m_fix[id / 2])
            {
                continue;
            }
            if (lam[s] < -tol * cond && lam[s] / cond < wl)
            {
                wl = lam[s] / cond;
                worst = 2 * g + s;
            }
        }
    }
    return worst;
}

// Допуск на множители по масштабу задачи (как в lib_bvls.cpp): ошибка
// округления g = A^T r порядка eps |A| (|A| |x| + |b|), запас ~ 500 раз.
//   A, b -- задача; m, n -- размеры; lb, ub -- границы
static double MultiplierTol(const double* A, const double* b, int m, int n,
                            const double* lb, const double* ub)
{
    double na = 0.0;
    double nb = 0.0;
    double nx = 1.0;
    for (int i = 0; i < m; i++)
    {
        nb = fabs(b[i]) > nb ? fabs(b[i]) : nb;
        for (int j = 0; j < n; j++)
        {
            na = fabs(A[i * n + j]) > na ? fabs(A[i * n + j]) : na;
        }
    }
    for (int j = 0; j < n; j++)
    {
        nx = fabs(lb[j]) > nx ? fabs(lb[j]) : nx;
        nx = fabs(ub[j]) > nx ? fabs(ub[j]) : nx;
    }
    return 1e-13 * na * (na * nx + nb);
}

int Lsc::solve(const double* A, const double* b, int m, double* x, int maxIter)
{
    int n = m_n;
    if (n < 1 || m < n || m > LSC_MAX_ROWS)
    {
        return -2;
    }
    double xc[LSC_MAX_COLS];
    for (int j = 0; j < n; j++)
    {
        xc[j] = x[j] < m_lb[j] ? m_lb[j] : (x[j] > m_ub[j] ? m_ub[j] : x[j]);
        xc[j] = m_fix[j] ? m_lb[j] : xc[j];
    }
    for (int c = 0; c < m_nc; c++)
    {
        int g = m_cg[c];
        double v = m_ca[c][0] * xc[m_gv[g][0]] + m_ca[c][1] * xc[m_gv[g][1]] - m_ch[c];
        if (v > 1e-9 * (1.0 + fabs(m_ch[c])))
        {
            return -3;
        }
    }
    memcpy(x, xc, sizeof(double) * (size_t)n);
    startSet(x);
    double tol = MultiplierTol(A, b, m, n, m_lb, m_ub);
    int dropped = -1;
    for (int it = 0; it < maxIter; it++)
    {
        if (!subproblem(A, b, m, x))
        {
            return -2;
        }
        double alpha = 1.0;
        int blk = -1;
        int bg = -1;
        for (int g = 0; g < n; g++)
        {
            int j = m_gv[g][0];
            if (j < 0)
            {
                continue;
            }
            int k = m_gv[g][1];
            int size = k >= 0 ? 2 : 1;
            if (m_wn[g] >= size)
            {
                continue;
            }
            double p0 = m_xs[j] - x[j];
            double p1 = k >= 0 ? m_xs[k] - x[k] : 0.0;
            double pn = sqrt(p0 * p0 + p1 * p1);
            if (pn <= 0.0)
            {
                continue;
            }
            int ids[4];
            int ni = 0;
            for (int s = 0; s < size; s++)
            {
                int v = m_gv[g][s];
                if (!m_fix[v])
                {
                    ids[ni++] = 2 * v;
                    ids[ni++] = 2 * v + 1;
                }
            }
            int c = k >= 0 ? m_first[g] : -1;
            for (int q = 0; q < ni || c >= 0; q++)
            {
                int id = 0;
                if (q < ni)
                {
                    id = ids[q];
                }
                else
                {
                    id = 2 * n + c;
                    c = m_next[c];
                }
                bool in = false;
                for (int s = 0; s < m_wn[g]; s++)
                {
                    in = in || m_ws[g][s] == id;
                }
                if (in)
                {
                    continue;
                }
                double a0 = 0.0;
                double a1 = 0.0;
                double h = 0.0;
                row(id, &a0, &a1, &h);
                double s = a0 * p0 + a1 * p1;
                if (s <= 1e-14 * pn || !independent(g, id))
                {
                    continue;
                }
                double xg1 = k >= 0 ? x[k] : 0.0;
                double a = (h - a0 * x[j] - a1 * xg1) / s;
                a = a < 0.0 ? 0.0 : a;
                if (a < alpha)
                {
                    alpha = a;
                    blk = id;
                    bg = g;
                }
            }
        }
        if (blk < 0)
        {
            memcpy(x, m_xs, sizeof(double) * (size_t)n);
            int w = worstMultiplier(A, b, m, x, tol);
            if (w < 0)
            {
                return it + 1;
            }
            int g = w / 2;
            int s = w % 2;
            dropped = m_ws[g][s];
            m_ws[g][s] = m_ws[g][m_wn[g] - 1];
            m_wn[g]--;
            continue;
        }
        if (blk == dropped && alpha <= 1e-12)
        {
            // Только что исключённое ограничение сразу же останавливает шаг:
            // знак его множителя -- шум округления (почти параллельные
            // ограничения вершины); x -- решение прежнего рабочего набора.
            m_ws[bg][m_wn[bg]++] = blk;
            return it + 1;
        }
        dropped = -1;
        for (int j = 0; j < n; j++)
        {
            x[j] += alpha * (m_xs[j] - x[j]);
        }
        if (blk < 2 * n)
        {
            int j = blk / 2;
            x[j] = (blk & 1) == 0 ? m_lb[j] : m_ub[j];
        }
        if (independent(bg, blk))
        {
            m_ws[bg][m_wn[bg]++] = blk;
        }
    }
    return -1;
}

} // namespace gnc

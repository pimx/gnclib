//============================================================================
// lib_bvls.cpp -- наименьшие квадраты с границами переменных (1.7.0), см.
// lib_bvls.h.
//============================================================================
#include "lib_bvls.h"

#include <math.h>
#include <string.h>

namespace gnc
{

Bvls::Bvls()
{
    m_nFree = 0;
    m_n = 0;
    memset(m_state, 0, sizeof(m_state));
}

int Bvls::boundCount() const
{
    int c = 0;
    for (int j = 0; j < m_n; j++)
    {
        if (m_state[j] == BVLS_LOWER || m_state[j] == BVLS_UPPER)
        {
            c++;
        }
    }
    return c;
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

bool Bvls::subproblem(const double* A, const double* b, int m, int n, const double* x)
{
    for (int i = 0; i < m; i++)
    {
        double s = b[i];
        for (int j = 0; j < n; j++)
        {
            if (m_state[j] != BVLS_FREE)
            {
                s -= A[i * n + j] * x[j];
            }
        }
        m_r[i] = s;
    }
    int nf = m_nFree;
    for (int k = 0; k < nf; k++)
    {
        for (int i = 0; i < m; i++)
        {
            m_W[k * m + i] = A[i * n + m_free[k]];
        }
    }
    for (int k = 0; k < nf; k++)
    {
        double beta = Householder(&m_W[k * m], m_u, k, m);
        if (beta == 0.0)
        {
            return false;
        }
        for (int j = k + 1; j < nf; j++)
        {
            Reflect(&m_W[j * m], m_u, beta, k, m);
        }
        Reflect(m_r, m_u, beta, k, m);
    }
    for (int k = nf - 1; k >= 0; k--)
    {
        double s = m_r[k];
        for (int j = k + 1; j < nf; j++)
        {
            s -= m_W[j * m + k] * m_z[j];
        }
        m_z[k] = s / m_W[k * m + k];
    }
    return true;
}

int Bvls::worstMultiplier(const double* A, const double* b, int m, int n,
                          const double* x, double tol) const
{
    double r[BVLS_MAX_ROWS];
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
    for (int j = 0; j < n; j++)
    {
        if (m_state[j] != BVLS_LOWER && m_state[j] != BVLS_UPPER)
        {
            continue;
        }
        double g = 0.0;
        for (int i = 0; i < m; i++)
        {
            g += A[i * n + j] * r[i];
        }
        double lambda = m_state[j] == BVLS_LOWER ? g : -g;
        if (lambda < wl)
        {
            wl = lambda;
            worst = j;
        }
    }
    return worst;
}

// Начальное состояние: границы, закреплённые переменные, тёплый старт.
//   lb, ub, x -- границы и приближение; n -- число; state -- состояния
static void StartSet(const double* lb, const double* ub, double* x, int n, int* state)
{
    for (int j = 0; j < n; j++)
    {
        if (ub[j] - lb[j] <= 1e-15 * (1.0 + fabs(lb[j])))
        {
            x[j] = lb[j];
            state[j] = BVLS_FIXED;
        }
        else if (x[j] <= lb[j])
        {
            x[j] = lb[j];
            state[j] = BVLS_LOWER;
        }
        else if (x[j] >= ub[j])
        {
            x[j] = ub[j];
            state[j] = BVLS_UPPER;
        }
        else
        {
            state[j] = BVLS_FREE;
        }
    }
}

// Допуск на множители по масштабу задачи: ошибка округления g = A^T r
// порядка eps |A| (|A| |x| + |b|); допуск -- с запасом ~ 500 раз.
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

int Bvls::solve(const double* A, const double* b, int m, int n,
                const double* lb, const double* ub, double* x, int maxIter)
{
    if (n < 1 || n > BVLS_MAX_COLS || m < n || m > BVLS_MAX_ROWS)
    {
        return -2;
    }
    m_n = n;
    StartSet(lb, ub, x, n, m_state);
    double tol = MultiplierTol(A, b, m, n, lb, ub);
    for (int it = 0; it < maxIter; it++)
    {
        m_nFree = 0;
        for (int j = 0; j < n; j++)
        {
            if (m_state[j] == BVLS_FREE)
            {
                m_free[m_nFree++] = j;
            }
        }
        if (m_nFree > 0 && !subproblem(A, b, m, n, x))
        {
            return -2;
        }
        double alpha = 1.0;
        int blk = -1;
        int side = 0;
        for (int k = 0; k < m_nFree; k++)
        {
            int j = m_free[k];
            double z = m_z[k];
            if (z < lb[j])
            {
                double a = (lb[j] - x[j]) / (z - x[j]);
                if (a < alpha)
                {
                    alpha = a;
                    blk = j;
                    side = BVLS_LOWER;
                }
            }
            else if (z > ub[j])
            {
                double a = (ub[j] - x[j]) / (z - x[j]);
                if (a < alpha)
                {
                    alpha = a;
                    blk = j;
                    side = BVLS_UPPER;
                }
            }
        }
        if (blk < 0)
        {
            for (int k = 0; k < m_nFree; k++)
            {
                x[m_free[k]] = m_z[k];
            }
            int w = worstMultiplier(A, b, m, n, x, tol);
            if (w < 0)
            {
                return it + 1;
            }
            m_state[w] = BVLS_FREE;
            continue;
        }
        for (int k = 0; k < m_nFree; k++)
        {
            int j = m_free[k];
            double v = x[j] + alpha * (m_z[k] - x[j]);
            x[j] = v < lb[j] ? lb[j] : (v > ub[j] ? ub[j] : v);
        }
        x[blk] = side == BVLS_LOWER ? lb[blk] : ub[blk];
        m_state[blk] = side;
    }
    return -1;
}

} // namespace gnc

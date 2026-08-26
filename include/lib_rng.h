/* lib_rng.h -- генератор случайных чисел (xorshift64*), без зависимостей.
 *
 * Детерминирован при фиксированном зерне -- важно для повторяемости
 * прогонов. Равномерное [0,1) и нормальное N(0,1) (Бокс-Мюллер).
 * Без STL, без исключений, без динамических аллокаций.
 */
#ifndef LIB_RNG_H
#define LIB_RNG_H

#include <math.h>

class Rng
{
public:
    explicit Rng(unsigned long long seed = 88172645463325252ULL)
        : m_state(seed ? seed : 1ULL), m_spare(0.0), m_hasSpare(false)
    {
    }

    void reseed(unsigned long long seed)
    {
        m_state = seed ? seed : 1ULL;
        m_hasSpare = false;
    }

    /* Равномерное [0,1). */
    double uniform()
    {
        m_state ^= m_state << 13;
        m_state ^= m_state >> 7;
        m_state ^= m_state << 17;
        unsigned long long v = m_state * 2685821657736338717ULL;
        return (double)(v >> 11) * (1.0 / 9007199254740992.0);
    }

    /* Равномерное [a,b). */
    double uniformRange(double a, double b)
    {
        return a + (b - a) * uniform();
    }

    /* Нормальное N(0,1), метод Бокса-Мюллера (с кэшем второго значения). */
    double gauss()
    {
        if (m_hasSpare)
        {
            m_hasSpare = false;
            return m_spare;
        }
        double u1 = uniform();
        double u2 = uniform();
        if (u1 < 1e-16)
        {
            u1 = 1e-16;
        }
        double r = sqrt(-2.0 * log(u1));
        double a = 6.283185307179586 * u2;
        m_spare = r * sin(a);
        m_hasSpare = true;
        return r * cos(a);
    }

private:
    unsigned long long m_state;
    double m_spare;
    bool m_hasSpare;
};

#endif /* LIB_RNG_H */

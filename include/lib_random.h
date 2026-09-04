//======================================================================
//  lib_random.h
//  RNG - pseudo-random source (xorshift64*), no dependencies.
//  Verbatim port of lib_rng.h (class Rng) into the library conventions:
//  the generated sequences are IDENTICAL to lib_rng for the same seed,
//  so earlier Monte-Carlo runs remain reproducible.
//  Deterministic for a fixed seed. Uniform [0, 1) and normal N(0, 1)
//  (Box-Muller with a cached spare). C++17, no STL, no heap.
//  Constexpr-constructible: usable as a constant-initialized static.
//======================================================================

#ifndef LIB_RANDOM_H
#define LIB_RANDOM_H

#include <math.h>

#include "lib_frames.h"

namespace gnc
{

class RNG
{
public:

    constexpr explicit RNG(unsigned long long seed = 88172645463325252ULL)
        : m_state(1ULL)
        , m_spare(0.0)
        , m_has(false)
    {
        if (seed != 0ULL)
        {
            m_state = seed;
        }
    }

    constexpr void Reseed(unsigned long long seed)
    {
        m_state = 1ULL;
        if (seed != 0ULL)
        {
            m_state = seed;
        }
        m_has = false;
    }

    //  Uniform deviate in [0, 1) with 53 significant bits.
    double Uniform()
    {
        m_state ^= m_state << 13;
        m_state ^= m_state >> 7;
        m_state ^= m_state << 17;
        unsigned long long v = m_state * 2685821657736338717ULL;
        return (double)(v >> 11) * (1.0 / 9007199254740992.0);
    }

    //  Uniform deviate in [a, b).
    double UniformRange(double a, double b)
    {
        return a + (b - a) * Uniform();
    }

    //  Standard normal deviate (Box-Muller, pair cached).
    double Gauss()
    {
        if (m_has)
        {
            m_has = false;
            return m_spare;
        }
        double u1 = Uniform();
        double u2 = Uniform();
        if (u1 < 1.0e-16)
        {
            u1 = 1.0e-16;
        }
        double r = sqrt(-2.0 * log(u1));
        double a = 2.0 * PI * u2;
        m_spare = r * sin(a);
        m_has = true;
        return r * cos(a);
    }

private:

    unsigned long long m_state;
    double m_spare;
    bool m_has;
};

}   // namespace gnc

#endif  // LIB_RANDOM_H

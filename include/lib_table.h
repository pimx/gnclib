//============================================================================
// lib_table.h -- интерполяция по таблицам с удержанием краёв (1.6.0).
//
// TableLocate -- строка и доля значения в возрастающей таблице;
// TableBilinear -- билинейная интерполяция по двумерной таблице (строки --
// первый аргумент, столбцы -- второй), хранимой построчно.
// C++17, без STL, без кучи, без исключений.
//============================================================================
#ifndef GNCLIB_LIB_TABLE_H
#define GNCLIB_LIB_TABLE_H

namespace gnc
{

// Номер строки и доля значения в возрастающей таблице, края удерживаются:
// ниже первой -- строка 0, доля 0; выше последней -- строка n - 2, доля 1.
//   tab -- возрастающие значения; n -- их число (>= 2)
//   x -- значение; f -- выход: доля к следующей строке
// Возвращает номер строки 0 .. n - 2.
inline int TableLocate(const double* tab, int n, double x, double* f)
{
    if (x <= tab[0])
    {
        *f = 0.0;
        return 0;
    }
    if (x >= tab[n - 1])
    {
        *f = 1.0;
        return n - 2;
    }
    int i = 0;
    while (i < n - 2 && x > tab[i + 1])
    {
        i = i + 1;
    }
    *f = (x - tab[i]) / (tab[i + 1] - tab[i]);
    return i;
}

// Билинейная интерполяция: сначала вдоль строки (столбцы iC, iC + 1),
// затем между строками iR, iR + 1.
//   tab -- таблица построчно, cols -- число столбцов
//   iR, fR -- строка и доля (TableLocate по первому аргументу)
//   iC, fC -- столбец и доля (TableLocate по второму аргументу)
inline double TableBilinear(const double* tab, int cols, int iR, double fR, int iC, double fC)
{
    const double* r0 = tab + iR * cols;
    const double* r1 = r0 + cols;
    double a = r0[iC] + fC * (r0[iC + 1] - r0[iC]);
    double b = r1[iC] + fC * (r1[iC + 1] - r1[iC]);
    return a + fR * (b - a);
}

} // namespace gnc

#endif // GNCLIB_LIB_TABLE_H

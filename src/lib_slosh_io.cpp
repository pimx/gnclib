//============================================================================
// lib_slosh_io.cpp -- таблицы параметров плескания в CSV (1.7.0): строка
// заголовка, 19 столбцов в порядке SloshRow через запятую:
//   q,zSurf,zLiq,m0r,z0,i0r,iRollr,mr1,mr2,mr3,z1,z2,z3,w2a1,w2a2,w2a3,
//   eta1,eta2,eta3
// Чтение принимает и пробелы. См. lib_slosh.h.
//============================================================================
#include "lib_slosh.h"

#include <stdio.h>
#include <stdlib.h>

namespace gnc
{

namespace
{

const int NCOL = 7 + 4 * SLOSH_NMODE;

// Закрыть файл и записать причину отказа.
bool Refuse(FILE* f, char* err, int errSize, const char* text)
{
    if (f != 0)
    {
        fclose(f);
    }
    if (err != 0 && errSize > 0)
    {
        snprintf(err, (size_t)errSize, "%s", text);
    }
    return false;
}

// Значения строки в порядке столбцов.
void RowToCols(double q, const SloshRow& r, double* c)
{
    c[0] = q;
    c[1] = r.zSurf;
    c[2] = r.zLiq;
    c[3] = r.m0r;
    c[4] = r.z0;
    c[5] = r.i0r;
    c[6] = r.iRollr;
    for (int n = 0; n < SLOSH_NMODE; ++n)
    {
        c[7 + n] = r.mr[n];
        c[7 + SLOSH_NMODE + n] = r.z[n];
        c[7 + 2 * SLOSH_NMODE + n] = r.w2a[n];
        c[7 + 3 * SLOSH_NMODE + n] = r.eta[n];
    }
}

// Строка из значений в порядке столбцов (c[0] -- q -- не входит в строку).
void ColsToRow(const double* c, SloshRow* r)
{
    r->zSurf = c[1];
    r->zLiq = c[2];
    r->m0r = c[3];
    r->z0 = c[4];
    r->i0r = c[5];
    r->iRollr = c[6];
    for (int n = 0; n < SLOSH_NMODE; ++n)
    {
        r->mr[n] = c[7 + n];
        r->z[n] = c[7 + SLOSH_NMODE + n];
        r->w2a[n] = c[7 + 2 * SLOSH_NMODE + n];
        r->eta[n] = c[7 + 3 * SLOSH_NMODE + n];
    }
}

// Пропустить пробелы и табуляции; true -- строка пустая или комментарий.
bool SkipLine(const char** s)
{
    while (**s == ' ' || **s == '\t')
    {
        ++*s;
    }
    char c = **s;
    return c == 0 || c == '\n' || c == '\r' || c == ';' || c == '#';
}

// Числа строки через запятые или пробелы, не больше NCOL + 1; возвращает
// их число (NCOL + 1 -- лишние числа).
//   s -- строка; c -- выход
int ParseLine(const char* s, double* c)
{
    int n = 0;
    const char* p = s;
    while (n <= NCOL)
    {
        while (*p == ' ' || *p == '\t' || *p == ',')
        {
            ++p;
        }
        char* end = 0;
        double v = strtod(p, &end);
        if (end == p)
        {
            break;
        }
        if (n < NCOL)
        {
            c[n] = v;
        }
        ++n;
        p = end;
    }
    return n;
}

} // namespace

//----------------------------------------------------------------------------
// Запись таблицы (см. lib_slosh.h), 10 значащих цифр.
//   path -- файл; tab -- таблица; err, errSize -- причина отказа
//----------------------------------------------------------------------------
bool SloshTableSaveCsv(const char* path, const SloshTable& tab, char* err, int errSize)
{
    char text[512];
    FILE* f = (path != 0) ? fopen(path, "w") : 0;
    if (f == 0)
    {
        snprintf(text, sizeof(text), "the slosh table %s cannot be written",
                 (path != 0) ? path : "(null)");
        return Refuse(0, err, errSize, text);
    }
    fprintf(f, "q,zSurf,zLiq,m0r,z0,i0r,iRollr,mr1,mr2,mr3,z1,z2,z3,"
               "w2a1,w2a2,w2a3,eta1,eta2,eta3\n");
    for (int i = 0; i < tab.n; ++i)
    {
        double c[NCOL];
        RowToCols(tab.q[i], tab.row[i], c);
        for (int k = 0; k < NCOL; ++k)
        {
            fprintf(f, (k + 1 < NCOL) ? "%.10g," : "%.10g\n", c[k]);
        }
    }
    bool ok = (ferror(f) == 0);
    fclose(f);
    if (!ok)
    {
        snprintf(text, sizeof(text), "an error while writing %s", path);
        return Refuse(0, err, errSize, text);
    }
    return true;
}

//----------------------------------------------------------------------------
// Чтение таблицы (см. lib_slosh.h). Заголовок -- строки до первой строки
// чисел; после неё строка не из чисел -- отказ.
//   path -- файл; out -- выход; err, errSize -- причина отказа
//----------------------------------------------------------------------------
bool SloshTableLoadCsv(const char* path, SloshTable* out, char* err, int errSize)
{
    if (out == 0)
    {
        return Refuse(0, err, errSize, "no output table");
    }
    char text[512];
    FILE* f = (path != 0) ? fopen(path, "r") : 0;
    if (f == 0)
    {
        snprintf(text, sizeof(text), "the slosh table %s cannot be opened",
                 (path != 0) ? path : "(null)");
        return Refuse(0, err, errSize, text);
    }
    SloshTable t;
    char line[1024];
    int lineNo = 0;
    while (fgets(line, (int)sizeof(line), f) != 0)
    {
        lineNo = lineNo + 1;
        const char* s = line;
        if (SkipLine(&s))
        {
            continue;
        }
        double c[NCOL];
        int n = ParseLine(s, c);
        if (n == 0 && t.n == 0)
        {
            continue; // заголовок
        }
        if (n != NCOL)
        {
            snprintf(text, sizeof(text), "line %d of %s is not %d numbers", lineNo, path, NCOL);
            return Refuse(f, err, errSize, text);
        }
        if (t.n >= SLOSH_NROW)
        {
            snprintf(text, sizeof(text), "%s holds more than %d rows", path, SLOSH_NROW);
            return Refuse(f, err, errSize, text);
        }
        if (t.n > 0 && c[0] <= t.q[t.n - 1])
        {
            snprintf(text, sizeof(text), "q of %s does not rise (line %d)", path, lineNo);
            return Refuse(f, err, errSize, text);
        }
        t.q[t.n] = c[0];
        ColsToRow(c, &t.row[t.n]);
        t.n = t.n + 1;
    }
    fclose(f);
    if (t.n < 1)
    {
        snprintf(text, sizeof(text), "%s holds no rows", path);
        return Refuse(0, err, errSize, text);
    }
    *out = t;
    return true;
}

} // namespace gnc

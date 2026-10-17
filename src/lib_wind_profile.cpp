//============================================================================
// lib_wind_profile.cpp -- чтение профиля ветра и ветер на высоте (1.6.0).
//============================================================================
#include "lib_wind_profile.h"

#include <stdio.h>

#include "lib_util.h"

namespace gnc
{

namespace
{

// Закрыть файл и записать причину отказа.
bool Refuse(FILE* f, char* err, int errSize, const char* text)
{
    if (f != 0)
    {
        fclose(f);
    }
    snprintf(err, (size_t)errSize, "%s", text);
    return false;
}

// Пропустить пробелы и табуляции; true -- строка пустая или комментарий.
bool SkipLine(char** s)
{
    while (**s == ' ' || **s == '\t')
    {
        ++*s;
    }
    char c = **s;
    return c == 0 || c == '\n' || c == '\r' || c == ';' || c == '#';
}

// Запятые -- в пробелы: файл с запятыми и файл с пробелами читаются одинаково.
void CommasToSpaces(char* s)
{
    for (char* c = s; *c != 0; ++c)
    {
        if (*c == ',')
        {
            *c = ' ';
        }
    }
}

} // namespace

bool ReadWindProfile(const char* path, WindProfile* out, char* err, int errSize)
{
    if (path == 0 || out == 0 || err == 0 || errSize < 2)
    {
        return false;
    }
    err[0] = 0;
    char text[512];
    FILE* f = fopen(path, "r");
    if (f == 0)
    {
        snprintf(text, sizeof(text), "the wind profile %s cannot be opened", path);
        return Refuse(0, err, errSize, text);
    }
    WindProfile p;
    char line[256];
    int lineNo = 0;
    while (fgets(line, (int)sizeof(line), f) != 0)
    {
        lineNo = lineNo + 1;
        char* s = line;
        if (SkipLine(&s))
        {
            continue;
        }
        double alt = 0.0;
        double speed = 0.0;
        double dir = 0.0;
        CommasToSpaces(s);
        if (sscanf(s, "%lf %lf %lf", &alt, &speed, &dir) != 3)
        {
            if (p.count == 0)
            {
                continue; // заголовок
            }
            snprintf(text, sizeof(text), "line %d of %s is not three numbers", lineNo, path);
            return Refuse(f, err, errSize, text);
        }
        if (p.count >= WIND_PROFILE_MAX)
        {
            snprintf(text, sizeof(text), "%s holds more than %d rows", path, WIND_PROFILE_MAX);
            return Refuse(f, err, errSize, text);
        }
        if (p.count > 0 && alt <= p.row[p.count - 1].alt)
        {
            snprintf(text, sizeof(text), "the altitudes of %s do not rise (line %d)", path, lineNo);
            return Refuse(f, err, errSize, text);
        }
        if (speed < 0.0)
        {
            snprintf(text, sizeof(text), "a negative speed in %s (line %d)", path, lineNo);
            return Refuse(f, err, errSize, text);
        }
        p.row[p.count].alt = alt;
        p.row[p.count].speed = speed;
        p.row[p.count].dir = Deg(dir);
        p.count = p.count + 1;
    }
    fclose(f);
    if (p.count < 2)
    {
        snprintf(text, sizeof(text), "%s holds fewer than two rows", path);
        return Refuse(0, err, errSize, text);
    }
    *out = p;
    return true;
}

void WindProfileAt(const WindProfile& profile, double alt, double* speed, double* from)
{
    if (speed == 0 || from == 0 || profile.count < 1)
    {
        return;
    }
    if (alt <= profile.row[0].alt)
    {
        *speed = profile.row[0].speed;
        *from = profile.row[0].dir;
        return;
    }
    int last = profile.count - 1;
    if (alt >= profile.row[last].alt)
    {
        *speed = profile.row[last].speed;
        *from = profile.row[last].dir;
        return;
    }
    int i = 0;
    while (i < last - 1 && alt > profile.row[i + 1].alt)
    {
        i = i + 1;
    }
    double span = profile.row[i + 1].alt - profile.row[i].alt;
    double f = span > 1.0e-12 ? (alt - profile.row[i].alt) / span : 0.0;
    *speed = profile.row[i].speed + f * (profile.row[i + 1].speed - profile.row[i].speed);
    *from = profile.row[i].dir + f * (profile.row[i + 1].dir - profile.row[i].dir);
}

} // namespace gnc

//============================================================================
// lib_wind_profile.h -- профиль ветра по высоте из файла (1.6.0).
//
// ФОРМАТ ФАЙЛА. Текст, строка -- три числа через запятые или пробелы:
//   высота, м   скорость, м/с   азимут, ОТКУДА дует, град
// Пустые строки и строки, начинающиеся с ';' или '#', пропускаются; так же
// пропускается заголовок до первой строки чисел. Высоты строго растут,
// строк от двух до WIND_PROFILE_MAX, скорость неотрицательна.
// Между строками -- линейно (скорость и азимут), за краями -- крайняя
// строка.
// C++17, без STL, без кучи, без исключений.
//============================================================================
#ifndef GNCLIB_LIB_WIND_PROFILE_H
#define GNCLIB_LIB_WIND_PROFILE_H

namespace gnc
{

inline constexpr int WIND_PROFILE_MAX = 64; // строк в профиле

struct WindProfileRow
{
    double alt = 0.0;   // высота, м
    double speed = 0.0; // скорость, м/с
    double dir = 0.0;   // азимут, откуда дует, рад
};

struct WindProfile
{
    int count = 0;
    WindProfileRow row[WIND_PROFILE_MAX] = { };
};

// Чтение профиля из файла.
//   path -- файл; out -- выход: профиль (меняется только при успехе)
//   err -- выход: причина отказа (ASCII); errSize -- ёмкость err
// Возвращает false, если файл не читается, строк меньше двух или больше
// WIND_PROFILE_MAX, высоты не растут, строка не из трёх чисел или
// скорость отрицательна.
bool ReadWindProfile(const char* path, WindProfile* out, char* err, int errSize);

// Ветер профиля на высоте: линейно между строками, за краями -- крайняя.
//   profile -- профиль; alt -- высота, м
//   speed -- выход: м/с; from -- выход: азимут, откуда дует, рад
void WindProfileAt(const WindProfile& profile, double alt, double* speed, double* from);

} // namespace gnc

#endif // GNCLIB_LIB_WIND_PROFILE_H

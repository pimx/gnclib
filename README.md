# gnclib

Статическая библиотека GNC: атмосфера, аэродинамика (полосовая модель,
Cx(M), экранирование ретро-струёй, экранный эффект), геодезия WGS-84
(ECEF/ENU, строгая сила тяжести J2 + центробежная), фреймы и
кватернионы, ветер, CLI, утилиты. C++17, без STL/исключений/RTTI/куч.

Сборка: `make` (g++, build/libgnc.a) или `build.bat` (VS2026,
build\gnc.lib). Заголовки -- include/, версия -- include/gnclib_version.h;
тег репозитория обязан совпадать с GNCLIB_VERSION.

Потребители подключают: -I<gnclib>/include и линкуют build/libgnc.a
(gnc.lib). Сборка потребителей против main запрещена -- только тег.

## 1.1.0 — интегрирование с силами на стадиях, вращающаяся Земля

- `lib_integrate.h` — `IntegrateRK4(s, aux, dt, fn, ctx)`: RK4 движения
  твёрдого тела с пересчётом сил на каждой стадии и вспомогательными
  состояниями (`AuxState`, до 32). `Integrate` (силы постоянны на шаге)
  остаётся лёгким вариантом.
- `lib_earth.h` — `EarthGravityLoc`, `EarthRateLoc`, `CoriolisAccLoc`,
  `SpecificForceLoc`: вращающаяся Земля в системе площадки (= ECEF).
- Проверка: `make test`.
- `lib_aero_strip.h` — `CylinderStripCT` и поправки Cd по числу Маха,
  ретро-плюм по коэффициенту тяги C_T.
- `lib_eskf.h` (+ `src/lib_eskf_*.cpp`) — общий навигационный ESKF
  `gnc::Eskf`; `lib_nav_types.h` — данные измерений; `lib_imu.h` —
  `ImuIdeal`, идеальные показания ИНС на вращающейся Земле с рычагом.
  Подключаются явно (не через `lib.h`).


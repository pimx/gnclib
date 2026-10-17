/*============================================================================
 * lib.h -- сводное включение основных заголовков библиотеки.
 *
 * Ядро (namespace gnc, именованные компоненты): lib_frames,
 * lib_rotation, lib_geodesy, lib_gravity, lib_atmosphere,
 * lib_aerodynamics, lib_wind, lib_random, lib_state, lib_dynamics,
 * lib_chain, lib_udp, lib_cli.
 *
 * Модули глобального пространства имён: lib_func, lib_task;
 * lib_linalg (Vec3/Mat3/Quat) НЕ включается сводным заголовком:
 * единственный потребитель безымянной линейной алгебры -- матричная
 * форма ESKF (lib_eskf_math.h). Типизация систем координат -- ТОЛЬКО
 * семейства gnc.
 *
 * Подключаются явно (не через lib.h): lib_eskf, lib_imu, lib_nav_types,
 * lib_util, lib_table, lib_eso, lib_notch, lib_gust, lib_wind_profile,
 * lib_orbit.
 *============================================================================*/
#ifndef LIB_H
#define LIB_H

#include "gnclib_version.h"

#include "lib_cli.h"			// разбор командной строки

#include "lib_aerodynamics.h"
#include "lib_atmosphere.h"
#include "lib_chain.h"
#include "lib_dynamics.h"
#include "lib_earth.h"			// вращающаяся Земля (1.1.0)
#include "lib_eskf.h"
#include "lib_eso.h"
#include "lib_integrate.h"		// RK4 с силами на стадиях (1.1.0)
#include "lib_frame_ops.h"		// операции над векторами площадки/BODY
#include "lib_frames.h"
#include "lib_func.h"
#include "lib_geodesy.h"
#include "lib_gravity.h"
#include "lib_random.h"
#include "lib_rotation.h"
#include "lib_state.h"
#include "lib_task.h"
#include "lib_udp.h"
#include "lib_wind.h"

// #include "lib_linalg.h"

#endif /* LIB_H */

/*============================================================================
 * lib.h -- сводное включение всех заголовков библиотеки (V60: cflib +
 * сохранённые модули libsim).
 *
 * Ядро cflib (namespace gnc, именованные компоненты): lib_frames,
 * lib_rotation, lib_geodesy, lib_gravity, lib_atmosphere,
 * lib_aerodynamics, lib_wind, lib_random, lib_state, lib_dynamics,
 * lib_chain, lib_udp, lib_cli.
 *
 * Сохранённые модули libsim (глобальное пространство имён): lib_func,
 * lib_task, lib_telemetry; lib_linalg (Vec3/Mat3/Quat) НЕ включается
 * сводным заголовком (V62): единственный потребитель безымянной
 * линейной алгебры -- матричная форма ESKF (src/eskf_math.h ->
 * src/navigation.h). Типизация систем координат -- ТОЛЬКО семейства
 * gnc (V61: прежние VecEnu/VecXBody/... удалены).
 *============================================================================*/
#ifndef LIB_H
#define LIB_H

#include "lib_frames.h"
#include "lib_rotation.h"
#include "lib_frame_ops.h" // операции над векторами площадки/BODY (V80)
#include "lib_geodesy.h"
#include "lib_gravity.h"
#include "lib_atmosphere.h"
#include "lib_aerodynamics.h"
#include "lib_wind.h"
#include "lib_random.h"
#include "lib_state.h"
#include "lib_dynamics.h"
#include "lib_chain.h"
#include "lib_udp.h"
#include "lib_cli.h"

#include "lib_func.h"
#include "lib_task.h"
#include "lib_telemetry.h"

#endif /* LIB_H */

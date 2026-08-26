/* lib_earth.h -- параметры Земли: WGS-84, геодезия, гравитация, ENU-фрейм.
 *
 * Состав:
 *   - константы эллипсоида WGS-84 (Wgs84);
 *   - статические преобразования: геодезические <-> ECEF, DCM ECEF <-> ENU,
 *     радиусы кривизны меридиана и первого вертикала (EarthGeodesy);
 *   - гравитация: центральное поле + J2, центробежное ускорение, полная
 *     сила тяжести вращающейся Земли, нормальная гравитация Сомильяны
 *     с высотной поправкой (EarthGravity);
 *   - class GeoFrame -- локальная ENU-площадка (x=East, y=North, z=Up),
 *     привязанная к геодезической точке; хранит константный вектор силы
 *     тяжести в ENU (гравитация J2 + центробежная -- статическая площадка
 *     на вращающейся Земле).
 *
 * Контур управления работает в ENU; ECEF нужен для привязки к ИНС/СНС.
 * Без STL, без исключений, без динамических аллокаций.
 */
#ifndef LIB_EARTH_H
#define LIB_EARTH_H

#include "lib_frames.h"

/* Геодезические координаты: широта [рад], долгота [рад], высота [м]. */
struct Geodetic
{
    double lat, lon, alt;
    Geodetic() : lat(0.0), lon(0.0), alt(0.0)
    {
    }
    Geodetic(double la, double lo, double al) : lat(la), lon(lo), alt(al)
    {
    }
};

/* Константы эллипсоида WGS-84. */
struct Wgs84
{
    static double a()
    {
        return 6378137.0;
    } /* большая полуось, м */
    static double f()
    {
        return 1.0 / 298.257223563;
    } /* сжатие */
    static double b()
    {
        return a() * (1.0 - f());
    } /* малая полуось, м */
    static double e2()
    {
        return f() * (2.0 - f());
    } /* квадрат эксцентриситета */
    static double mu()
    {
        return 3.986004418e14;
    } /* GM, м^3/с^2 */
    static double j2()
    {
        return 1.08262668355e-3;
    } /* вторая зональная гармоника */
    static double omega()
    {
        return 7.2921159e-5;
    } /* вращение Земли, рад/с */
};

/* Статические геодезические преобразования (не зависят от площадки). */
class EarthGeodesy
{
public:
    /* Геодезические -> ECEF через радиус кривизны первого вертикала. */
    static VecEcef geodeticToEcef(Geodetic g);

    /* ECEF -> геодезические (итерации Боуринга по широте, 6 проходов). */
    static Geodetic ecefToGeodetic(VecEcef r);

    /* Матрица R такая, что v_enu = R * v_ecef в точке (lat, lon). */
    static Mat3 ecefToEnuDcm(double lat, double lon);

    /* Радиус кривизны первого вертикала N(lat), м. */
    static double primeVerticalRadius(double lat);

    /* Радиус кривизны меридиана M(lat), м. */
    static double meridianRadius(double lat);

    /* Геоцентрическая широта по геодезической (на поверхности эллипсоида). */
    static double geocentricLat(double geodeticLat);
};

/* Модели гравитации. */
class EarthGravity
{
public:
    /* Гравитационное ускорение в ECEF: центральное поле + J2
     * (БЕЗ центробежного члена). */
    static VecEcef gravitationEcefJ2(VecEcef rEcef);

    /* Центробежное ускорение точки, неподвижной в ECEF: w x (w x r),
     * взятое с обратным знаком (кажущаяся сила во вращающейся системе). */
    static VecEcef centrifugalEcef(VecEcef rEcef);

    /* Полная сила тяжести на вращающейся Земле (гравитация J2 +
     * центробежная) -- то, что измеряет неподвижный акселерометр. */
    static VecEcef gravityEcef(VecEcef rEcef);

    /* Нормальная гравитация Сомильяны на эллипсоиде + высотная поправка
     * (свободный воздух 2-го порядка), м/с^2. Скаляр вдоль нормали. */
    static double normal(double lat, double h);
};

/* Локальная площадка ENU, привязанная к геодезической точке WGS-84. */
class GeoFrame
{
public:
    GeoFrame(); /* заглушка; перед работой вызвать reset() */
    explicit GeoFrame(Geodetic site); /* привязка при построении */

    /* Перепривязка фрейма к новой площадке: начало, повороты и константный
     * вектор силы тяжести в ENU (J2 + центробежная, повёрнут в ENU). */
    void reset(Geodetic site);

    /* Преобразования ПОЗИЦИЙ (перенос + поворот). */
    VecEnu ecefPosToEnu(VecEcef rEcef) const;
    VecEcef enuPosToEcef(VecEnu rEnu) const;

    /* Преобразования свободных ВЕКТОРОВ (только поворот):
     * скорости, силы, направления. */
    VecEnu ecefVecToEnu(VecEcef vEcef) const;
    VecEcef enuVecToEcef(VecEnu vEnu) const;

    /* Те же преобразования в конвенции ГОСТ 20058-80 (ГОСК): системы
     * ECEF и ГОСК эквивалентны, переход -- перестановка осей. */
    VecEnu goskPosToEnu(VecGosk rGosk) const
    {
        return ecefPosToEnu(rGosk.ecef());
    }
    VecGosk enuPosToGosk(VecEnu rEnu) const
    {
        return enuPosToEcef(rEnu).gosk();
    }
    VecEnu goskVecToEnu(VecGosk vGosk) const
    {
        return ecefVecToEnu(vGosk.ecef());
    }
    VecGosk enuVecToGosk(VecEnu vEnu) const
    {
        return enuVecToEcef(vEnu).gosk();
    }

    /* Константный вектор силы тяжести площадки в ENU (~(0,0,-9.81)). */
    VecEnu gravityEnu() const
    {
        return m_gEnu;
    }

    /* Начало площадки в ECEF и её широта/долгота. */
    VecEcef originEcef() const
    {
        return m_originEcef;
    }
    double lat() const
    {
        return m_lat;
    }
    double lon() const
    {
        return m_lon;
    }

private:
    VecEcef m_originEcef;
    Mat3 m_Ref; /* ECEF -> ENU */
    Mat3 m_Rfe; /* ENU -> ECEF */
    VecEnu m_gEnu;
    double m_lat, m_lon;
};

#endif /* LIB_EARTH_H */

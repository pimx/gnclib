/* lib_aerodynamics.h -- модели аэродинамических сил.
 *
 * class Aerodynamics вычисляет аэродинамические силы для всех режимов
 * полёта линейки аппаратов:
 *   - сопротивление с постоянным Cd (стек на подъёме, кувыркающийся
 *     корпус);
 *   - посадочный закон сопротивления: Cx(M) с трансзвуковым горбом роста
 *     (~1.70 у M=1), выходящим на сверхзвуковое плато (~1.25), с
 *     опциональным снижением экранированием ретро-струёй (выхлоп навстречу
 *     потоку прикрывает донную область и снижает Cx при работе двигателя);
 *   - маховская поправка осевого Cx (тождественна при M < 0.3 --
 *     дозвуковые аппараты);
 *   - полосовая (strip) модель поперечного обтекания цилиндра: местная
 *     скорость каждой полосы с учётом вращения (v + omega x r) и
 *     квадратичное сопротивление -- автоматически даёт поперечную силу,
 *     опрокидывающий момент от ветра и демпфирование вращения;
 *   - потенциал решётчатых рулей (поперечная сила пары на отклонение при
 *     текущем скоростном напоре);
 *   - экранный эффект (Чизмен-Беннетт): прирост эффективной тяги вблизи
 *     поверхности.
 *
 * Класс держит ссылку на Atmosphere (среда), а параметры конкретного
 * аппарата принимает явными аргументами методов -- один экземпляр
 * Aerodynamics обслуживает любое число тел, полностью переиспользуем.
 * Без STL, без исключений, без динамических аллокаций.
 */
#ifndef LIB_AERODYNAMICS_H
#define LIB_AERODYNAMICS_H

#include "lib_frames.h"
#include "lib_atmosphere.h"

class Aerodynamics
{
public:
    /* Привязка к модели атмосферы (хранится ссылкой; владеет вызывающий). */
    explicit Aerodynamics(const Atmosphere& atmo);

    /* Сопротивление с постоянным Cd [Н], против воздушной скорости:
     * F = -1/2 rho |v| v Cd A. Сила КОЛЛИНЕАРНА воздушной скорости,
     * поэтому даётся ДВУМЯ перегрузками -- в связанной системе и в
     * системе площадки: тип результата совпадает с типом аргумента,
     * смешать системы нельзя. Общее скалярное ядро -- dragMagnitude. */
    VecXBody dragConstantCd(double h, VecXBody vAir, double Cd,
                           double aRef) const;
    VecEnu dragConstantCd(double h, VecEnu vAir, double Cd, double aRef) const;

    /* Модуль силы сопротивления с постоянным Cd [Н] при воздушной
     * скорости V (ядро перегрузок выше; для скалярных прогнозов). */
    double dragMagnitude(double h, double V, double Cd, double aRef) const;

    /* Посадочный коэффициент сопротивления от Маха:
     * Cx(M) = 1.4 + 0.25 exp(-(M-1)^2/0.12) - 0.15 tanh(0.8 (M-1.4)). */
    static double cxMachLanding(double M);

    /* Относительная маховская поправка осевого Cx: 1.0 ниже M = 0.3
     * (несжимаемое обтекание), выше -- рост по закону cxMachLanding,
     * нормированному к дозвуковому значению. */
    static double machFactorAxial(double M);

    /* Снижение посадочного Cx ретро-струёй: Cx_eff = Cx (1 - k thrFrac),
     * пол 5% от Cx; k <= 0 или двигатель выключен -- возврат Cx. */
    static double cxRetroPlume(double Cx, double k, double thrFrac);

    /* Поправка сжимаемости Прандтля-Глауэрта для дозвуковых коэффициентов
     * подъёмной силы: 1/sqrt(1 - M^2); кламп M <= 0.95. */
    static double prandtlGlauert(double M);

    /* Экранный эффект по Чизмену-Беннетту: отношение тяги
     * T_ige/T_oge = 1/(1 - (R/(4h))^2), R -- эффективный радиус диска,
     * h -- высота среза сопел над поверхностью. Кламп сверху maxFactor;
     * при h <= R/4 возврат maxFactor; R <= 0 отключает эффект (1.0). */
    static double groundEffectFactor(double diskRadius, double h,
                                     double maxFactor);

    /* Посадочная сила сопротивления [Н]: Cx(M) со снижением ретро-струёй,
     * против воздушной скорости. thrFrac = тяга/максимум в [0,1].
     * Две перегрузки по системе координат (см. dragConstantCd); общее
     * скалярное ядро -- dragLandingMagnitude. */
    VecXBody dragLanding(double h, VecXBody vAir, double aRef, double retroPlumeK,
                        double thrFrac) const;
    VecEnu dragLanding(double h, VecEnu vAir, double aRef, double retroPlumeK,
                       double thrFrac) const;

    /* Модуль посадочной силы сопротивления [Н] при воздушной скорости V
     * (ядро перегрузок выше). */
    double dragLandingMagnitude(double h, double V, double aRef,
                                double retroPlumeK, double thrFrac) const;

    /* Скалярное посадочное ЗАМЕДЛЕНИЕ [м/с^2] для точечных прогнозов:
     * a = q Cx(M) A / m (со снижением струёй). Для планировщика/борта. */
    double dragDecelLanding(double h, double V, double aRef, double m,
                            double retroPlumeK, double thrFrac) const;

    /* Максимальная поперечная сила ПАРЫ рулей на полном отклонении [Н]
     * при данном режиме (две эффективные консоли на поперечную ось). */
    double finForceMax(double h, double V, double finArea, double finClDelta,
                       double finMaxDeg) const;

    /* Полосовая модель цилиндрического корпуса в осях VecXBody (ось тела
     * X к носу, y и z поперечные). Поперечное обтекание по nStrips полосам
     * + осевое сопротивление по миделю с маховской поправкой и
     * экранированием ретро-струёй при спуске (vRelBody.x < 0).
     *   vRelBody  -- скорость аппарата относительно воздуха в BODY, м/с;
     *   omegaBody -- угловая скорость, рад/с;
     *   rho       -- плотность воздуха, кг/м^3;
     *   mach      -- число Маха воздушной скорости;
     *   thrFrac   -- тяга/максимум [0..1] для экранирования струёй;
     *   length, radius -- геометрия цилиндра, м;
     *   cmFromTail -- ЦМ от хвоста, м;
     *   cdLateral, cdAxial -- коэффициенты сопротивления;
     *   retroPlumeK -- коэффициент экранирования (0 -- выключено);
     *   nStrips   -- число полос (кламп [1..32]).
     * Возвращает силу и момент (относительно ЦМ) в осях BODY. */
    void cylinderStrip(VecXBody vRelBody, VecXBody omegaBody, double rho,
                       double mach, double thrFrac, double length,
                       double radius, double cmFromTail, double cdLateral,
                       double cdAxial, double retroPlumeK, int nStrips,
                       VecXBody* forceBody, VecXBody* momentBody) const;

    /* Доступ к привязанной атмосфере (для прямых запросов rho/Маха). */
    const Atmosphere& atmosphere() const
    {
        return m_atmo;
    }

private:
    const Atmosphere& m_atmo;
};

#endif /* LIB_AERODYNAMICS_H */

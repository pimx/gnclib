//============================================================================
// lib_slosh.h -- плескание топлива: пружинно-массовый аналог в связанных
// осях, связанный с движением твёрдого тела через IntegrateRK4 (1.7.0).
//
// Параметры по заполнению и осевой удельной силе (lib_slosh_params.cpp):
//   цилиндр с плоским днищем радиуса R, глубина h (Abramson NASA SP-106,
//   Dodge SwRI 2000), xi_n -- нули J1'(x):
//     w_n^2 = |a| (xi_n / R) th(xi_n h / R)
//     m_n   = mL 2 R th(xi_n h / R) / (xi_n (xi_n^2 - 1) h)
//     l_n   = h/2 - (2 R / xi_n) th(xi_n h / (2 R))   (выше ЦМ жидкости)
//     m0 = mL - sum m_n, плечо -- из условия ЦМ; I0 -- из импульсного
//     потенциала Стокса-Жуковского (ряд Бесселя-Дини, 100 членов; тона
//     сверх nModes присоединяются к неподвижной массе);
//   демпфирование: гладкая стенка (Stephens и др. / Dodge) -- часть, не
//   зависящая от амплитуды; кольцевые перегородки -- по Майлсу (зависит от
//   амплитуды, первый тон).
//
// Динамика (lib_slosh_eval.cpp). Интегрируется ЖЁСТКАЯ ЧАСТЬ (сухая
// конструкция + неподвижные массы m0 всех баков, AssembleRigid); её ЦМ --
// ЦМ состояния VehState. Масса плескания j (плечо r_j = покой + xi_j,
// rho_j = r_j - cm):
//   xi_j'' = -w_j^2 xi_j - 2 zeta_j w_j v_j
//            - [ f_r + eps x rho_j + w x (w x rho_j) + 2 w x v_j ]_поперёк
// f_r -- удельная сила ЦМ жёсткой части. Пружина и демпфер действуют на
// жёсткую часть поперёк в точке r_j; вдоль оси массу несёт бак:
//   N_j = m_j [ f_r + eps x rho_j + w x (w x rho_j) + 2 w x v_j ]_bu,
// отсюда момент m a xi маятника. N_j зависит от ускорений; зависимость
// разрешена ТОЧНО (осевые суммы -- по полной массе, часть от eps -- как
// приращение инерции dJ в возвращаемом INR с гироскопической поправкой
// момента), поэтому ускорения, которые IntegrateRK4 считает по
// возвращённым силе, моменту и INR, -- те же, что взяты в производных
// плескания. Импульс, момент импульса и энергия (zeta = 0) связанной
// системы сохраняются с точностью интегрирования (tests/test_slosh.cpp,
// проверка 8).
//
// Применение:
//   model.Init(tanks, nTank, auxOffset); model.ResetAux(&aux);
//   на каждом шаге: model.Update(mLiquid, aAxial), aAxial = F_ng.bu / M_полн;
//              model.AssembleRigid(dryInr, dryCm, &rigidInr, &rigidCm);
//              state.SetCmArm(rigidCm);
//   в функции производных fn(ctx, s, aux, d):
//     SloshIn in: fNg -- тяга + аэродинамика (связанные, БЕЗ тяжести),
//     mExt -- их момент относительно ЦМ ЖЁСТКОЙ части, rigid = rigidInr,
//     cm = s.CmArm(), w = s.OmgBdw();
//     model.Eval(in, aux, &out, d->xdot);
//     d->f = B2L(fNg + out.force) + rigidInr.m * g_площадки;
//     d->mom = out.mom (mExt уже внутри); d->inr = out.inr.
// Модель линейна по поперечному смещению (малые амплитуды); при малой |a|
// (пассивный участок, число Бонда) неприменима -- признак valid.
// Подключается явно (не через lib.h). C++17, без STL, без кучи, без
// исключений, без виртуальных функций.
//============================================================================
#ifndef GNCLIB_LIB_SLOSH_H
#define GNCLIB_LIB_SLOSH_H

#include "lib_slosh_types.h"
#include "lib_integrate.h"

namespace gnc
{

//-- параметры (lib_slosh_params.cpp) ----------------------------------------

// n-й нуль J1'(x) (n = 1, 2, ...): таблица до n = 10, далее разложение Макмагона.
//   n -- номер нуля
double SloshBesselZero(int n);

// Приведённый аналог цилиндра с плоским днищем при глубине жидкости h.
//   radius -- R, м; depth -- h, м (> 0); nModes -- явных тонов (1..3,
//   лишние записи обнулены); row -- выход
// Возвращает false при неверных аргументах.
bool SloshCylinderRow(double radius, double depth, int nModes, SloshRow* row);

// Таблица строк цилиндра с плоским днищем на всё заполнение.
//   radius, length -- R, L, м; nModes -- хранимых тонов (1..3);
//   nRows -- 2..SLOSH_NROW; tab -- выход (q гуще при малом заполнении)
// Возвращает false при неверных аргументах.
bool SloshTableCylinder(double radius, double length, int nModes, int nRows, SloshTable* tab);

// Строка таблицы при доле заполнения q (линейно, края удерживаются).
//   tab -- таблица (n >= 1); q -- доля заполнения; row -- выход
void SloshTableRow(const SloshTable& tab, double q, SloshRow* row);

// Параметры бака по приведённой строке: тона сверх tank.nModes
// присоединяются к неподвижной массе, масштаб по mL, положения по знаку a,
// демпфирование и признак применимости.
//   tank -- описание; row -- приведённый аналог; mL -- масса жидкости, кг;
//   a -- осевая удельная сила (со знаком; |a| не меньше aMin), м/с^2;
//   aMin, bondMin -- пределы применимости; p -- выход
void SloshRowToParams(const SloshTank& tank, const SloshRow& row, double mL, double a,
                      double aMin, double bondMin, SloshTankParams* p);

// Демпфирование первого тона гладкой стенкой (Stephens и др., Dodge 2000):
// zeta = 0.79 sqrt(nu) |a|^-1/4 R^-3/4 [1 + 0.318 / sh(1.84 h/R)
// (1 + (1 - h/R) / ch(1.84 h/R))].
//   nu -- кинематическая вязкость, м^2/с; a -- |a|, м/с^2; radius, depth -- м
double SloshZetaWall(double nu, double a, double radius, double depth);

// Демпфирование кольцевой перегородкой по Майлсу: zeta = 2.83 exp(-4.6 d/R)
// beta^1.5 sqrt(wave / R), beta -- перекрытая доля сечения.
//   depth -- d, заглубление пластины под свободной поверхностью, м (<= 0 -- 0);
//   width -- ширина кольца, м; radius -- R, м; wave -- высота волны у стенки, м
double SloshZetaBaffle(double depth, double width, double radius, double wave);

//-- таблицы в CSV (lib_slosh_io.cpp) ----------------------------------------

// Запись таблицы: строка заголовка, затем столбцы через запятую
// q,zSurf,zLiq,m0r,z0,i0r,iRollr,mr1..3,z1..3,w2a1..3,eta1..3.
//   path -- файл; tab -- таблица; err, errSize -- причина отказа (ASCII)
// Возвращает false, если файл не записывается.
bool SloshTableSaveCsv(const char* path, const SloshTable& tab, char* err, int errSize);

// Чтение таблицы в том же порядке столбцов (через запятые или пробелы;
// пустые строки, строки с ';' или '#' и заголовок до первой строки чисел
// пропускаются).
//   path -- файл; out -- выход (меняется только при успехе);
//   err, errSize -- причина отказа (ASCII)
// Возвращает false, если файл не читается, строка не из 19 чисел, q не
// растёт, строк нет или больше SLOSH_NROW.
bool SloshTableLoadCsv(const char* path, SloshTable* out, char* err, int errSize);

//-- связанная модель (lib_slosh_model.cpp, lib_slosh_eval.cpp) --------------

// Вход SloshModel::Eval на стадии.
struct SloshIn
{
    CfBodyVector fNg { };       // внешняя негравитационная сила, связанные, Н
    CfBodyMoment mExt { };      // её момент относительно ЦМ жёсткой части, Н м
    INR rigid { };              // жёсткая часть (AssembleRigid), относительно её ЦМ
    CfBodyVector cm { };        // плечо ЦМ жёсткой части (s.CmArm()), м
    CfBodyOmega w { };          // угловая скорость (s.OmgBdw()), рад/с
};

// Выход SloshModel::Eval.
struct SloshOut
{
    CfBodyVector force { };     // сила плескания на жёсткую часть, связанные, Н
    CfBodyMoment mom { };       // ПОЛНЫЙ момент для DynDeriv (mExt + плескание + dJ)
    INR inr { };                // действующий INR для DynDeriv (жёсткая часть + dJ)
    CfBodyOmegaAcc eps { };     // угловое ускорение, рад/с^2
    CfBodyVector fRigid { };    // удельная сила ЦМ жёсткой части, связанные, м/с^2
};

class SloshModel
{
public:
    static constexpr int NTANK = 4;
    static constexpr int NJ = NTANK * SLOSH_NMODE;
    static constexpr int NAUX_MODE = 4; // xi_bn, xi_be, v_bn, v_be

    // Пустая модель (без баков), пределы aMin = 0.1 м/с^2, bondMin = 100.
    SloshModel();

    // Настройка баков.
    //   tanks -- описания (копируются; таблицы -- по указателю, должны
    //   жить); nTank -- 1..NTANK; auxOffset -- первое вспомогательное состояние
    // Возвращает число занятых вспомогательных состояний; -1 при неверных
    // аргументах или если auxOffset + 4 * (сумма тонов) больше AuxState::NMAX.
    int Init(const SloshTank* tanks, int nTank, int auxOffset);

    // Пределы применимости: |a| меньше aMin или число Бонда меньше bondMin
    // (при sigma > 0) снимают признак valid; |a| поднимается до aMin.
    //   aMin -- м/с^2; bondMin -- число Бонда
    void SetLimits(double aMin, double bondMin);

    // Обнуление состояний плескания; aux->n поднимается до их конца.
    //   aux -- вспомогательные состояния
    void ResetAux(AuxState* aux) const;

    // Параметры всех баков на шаг.
    //   mLiquid -- масса жидкости каждого бака, кг; aAxial -- осевая
    //   удельная сила F_ng.bu / M_полн, м/с^2 (> 0 при тяге)
    void Update(const double* mLiquid, double aAxial);

    // Жёсткая часть = сухая часть + неподвижные массы всех баков (по
    // теореме о параллельных осях).
    //   dry, dryCm -- сухая часть и плечо её ЦМ; rigid, rigidCm -- выход
    void AssembleRigid(const INR& dry, const CfBodyVector& dryCm, INR* rigid,
                       CfBodyVector* rigidCm) const;

    // Силы плескания и производные на стадии (см. шапку файла).
    //   in -- вход стадии; aux -- вспомогательные состояния; out -- выход;
    //   xdot -- производные всех вспомогательных состояний (пишутся только
    //   состояния плескания)
    void Eval(const SloshIn& in, const AuxState& aux, SloshOut* out, double* xdot) const;

    // Сумма масс плескания, кг (полная масса = жёсткая + MassSlosh()).
    double MassSlosh() const;

    // ЦМ всей системы (жёсткая часть + смещённые массы плескания).
    //   aux -- состояния; mRigid, rigidCm -- жёсткая часть
    // Возвращает плечо, м.
    CfBodyVector SystemCm(const AuxState& aux, double mRigid, const CfBodyVector& rigidCm) const;

    // Положение массы плескания в покое (плечо по bu и смещение оси бака), м.
    //   tank, mode -- номера (нулевой вектор, если тона нет)
    CfBodyVector RestArm(int tank, int mode) const;

    // Телеметрия тона.
    //   aux -- состояния; tank, mode -- номера; st -- выход
    // Возвращает false, если тона нет.
    bool ModeState(const AuxState& aux, int tank, int mode, SloshModeState* st) const;

    // Действующее демпфирование тона: zeta0 и для первого тона --
    // затопленные кольцевые перегородки при текущей амплитуде (Майлс), <= 1.
    //   aux -- состояния; tank, mode -- номера
    double ModeZeta(const AuxState& aux, int tank, int mode) const;

    // Число баков.
    int NumTanks() const
    {
        return m_nTank;
    }

    // Число вспомогательных состояний тонов плескания.
    int NumAux() const
    {
        return m_nAux;
    }

    // Первое вспомогательное состояние тона (далее xi_bn, xi_be, v_bn, v_be).
    //   tank, mode -- номера
    // Возвращает -1, если тона нет.
    int AuxIndex(int tank, int mode) const;

    // Описание бака.
    //   i -- номер бака, 0..NumTanks() - 1
    const SloshTank& Tank(int i) const
    {
        return m_tank[i];
    }

    // Параметры бака по последнему Update.
    //   i -- номер бака, 0..NumTanks() - 1
    const SloshTankParams& Params(int i) const
    {
        return m_par[i];
    }

private:
    SloshTank m_tank[NTANK];
    SloshTankParams m_par[NTANK];
    int m_nTank;
    int m_aux0;
    int m_nAux;
    double m_aMin;
    double m_bondMin;
};

} // namespace gnc

#endif // GNCLIB_LIB_SLOSH_H

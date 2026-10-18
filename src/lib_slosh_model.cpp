//============================================================================
// lib_slosh_model.cpp -- SloshModel (1.7.0): настройка, параметры шага,
// жёсткая часть, телеметрия. Связанная динамика (Eval) -- в
// lib_slosh_eval.cpp. См. lib_slosh.h.
//============================================================================
#include "lib_slosh.h"

#include <math.h>

#include "lib_util.h"

namespace gnc
{

namespace
{

// Тела, накопленные относительно начала связанной системы: статические
// моменты и тензор инерции относительно начала (элементы тензора, как в INR).
struct MassAcc
{
    double m = 0.0;
    double sBn = 0.0;
    double sBe = 0.0;
    double sBu = 0.0;
    double jNn = 0.0;
    double jEe = 0.0;
    double jUu = 0.0;
    double jNe = 0.0;
    double jNu = 0.0;
    double jEu = 0.0;
};

// Добавление тела: масса m, плечо ЦМ c, тензор jOwn относительно своего ЦМ.
void AccAdd(MassAcc* a, double m, const CfBodyVector& c, const INR& jOwn)
{
    a->m += m;
    a->sBn += m * c.bn;
    a->sBe += m * c.be;
    a->sBu += m * c.bu;
    a->jNn += jOwn.j_bnbn + m * (c.be * c.be + c.bu * c.bu);
    a->jEe += jOwn.j_bebe + m * (c.bn * c.bn + c.bu * c.bu);
    a->jUu += jOwn.j_bubu + m * (c.bn * c.bn + c.be * c.be);
    a->jNe += jOwn.j_bnbe - m * c.bn * c.be;
    a->jNu += jOwn.j_bnbu - m * c.bn * c.bu;
    a->jEu += jOwn.j_bebu - m * c.be * c.bu;
}

// Массово-инерционные характеристики относительно общего ЦМ.
void AccResult(const MassAcc& a, INR* inr, CfBodyVector* cm)
{
    INR r { };
    CfBodyVector c { };
    r.m = a.m;
    if (a.m > 0.0)
    {
        c.bn = a.sBn / a.m;
        c.be = a.sBe / a.m;
        c.bu = a.sBu / a.m;
    }
    r.j_bnbn = a.jNn - a.m * (c.be * c.be + c.bu * c.bu);
    r.j_bebe = a.jEe - a.m * (c.bn * c.bn + c.bu * c.bu);
    r.j_bubu = a.jUu - a.m * (c.bn * c.bn + c.be * c.be);
    r.j_bnbe = a.jNe + a.m * c.bn * c.be;
    r.j_bnbu = a.jNu + a.m * c.bn * c.bu;
    r.j_bebu = a.jEu + a.m * c.be * c.bu;
    *inr = r;
    *cm = c;
}

} // namespace

//----------------------------------------------------------------------------
// Пустая модель с пределами применимости по умолчанию.
//----------------------------------------------------------------------------
SloshModel::SloshModel()
    : m_tank()
    , m_par()
    , m_nTank(0)
    , m_aux0(0)
    , m_nAux(0)
    , m_aMin(0.1)
    , m_bondMin(100.0)
{
}

//----------------------------------------------------------------------------
// Настройка (см. lib_slosh.h). Неположительные размеры или плотность --
// отказ.
//   tanks -- описания; nTank -- 1..NTANK; auxOffset -- первое состояние
//----------------------------------------------------------------------------
int SloshModel::Init(const SloshTank* tanks, int nTank, int auxOffset)
{
    if (tanks == 0 || nTank < 1 || nTank > NTANK || auxOffset < 0)
    {
        return -1;
    }
    int nAux = 0;
    for (int t = 0; t < nTank; ++t)
    {
        const SloshTank& k = tanks[t];
        if (k.radius <= 0.0 || k.length <= 0.0 || k.rho <= 0.0)
        {
            return -1;
        }
        if (k.nModes < 1 || k.nModes > SLOSH_NMODE)
        {
            return -1;
        }
        if (k.nBaffle < 0 || k.nBaffle > SLOSH_NBAFFLE)
        {
            return -1;
        }
        nAux += NAUX_MODE * k.nModes;
    }
    if (auxOffset + nAux > AuxState::NMAX)
    {
        return -1;
    }
    for (int t = 0; t < NTANK; ++t)
    {
        m_par[t] = SloshTankParams();
    }
    for (int t = 0; t < nTank; ++t)
    {
        m_tank[t] = tanks[t];
    }
    m_nTank = nTank;
    m_aux0 = auxOffset;
    m_nAux = nAux;
    return nAux;
}

//----------------------------------------------------------------------------
// Пределы применимости (см. lib_slosh.h).
//   aMin -- м/с^2 (> 0); bondMin -- число Бонда
//----------------------------------------------------------------------------
void SloshModel::SetLimits(double aMin, double bondMin)
{
    m_aMin = (aMin > 1.0e-6) ? aMin : 1.0e-6;
    m_bondMin = bondMin;
}

//----------------------------------------------------------------------------
// Обнуление состояний плескания (см. lib_slosh.h).
//   aux -- вспомогательные состояния
//----------------------------------------------------------------------------
void SloshModel::ResetAux(AuxState* aux) const
{
    if (aux == 0)
    {
        return;
    }
    for (int i = m_aux0; i < m_aux0 + m_nAux; ++i)
    {
        aux->x[i] = 0.0;
    }
    if (aux->n < m_aux0 + m_nAux)
    {
        aux->n = m_aux0 + m_nAux;
    }
}

//----------------------------------------------------------------------------
// Параметры шага: строка цилиндра при глубине h = mL / (rho pi R^2)
// (в пределах 1e-4 R .. 0.999 L) или строка таблицы при q = mL / massFull.
//   mLiquid -- кг на бак; aAxial -- м/с^2
//----------------------------------------------------------------------------
void SloshModel::Update(const double* mLiquid, double aAxial)
{
    for (int t = 0; t < m_nTank; ++t)
    {
        const SloshTank& k = m_tank[t];
        double mL = (mLiquid != 0) ? mLiquid[t] : 0.0;
        SloshRow row;
        if (k.table != 0)
        {
            double mFull = k.massFull;
            if (mFull <= 0.0)
            {
                mFull = k.rho * PI * k.radius * k.radius * k.length;
            }
            SloshTableRow(*k.table, (mFull > 0.0) ? mL / mFull : 0.0, &row);
        }
        else
        {
            double h = mL / (k.rho * PI * k.radius * k.radius);
            h = Clamp(h, 1.0e-4 * k.radius, 0.999 * k.length);
            SloshCylinderRow(k.radius, h, k.nModes, &row);
        }
        SloshRowToParams(k, row, mL, aAxial, m_aMin, m_bondMin, &m_par[t]);
    }
}

//----------------------------------------------------------------------------
// Жёсткая часть (см. lib_slosh.h): неподвижная масса бака -- на оси бака,
// I0 вокруг bn и be, kRoll * замороженный момент крена вокруг bu.
//   dry, dryCm -- сухая часть; rigid, rigidCm -- выход
//----------------------------------------------------------------------------
void SloshModel::AssembleRigid(const INR& dry, const CfBodyVector& dryCm, INR* rigid,
                               CfBodyVector* rigidCm) const
{
    MassAcc acc;
    AccAdd(&acc, dry.m, dryCm, dry);
    for (int t = 0; t < m_nTank; ++t)
    {
        const SloshTankParams& p = m_par[t];
        if (p.m0 <= 0.0)
        {
            continue;
        }
        INR j { };
        j.m = p.m0;
        j.j_bnbn = p.i0;
        j.j_bebe = p.i0;
        j.j_bubu = p.iRoll;
        CfBodyVector c { };
        c.bn = m_tank[t].offBn;
        c.be = m_tank[t].offBe;
        c.bu = p.arm0;
        AccAdd(&acc, p.m0, c, j);
    }
    AccResult(acc, rigid, rigidCm);
}

//----------------------------------------------------------------------------
// Сумма масс плескания, кг.
//----------------------------------------------------------------------------
double SloshModel::MassSlosh() const
{
    double m = 0.0;
    for (int t = 0; t < m_nTank; ++t)
    {
        for (int n = 0; n < m_par[t].nModes; ++n)
        {
            m += m_par[t].mode[n].m;
        }
    }
    return m;
}

//----------------------------------------------------------------------------
// Первое вспомогательное состояние тона (тона уложены по бакам подряд),
// -1, если тона нет.
//   tank, mode -- номера
//----------------------------------------------------------------------------
int SloshModel::AuxIndex(int tank, int mode) const
{
    if (tank < 0 || tank >= m_nTank || mode < 0 || mode >= m_tank[tank].nModes)
    {
        return -1;
    }
    int idx = m_aux0;
    for (int t = 0; t < tank; ++t)
    {
        idx += NAUX_MODE * m_tank[t].nModes;
    }
    return idx + NAUX_MODE * mode;
}

//----------------------------------------------------------------------------
// Плечо массы плескания в покое (нулевой вектор, если тона нет).
//   tank, mode -- номера
//----------------------------------------------------------------------------
CfBodyVector SloshModel::RestArm(int tank, int mode) const
{
    CfBodyVector r { };
    if (AuxIndex(tank, mode) < 0)
    {
        return r;
    }
    r.bn = m_tank[tank].offBn;
    r.be = m_tank[tank].offBe;
    r.bu = m_par[tank].mode[mode].arm;
    return r;
}

//----------------------------------------------------------------------------
// Действующее демпфирование (см. lib_slosh.h). Заглубление перегородки под
// поверхностью: sgn (armSurf - arm) -- положительно для затопленной
// пластины при обоих положениях жидкости.
//   aux -- состояния; tank, mode -- номера
//----------------------------------------------------------------------------
double SloshModel::ModeZeta(const AuxState& aux, int tank, int mode) const
{
    int idx = AuxIndex(tank, mode);
    if (idx < 0)
    {
        return 0.0;
    }
    const SloshTank& k = m_tank[tank];
    const SloshTankParams& p = m_par[tank];
    const SloshModeParams& md = p.mode[mode];
    double z = md.zeta0;
    if (mode == 0 && k.nBaffle > 0 && md.w > 0.0)
    {
        double xi2 = aux.x[idx] * aux.x[idx] + aux.x[idx + 1] * aux.x[idx + 1];
        double v2 = aux.x[idx + 2] * aux.x[idx + 2] + aux.x[idx + 3] * aux.x[idx + 3];
        double wave = md.eta * sqrt(xi2 + v2 / (md.w * md.w));
        for (int b = 0; b < k.nBaffle; ++b)
        {
            double depth = p.sgn * (p.armSurf - k.baffle[b].arm);
            z += SloshZetaBaffle(depth, k.baffle[b].width, k.radius, wave);
        }
    }
    return (z > 1.0) ? 1.0 : z;
}

//----------------------------------------------------------------------------
// Телеметрия тона (см. lib_slosh.h).
//   aux -- состояния; tank, mode -- номера; st -- выход
//----------------------------------------------------------------------------
bool SloshModel::ModeState(const AuxState& aux, int tank, int mode, SloshModeState* st) const
{
    int idx = AuxIndex(tank, mode);
    if (idx < 0 || st == 0)
    {
        return false;
    }
    const SloshModeParams& md = m_par[tank].mode[mode];
    SloshModeState s;
    s.xiBn = aux.x[idx];
    s.xiBe = aux.x[idx + 1];
    s.vBn = aux.x[idx + 2];
    s.vBe = aux.x[idx + 3];
    double v2 = s.vBn * s.vBn + s.vBe * s.vBe;
    double w2 = (md.w > 0.0) ? md.w * md.w : 1.0;
    s.amp = sqrt(s.xiBn * s.xiBn + s.xiBe * s.xiBe + v2 / w2);
    s.wave = md.eta * s.amp;
    s.zeta = ModeZeta(aux, tank, mode);
    *st = s;
    return true;
}

//----------------------------------------------------------------------------
// ЦМ всей системы (см. lib_slosh.h).
//   aux -- состояния; mRigid, rigidCm -- жёсткая часть
//----------------------------------------------------------------------------
CfBodyVector SloshModel::SystemCm(const AuxState& aux, double mRigid,
                                  const CfBodyVector& rigidCm) const
{
    double m = mRigid;
    CfBodyVector s = Scale(rigidCm, mRigid);
    for (int t = 0; t < m_nTank; ++t)
    {
        for (int n = 0; n < m_par[t].nModes; ++n)
        {
            int idx = AuxIndex(t, n);
            double mj = m_par[t].mode[n].m;
            CfBodyVector r = RestArm(t, n);
            r.bn += aux.x[idx];
            r.be += aux.x[idx + 1];
            s = Add(s, Scale(r, mj));
            m += mj;
        }
    }
    if (m > 0.0)
    {
        s = Scale(s, 1.0 / m);
    }
    return s;
}

} // namespace gnc

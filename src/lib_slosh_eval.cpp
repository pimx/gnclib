//============================================================================
// lib_slosh_eval.cpp -- SloshModel::Eval (1.7.0): силы масс плескания на
// жёсткую часть и производные состояний плескания на стадии; осевая связь
// разрешается точно. См. lib_slosh.h.
//============================================================================
#include "lib_slosh.h"

namespace gnc
{

namespace
{

// Рабочие данные одной массы плескания на стадии.
struct Work
{
    int idx = 0;            // первое вспомогательное состояние
    double m = 0.0;         // масса, кг
    double k = 0.0;         // жёсткость, Н/м
    double c = 0.0;         // демпфер, Н с/м (2 zeta w m)
    double w2 = 0.0;        // w^2, 1/с^2
    double cz = 0.0;        // 2 zeta w, 1/с
    double xiBn = 0.0;      // смещение от покоя, м
    double xiBe = 0.0;
    double vBn = 0.0;       // скорость относительно бака, м/с
    double vBe = 0.0;
    CfBodyVector rho { };   // плечо от ЦМ жёсткой части, м
    CfBodyVector wwr { };   // w x (w x rho), м/с^2
    CfBodyVector wv2 { };   // 2 w x v, м/с^2
    double gBn = 0.0;       // (eps x rho).bu = gBn eps.bn + gBe eps.be
    double gBe = 0.0;
    double n0 = 0.0;        // осевая сила на массу без части от eps, Н
};

// a x r в связанных осях, правая тройка (bu, be, bn): bu x be = bn.
CfBodyVector CrossB(double aBn, double aBe, double aBu, const CfBodyVector& r)
{
    CfBodyVector v { };
    v.bu = aBe * r.bn - aBn * r.be;
    v.be = aBn * r.bu - aBu * r.bn;
    v.bn = aBu * r.be - aBe * r.bu;
    return v;
}

// Рабочие данные всех масс плескания; возвращает их число.
//   md -- модель; in -- вход стадии; aux -- состояния; wk -- выход
int Gather(const SloshModel& md, const SloshIn& in, const AuxState& aux, Work* wk)
{
    int nj = 0;
    const CfBodyOmega& w = in.w;
    for (int t = 0; t < md.NumTanks(); ++t)
    {
        const SloshTankParams& p = md.Params(t);
        for (int n = 0; n < p.nModes; ++n)
        {
            Work& j = wk[nj];
            j.idx = md.AuxIndex(t, n);
            const SloshModeParams& mp = p.mode[n];
            double zeta = md.ModeZeta(aux, t, n);
            j.m = mp.m;
            j.k = mp.k;
            j.w2 = mp.w * mp.w;
            j.cz = 2.0 * zeta * mp.w;
            j.c = j.cz * mp.m;
            j.xiBn = aux.x[j.idx];
            j.xiBe = aux.x[j.idx + 1];
            j.vBn = aux.x[j.idx + 2];
            j.vBe = aux.x[j.idx + 3];
            CfBodyVector r = md.RestArm(t, n);
            r.bn += j.xiBn;
            r.be += j.xiBe;
            j.rho = Sub(r, in.cm);
            CfBodyVector wr = CrossB(w.bn, w.be, w.bu, j.rho);
            j.wwr = CrossB(w.bn, w.be, w.bu, wr);
            CfBodyVector v { };
            v.bn = j.vBn;
            v.be = j.vBe;
            j.wv2 = Scale(CrossB(w.bn, w.be, w.bu, v), 2.0);
            j.gBn = -j.rho.be;
            j.gBe = j.rho.bn;
            j.n0 = 0.0;
            ++nj;
        }
    }
    return nj;
}

} // namespace

//----------------------------------------------------------------------------
// Силы и производные на стадии (см. lib_slosh.h).
// Алгоритм. N_j -- осевая сила бака на массу j (вдоль bu):
//   N_j = m_j [ f_r.bu + (eps x rho_j).bu + X_j ],  X_j = (w x (w x rho_j)
//   + 2 w x v_j).bu,  (eps x rho_j).bu = g_j . eps,  g_j = (-rho.be, rho.bn).
// Осевое уравнение жёсткой части даёт f_r.bu через ПОЛНУЮ массу M:
//   f_r.bu = (F.bu - sum m X) / M - gm . eps,  gm = sum m g / M, откуда
//   N_j = N0_j + m_j (g_j - gm) . eps,  N0_j = m_j (F.bu / M - sum m X / M + X_j).
// Часть момента сил -N_j в точках rho_j от eps равна -dJ eps,
//   dJ = sum m g g^T - M gm gm^T (блок bn-be, неотрицательно определён),
// так что (J + dJ) eps = M_ext + M0 - w x J w; в IntegrateRK4 уходят INR
// J + dJ и момент M_ext + M0 + w x (dJ w) (EpsOf вычитает w x (J + dJ) w).
// Здесь eps находится тем же EpsOf, затем -- точные N_j, сила плескания,
// f_r и поперечные ускорения масс.
//   in -- вход стадии; aux -- состояния; out -- выход; xdot -- производные
//----------------------------------------------------------------------------
void SloshModel::Eval(const SloshIn& in, const AuxState& aux, SloshOut* out, double* xdot) const
{
    Work wk[NJ];
    int nj = Gather(*this, in, aux, wk);
    double mR = (in.rigid.m > 1.0e-9) ? in.rigid.m : 1.0e-9;
    double mT = mR;
    double sx = 0.0;
    double gmBn = 0.0;
    double gmBe = 0.0;
    for (int i = 0; i < nj; ++i)
    {
        mT += wk[i].m;
        sx += wk[i].m * (wk[i].wwr.bu + wk[i].wv2.bu);
        gmBn += wk[i].m * wk[i].gBn;
        gmBe += wk[i].m * wk[i].gBe;
    }
    gmBn /= mT;
    gmBe /= mT;
    double fax0 = (in.fNg.bu - sx) / mT;
    CfBodyVector fs { };
    CfBodyMoment m0 = in.mExt;
    double dNn = -mT * gmBn * gmBn;
    double dEe = -mT * gmBe * gmBe;
    double dNe = -mT * gmBn * gmBe;
    for (int i = 0; i < nj; ++i)
    {
        Work& j = wk[i];
        j.n0 = j.m * (fax0 + j.wwr.bu + j.wv2.bu);
        CfBodyVector r { };
        r.bn = j.k * j.xiBn + j.c * j.vBn;
        r.be = j.k * j.xiBe + j.c * j.vBe;
        r.bu = -j.n0;
        fs.bn += r.bn;
        fs.be += r.be;
        m0 = Add(m0, Torque(j.rho, r));
        dNn += j.m * j.gBn * j.gBn;
        dEe += j.m * j.gBe * j.gBe;
        dNe += j.m * j.gBn * j.gBe;
    }
    INR je = in.rigid;
    je.j_bnbn += dNn;
    je.j_bebe += dEe;
    je.j_bnbe += dNe;
    const CfBodyOmega& w = in.w;
    CfBodyVector h { };
    h.bn = dNn * w.bn + dNe * w.be;
    h.be = dNe * w.bn + dEe * w.be;
    CfBodyVector wh = CrossB(w.bn, w.be, w.bu, h);
    CfBodyMoment mom = m0;
    mom.bn += wh.bn;
    mom.be += wh.be;
    mom.bu += wh.bu;
    CfBodyOmegaAcc eps = EpsOf(je, mom, w);
    for (int i = 0; i < nj; ++i)
    {
        const Work& j = wk[i];
        double nAx = j.n0 + j.m * ((j.gBn - gmBn) * eps.bn + (j.gBe - gmBe) * eps.be);
        fs.bu -= nAx;
    }
    CfBodyVector fr = Scale(Add(in.fNg, fs), 1.0 / mR);
    for (int i = 0; i < nj; ++i)
    {
        const Work& j = wk[i];
        CfBodyVector er = CrossB(eps.bn, eps.be, eps.bu, j.rho);
        double aBn = fr.bn + er.bn + j.wwr.bn + j.wv2.bn;
        double aBe = fr.be + er.be + j.wwr.be + j.wv2.be;
        xdot[j.idx] = j.vBn;
        xdot[j.idx + 1] = j.vBe;
        xdot[j.idx + 2] = -j.w2 * j.xiBn - j.cz * j.vBn - aBn;
        xdot[j.idx + 3] = -j.w2 * j.xiBe - j.cz * j.vBe - aBe;
    }
    out->force = fs;
    out->mom = mom;
    out->inr = je;
    out->eps = eps;
    out->fRigid = fr;
}

} // namespace gnc

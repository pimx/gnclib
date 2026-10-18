// test_slosh_sim.h -- gnclib 1.7.0: обвязка проверок плескания -- изделие
// (сухая часть + баки), функция производных для IntegrateRK4, инварианты
// (энергия, импульс, момент импульса), частота плескания по нулям.
#ifndef TEST_SLOSH_SIM_H
#define TEST_SLOSH_SIM_H

#include "lib_slosh.h"

#include <math.h>
#include <stdio.h>

namespace sltest
{

using namespace gnc;

const double G0 = 9.80665;

// Число непройденных проверок программы.
inline int& BadCount()
{
    static int bad = 0;
    return bad;
}

// Строка проверки: |val| <= lim -- OK.
//   name -- описание (ASCII); val -- значение; lim -- предел
inline void Check(const char* name, double val, double lim)
{
    bool ok = fabs(val) <= lim;
    printf("%-60s %12.3e  (limit %.1e)  %s\n", name, val, lim, ok ? "OK" : "FAIL");
    BadCount() += ok ? 0 : 1;
}

// dir/name в buf.
//   buf, size -- выход; dir -- каталог; name -- имя файла
inline const char* PathJoin(char* buf, int size, const char* dir, const char* name)
{
    snprintf(buf, (size_t)size, "%s/%s", dir, name);
    return buf;
}

// Пример бака: R = 1.5 м, L = 8 м, заднее днище в 6 м от хвоста, жидкость
// плотности 1141 кг/м^3, 3 тона, демпфирование гладкой стенкой, без перегородок.
inline SloshTank TankExample()
{
    SloshTank k;
    k.radius = 1.5;
    k.length = 8.0;
    k.armAft = 6.0;
    k.rho = 1141.0;
    k.nu = 1.7e-7;
    k.sigma = 0.0132;
    k.nModes = 3;
    return k;
}

// Полная масса жидкости бака, кг.
//   k -- бак
inline double MassFull(const SloshTank& k)
{
    return k.rho * PI * k.radius * k.radius * k.length;
}

// Изделие проверок: сухая часть, баки, тяга вдоль +bu в хвосте.
struct Vehicle
{
    SloshModel model;
    INR dry { };
    CfBodyVector dryCm { };
    INR rigid { };
    CfBodyVector rigidCm { };
    double thrust = 0.0;
};

// Параметры шага и жёсткая часть.
//   v -- изделие; mL -- массы жидкости, кг; a -- осевая удельная сила, м/с^2
inline void VehUpdate(Vehicle* v, const double* mL, double a)
{
    v->model.Update(mL, a);
    v->model.AssembleRigid(v->dry, v->dryCm, &v->rigid, &v->rigidCm);
}

// Функция производных: тяга в хвосте (начало связанной системы), без тяжести.
inline void FnVeh(void* c, const VehState& s, const AuxState& aux, DynDeriv* d)
{
    Vehicle* v = (Vehicle*)c;
    SloshIn in;
    in.fNg.bu = v->thrust;
    CfBodyVector tail { };
    in.mExt = MomentAboutCm(in.fNg, tail, s.CmArm());
    in.rigid = v->rigid;
    in.cm = s.CmArm();
    in.w = s.OmgBdw();
    SloshOut out;
    v->model.Eval(in, aux, &out, d->xdot);
    d->f = s.Att().B2L(Add(in.fNg, out.force));
    d->mom = out.mom;
    d->inr = out.inr;
}

// Начальное состояние: ЦМ жёсткой части в начале площадки, единичная ориентация.
//   v -- изделие; w -- угловая скорость; vel -- скорость, оси площадки
inline VehState InitState(const Vehicle& v, const CfBodyOmega& w, const CfLocalVector& vel)
{
    VehState s;
    CfLocalPoint p { };
    CfLocalVector a { };
    s.SetLinear(p, vel, a);
    CfBodyOmegaAcc e { };
    s.SetAngular(s.Att(), w, e);
    s.SetCmArm(v.rigidCm);
    return s;
}

// a x b в осях площадки, правая тройка (u, e, n).
inline CfLocalVector CrossL(const CfLocalVector& a, const CfLocalVector& b)
{
    CfLocalVector c { };
    c.n = a.u * b.e - a.e * b.u;
    c.u = a.e * b.n - a.n * b.e;
    c.e = a.n * b.u - a.u * b.n;
    return c;
}

struct Invariants
{
    double energy = 0.0;        // кинетическая + энергия пружин, Дж
    CfLocalVector p { };        // импульс, кг м/с
    CfLocalVector h { };        // момент импульса относительно начала площадки
};

// Инварианты всей системы (жёсткая часть + массы плескания).
//   v -- изделие; s -- состояние; aux -- состояния плескания
inline Invariants Invar(const Vehicle& v, const VehState& s, const AuxState& aux)
{
    const INR& j = v.rigid;
    const CfBodyOmega& w = s.OmgBdw();
    CfBodyVector hb { };
    hb.bn = j.j_bnbn * w.bn + j.j_bnbe * w.be + j.j_bnbu * w.bu;
    hb.be = j.j_bnbe * w.bn + j.j_bebe * w.be + j.j_bebu * w.bu;
    hb.bu = j.j_bnbu * w.bn + j.j_bebu * w.be + j.j_bubu * w.bu;
    CfLocalVector rr { };
    rr.n = s.PosLop().n;
    rr.e = s.PosLop().e;
    rr.u = s.PosLop().u;
    const CfLocalVector& vr = s.VelLoc();
    Invariants r;
    r.energy = 0.5 * j.m * Norm2(vr) + 0.5 * (w.bn * hb.bn + w.be * hb.be + w.bu * hb.bu);
    r.p = Scale(vr, j.m);
    r.h = Add(s.Att().B2L(hb), Scale(CrossL(rr, vr), j.m));
    for (int t = 0; t < v.model.NumTanks(); ++t)
    {
        for (int n = 0; n < v.model.Params(t).nModes; ++n)
        {
            int idx = v.model.AuxIndex(t, n);
            double m = v.model.Params(t).mode[n].m;
            double k = v.model.Params(t).mode[n].k;
            CfBodyVector rho = Sub(v.model.RestArm(t, n), s.CmArm());
            rho.bn += aux.x[idx];
            rho.be += aux.x[idx + 1];
            CfBodyVector vrel { };
            vrel.bn = aux.x[idx + 2];
            vrel.be = aux.x[idx + 3];
            CfLocalVector pj = Add(rr, s.Att().B2L(rho));
            CfLocalVector vj = Add(vr, s.Att().B2L(Add(OmegaCrossR(w, rho), vrel)));
            double xi2 = aux.x[idx] * aux.x[idx] + aux.x[idx + 1] * aux.x[idx + 1];
            r.energy += 0.5 * m * Norm2(vj) + 0.5 * k * xi2;
            r.p = Add(r.p, Scale(vj, m));
            r.h = Add(r.h, Scale(CrossL(pj, vj), m));
        }
    }
    return r;
}

// Частота xi_bn первого тона бака 0 после начального смещения, свободный
// полёт с тягой изделия; полупериоды -- по линейно уточнённым нулям.
//   v -- изделие (параметры обновлены); fGuess -- ожидаемая частота, Гц;
//   xi0 -- начальное смещение, м; periods -- периодов для замера
// Возвращает частоту, Гц (0 -- мало нулей).
inline double SimFrequency(Vehicle* v, double fGuess, double xi0, int periods)
{
    AuxState aux { };
    v->model.ResetAux(&aux);
    int idx = v->model.AuxIndex(0, 0);
    aux.x[idx] = xi0;
    CfBodyOmega w0 { };
    CfLocalVector v0 { };
    VehState s = InitState(*v, w0, v0);
    double dt = 1.0 / (200.0 * fGuess);
    int nMax = periods * 200 * 4;
    int nCross = 0;
    int nNeed = 2 * periods + 1;
    double tFirst = 0.0;
    double tLast = 0.0;
    double prev = aux.x[idx];
    for (int i = 0; i < nMax && nCross < nNeed; ++i)
    {
        double t0 = s.Time();
        s = IntegrateRK4(s, &aux, dt, FnVeh, v);
        double cur = aux.x[idx];
        if ((prev < 0.0 && cur >= 0.0) || (prev > 0.0 && cur <= 0.0))
        {
            double tc = t0 + dt * prev / (prev - cur);
            if (nCross == 0)
            {
                tFirst = tc;
            }
            tLast = tc;
            ++nCross;
        }
        prev = cur;
    }
    if (nCross < 3)
    {
        return 0.0;
    }
    return (nCross - 1) / (2.0 * (tLast - tFirst));
}

} // namespace sltest

#endif // TEST_SLOSH_SIM_H

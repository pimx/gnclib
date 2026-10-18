// test_slosh_runs.h -- gnclib 1.7.0: прогоны проверок плескания
// (сохранение, частоты, параметры по заполнению, процесс во времени);
// CSV -- через запятую, одна строка заголовка.
#ifndef TEST_SLOSH_RUNS_H
#define TEST_SLOSH_RUNS_H

#include "test_slosh_sim.h"

namespace sltest
{

// Относительный уход инвариантов за прогон.
struct Drift
{
    double e = 0.0;
    double p = 0.0;
    double h = 0.0;
};

// Два бака (второй вне оси), по 3 тона, zeta = 0, без внешних сил, 20 с:
// уход энергии, импульса и момента импульса.
//   dt -- шаг, с; tumble -- 1: начальная угловая скорость (0.3, -0.2, 0.5)
//   рад/с, 0: покой (движение только от плескания)
inline Drift ConsRun(double dt, int tumble)
{
    SloshTank tk[2];
    tk[0] = TankExample();
    tk[0].radius = 1.4;
    tk[0].length = 6.0;
    tk[0].armAft = 9.0;
    tk[1] = TankExample();
    tk[1].radius = 0.8;
    tk[1].length = 5.0;
    tk[1].armAft = 2.0;
    tk[1].offBn = 0.6;
    tk[1].offBe = -0.4;
    tk[1].rho = 810.0;
    for (int n = 0; n < SLOSH_NMODE; ++n)
    {
        tk[0].zetaFix[n] = 1.0e-30;
        tk[1].zetaFix[n] = 1.0e-30;
    }
    Vehicle v;
    int nAux = v.model.Init(tk, 2, 0);
    v.dry.m = 3000.0;
    v.dry.j_bubu = 2500.0;
    v.dry.j_bebe = 50000.0;
    v.dry.j_bnbn = 52000.0;
    v.dry.j_bnbe = 300.0;
    v.dry.j_bnbu = -150.0;
    v.dry.j_bebu = 200.0;
    v.dryCm.bn = 0.05;
    v.dryCm.be = -0.03;
    v.dryCm.bu = 8.0;
    double mL[2];
    mL[0] = 0.6 * MassFull(tk[0]);
    mL[1] = 0.3 * MassFull(tk[1]);
    VehUpdate(&v, mL, 20.0);
    v.thrust = 0.0;
    AuxState aux { };
    v.model.ResetAux(&aux);
    for (int i = 0; i < nAux; ++i)
    {
        double sgn = (i % 2 == 0) ? 1.0 : -1.0;
        aux.x[i] = sgn * ((i % 4 < 2) ? 0.02 + 0.003 * i : 0.05 + 0.01 * i);
    }
    CfBodyOmega w0 { };
    if (tumble != 0)
    {
        w0.bn = 0.3;
        w0.be = -0.2;
        w0.bu = 0.5;
    }
    CfLocalVector vel { };
    vel.n = 1.0;
    vel.e = 2.0;
    vel.u = 3.0;
    VehState s = InitState(v, w0, vel);
    Invariants a = Invar(v, s, aux);
    int nStep = (int)(20.0 / dt + 0.5);
    for (int i = 0; i < nStep; ++i)
    {
        s = IntegrateRK4(s, &aux, dt, FnVeh, &v);
    }
    Invariants b = Invar(v, s, aux);
    Drift d;
    d.e = (b.energy - a.energy) / a.energy;
    d.p = Norm(Sub(b.p, a.p)) / Norm(a.p);
    d.h = Norm(Sub(b.h, a.h)) / Norm(a.h);
    return d;
}

// Сохранение: движение от плескания и кувыркание -- точность RK4 (с 1.7.1
// ориентация в IntegrateRK4 4-го порядка; до 1.7.0 уход при кувыркании
// убывал лишь как dt^2).
inline void RunConservation()
{
    SloshTank k = TankExample();
    SloshModel m;
    SloshTank two[2];
    two[0] = k;
    two[1] = k;
    Check("(8) Init of 2 tanks x 3 modes: aux states - 24", m.Init(two, 2, 0) - 24, 0.0);
    Drift r = ConsRun(0.002, 0);
    Check("(8) slosh-driven, dt 2 ms, 20 s: relative energy drift", r.e, 1.0e-10);
    Check("(8) slosh-driven: relative drift of the linear momentum", r.p, 1.0e-11);
    Check("(8) slosh-driven: relative drift of the angular momentum", r.h, 1.0e-10);
    Drift a = ConsRun(0.002, 1);
    Check("(8) tumbling, dt 2 ms, 20 s: relative energy drift", a.e, 1.0e-10);
    Check("(8) tumbling: relative drift of the linear momentum", a.p, 1.0e-11);
    Check("(8) tumbling: relative drift of the angular momentum", a.h, 1.0e-10);
}

// Изделие прогонов частот: пример бака, nModes тонов, демпфирование zeta,
// сухая часть с множителем scale (1e9 -- "зажатый" бак на вибростенде).
//   v -- выход; nModes -- тонов; zeta -- демпфирование; scale -- множитель
inline void MakeFreqVehicle(Vehicle* v, int nModes, double zeta, double scale)
{
    SloshTank k = TankExample();
    k.nModes = nModes;
    for (int n = 0; n < SLOSH_NMODE; ++n)
    {
        k.zetaFix[n] = zeta;
    }
    v->model.Init(&k, 1, 0);
    INR d { };
    d.m = 4000.0 * scale;
    d.j_bubu = 3000.0 * scale;
    d.j_bebe = 60000.0 * scale;
    d.j_bnbn = 60000.0 * scale;
    v->dry = d;
    CfBodyVector c { };
    c.bu = 9.0;
    v->dryCm = c;
}

// Проверка частот: бак отдельно (аналитика) против "зажатого" 6DOF,
// связанная теория w^2 (1 + m/M + m x (x + L)/J) против свободного 6DOF;
// slosh_freq.csv.
//   dir -- каталог вывода
inline void RunFrequency(const char* dir)
{
    char path[512];
    FILE* f = fopen(PathJoin(path, sizeof(path), dir, "slosh_freq.csv"), "w");
    if (f == 0)
    {
        Check("(9) open slosh_freq.csv", 1.0, 0.0);
        return;
    }
    fprintf(f, "nx,q,h_m,mL_kg,m1_kg,arm1_m,L1_m,cmRigid_m,f_tank_Hz,f_sim_clamped_Hz,"
               "err_clamped,f_coupled_Hz,f_sim_free_Hz,err_free\n");
    const double nxs[3] = { 1.0, 2.0, 4.0 };
    double errC = 0.0;
    double errF = 0.0;
    int nRows = 0;
    Vehicle v;
    Vehicle vc;
    MakeFreqVehicle(&v, 1, 1.0e-30, 1.0);
    MakeFreqVehicle(&vc, 1, 1.0e-30, 1.0e9);
    double mFull = MassFull(v.model.Tank(0));
    for (int i = 0; i < 3; ++i)
    {
        for (int iq = 1; iq <= 19; ++iq)
        {
            double q = 0.05 * iq;
            double mL = q * mFull;
            double a = nxs[i] * G0;
            VehUpdate(&v, &mL, a);
            VehUpdate(&vc, &mL, a);
            v.thrust = a * (v.rigid.m + v.model.MassSlosh());
            vc.thrust = a * (vc.rigid.m + vc.model.MassSlosh());
            const SloshModeParams& p = v.model.Params(0).mode[0];
            double fTank = p.w / (2.0 * PI);
            double x = p.arm - v.rigidCm.bu;
            double om2 = p.w * p.w * (1.0 + p.m / v.rigid.m
                                      + p.m * x * (x + p.lPend) / v.rigid.j_bebe);
            double fCoup = (om2 > 0.0) ? sqrt(om2) / (2.0 * PI) : 0.0;
            double fC = SimFrequency(&vc, fTank, 0.01, 20);
            double fF = (fCoup > 0.0) ? SimFrequency(&v, fCoup, 0.01, 20) : 0.0;
            double eC = fC / fTank - 1.0;
            double eF = (fCoup > 0.0) ? fF / fCoup - 1.0 : 0.0;
            errC = (fabs(eC) > errC) ? fabs(eC) : errC;
            errF = (fabs(eF) > errF) ? fabs(eF) : errF;
            double h = v.model.Params(0).armSurf - v.model.Tank(0).armAft;
            fprintf(f, "%g,%.2f,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.9g,%.9g,%.3e,%.9g,%.9g,%.3e\n",
                    nxs[i], q, h, mL, p.m, p.arm, p.lPend, v.rigidCm.bu, fTank, fC, eC, fCoup, fF,
                    eF);
            ++nRows;
        }
    }
    fclose(f);
    printf("    slosh_freq.csv: %d cases (nx = 1, 2, 4; q = 0.05 .. 0.95)\n", nRows);
    Check("(9) clamped 6DOF vs tank frequency: max relative error", errC, 1.0e-5);
    Check("(9) free 6DOF vs coupled theory: max relative error", errF, 1.0e-5);
}

// Параметры по заполнению при nx = 1, 2, 4 через SloshModel::Update;
// slosh_params.csv.
//   dir -- каталог вывода
inline void RunParamsSweep(const char* dir)
{
    char path[512];
    FILE* f = fopen(PathJoin(path, sizeof(path), dir, "slosh_params.csv"), "w");
    if (f == 0)
    {
        Check("(10) open slosh_params.csv", 1.0, 0.0);
        return;
    }
    fprintf(f, "nx,q,h_m,mL_kg,f1_Hz,f2_Hz,f3_Hz,m1_kg,m2_kg,m3_kg,arm1_m,arm2_m,arm3_m,"
               "L1_m,m0_kg,arm0_m,I0_kgm2,zeta1,eta1,bond,valid\n");
    Vehicle v;
    MakeFreqVehicle(&v, 3, 0.0, 1.0);
    double mFull = MassFull(v.model.Tank(0));
    const double nxs[3] = { 1.0, 2.0, 4.0 };
    for (int i = 0; i < 3; ++i)
    {
        for (int iq = 1; iq <= 49; ++iq)
        {
            double q = 0.02 * iq;
            double mL = q * mFull;
            VehUpdate(&v, &mL, nxs[i] * G0);
            const SloshTankParams& p = v.model.Params(0);
            double h = p.armSurf - v.model.Tank(0).armAft;
            fprintf(f, "%g,%.2f,%.6g,%.6g", nxs[i], q, h, mL);
            for (int n = 0; n < 3; ++n)
            {
                fprintf(f, ",%.6g", p.mode[n].w / (2.0 * PI));
            }
            for (int n = 0; n < 3; ++n)
            {
                fprintf(f, ",%.6g", p.mode[n].m);
            }
            for (int n = 0; n < 3; ++n)
            {
                fprintf(f, ",%.6g", p.mode[n].arm);
            }
            fprintf(f, ",%.6g,%.6g,%.6g,%.6g,%.4e,%.6g,%.4e,%d\n", p.mode[0].lPend, p.m0, p.arm0,
                    p.i0, p.mode[0].zeta0, p.mode[0].eta, p.bond, p.valid);
        }
    }
    fclose(f);
    printf("    slosh_params.csv: 147 rows (nx = 1, 2, 4; q = 0.02 .. 0.98)\n");
}

// Свободный полёт с тягой (nx = 2, q = 0.5), 3 тона, гладкая стенка и две
// кольцевые перегородки (одна затоплена), начальные смещение и поперечная
// скорость; slosh_time.csv. Проверка затухания первого тона.
//   dir -- каталог вывода
inline void RunTimeHistory(const char* dir)
{
    SloshTank k = TankExample();
    k.nBaffle = 2;
    k.baffle[0].arm = 9.5;
    k.baffle[0].width = 0.15;
    k.baffle[1].arm = 11.5;
    k.baffle[1].width = 0.15;
    Vehicle v;
    v.model.Init(&k, 1, 0);
    v.dry.m = 4000.0;
    v.dry.j_bubu = 3000.0;
    v.dry.j_bebe = 60000.0;
    v.dry.j_bnbn = 60000.0;
    v.dryCm.bu = 9.0;
    double mL = 0.5 * MassFull(k);
    double a = 2.0 * G0;
    VehUpdate(&v, &mL, a);
    v.thrust = a * (v.rigid.m + v.model.MassSlosh());
    AuxState aux { };
    v.model.ResetAux(&aux);
    int i1 = v.model.AuxIndex(0, 0);
    aux.x[i1] = 0.05;
    aux.x[i1 + 3] = 0.03;
    CfBodyOmega w0 { };
    CfLocalVector v0 { };
    VehState s = InitState(v, w0, v0);
    char path[512];
    FILE* f = fopen(PathJoin(path, sizeof(path), dir, "slosh_time.csv"), "w");
    if (f == 0)
    {
        Check("(11) open slosh_time.csv", 1.0, 0.0);
        return;
    }
    fprintf(f, "t_s,xi1_bn_m,xi1_be_m,amp1_m,wave1_m,zeta1,xi2_bn_m,xi3_bn_m,"
               "w_bn_rads,w_be_rads,w_bu_rads,cmSys_bn_m,cmSys_be_m\n");
    SloshModeState st0;
    v.model.ModeState(aux, 0, 0, &st0);
    SloshModeState st;
    const double dt = 0.005;
    for (int i = 0; i <= 12000; ++i)
    {
        if (i % 5 == 0)
        {
            v.model.ModeState(aux, 0, 0, &st);
            CfBodyVector cs = v.model.SystemCm(aux, v.rigid.m, v.rigidCm);
            fprintf(f, "%.3f,%.6e,%.6e,%.6e,%.6e,%.5f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                    s.Time(), st.xiBn, st.xiBe, st.amp, st.wave, st.zeta,
                    aux.x[v.model.AuxIndex(0, 1)], aux.x[v.model.AuxIndex(0, 2)], s.OmgBdw().bn,
                    s.OmgBdw().be, s.OmgBdw().bu, cs.bn, cs.be);
        }
        s = IntegrateRK4(s, &aux, dt, FnVeh, &v);
    }
    fclose(f);
    printf("    slosh_time.csv: 60 s, zeta1 %.4f at the start, %.4f at the end, amplitude "
           "%.4f -> %.4f m\n", st0.zeta, st.zeta, st0.amp, st.amp);
    Check("(11) baffle damping: amplitude ratio after 60 s (< 0.5)", st.amp / st0.amp, 0.5);
}

} // namespace sltest

#endif // TEST_SLOSH_RUNS_H

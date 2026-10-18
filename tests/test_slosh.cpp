// test_slosh.cpp -- gnclib 1.7.0: модель плескания топлива (lib_slosh.h).
// Сборка: g++ -std=c++17 -O2 -Iinclude tests/test_slosh.cpp build/libgnc.a
//         -o build/test_slosh
// Запуск: build/test_slosh [каталог для CSV, по умолчанию .]
#include "test_slosh_runs.h"

using namespace gnc;
using namespace sltest;

namespace
{

// Наибольшее относительное отклонение строки от эталона (3 тона).
double RowDev(const SloshRow& r, const double* mr, const double* z, const double* w2a,
              const double* eta, double m0r, double z0, double i0r)
{
    double d = 0.0;
    double v[4];
    for (int n = 0; n < 3; ++n)
    {
        v[0] = fabs(r.mr[n] / mr[n] - 1.0);
        v[1] = fabs(r.z[n] / z[n] - 1.0);
        v[2] = fabs(r.w2a[n] / w2a[n] - 1.0);
        v[3] = fabs(r.eta[n] / eta[n] - 1.0);
        for (int k = 0; k < 4; ++k)
        {
            d = (v[k] > d) ? v[k] : d;
        }
    }
    v[0] = fabs(r.m0r / m0r - 1.0);
    v[1] = fabs(r.z0 / z0 - 1.0);
    v[2] = fabs(r.i0r / i0r - 1.0);
    for (int k = 0; k < 3; ++k)
    {
        d = (v[k] > d) ? v[k] : d;
    }
    return d;
}

// (1) строки против эталона для потенциального течения (Python/scipy,
// 2000 членов; импульсная инерция сверена с МКР); (2) масса и ЦМ жидкости
// сохраняются.
void TestReference()
{
    const double mrA[3] = { 0.432196686801, 0.0136782221906, 0.00326001567547 };
    const double zA[3] = { 0.211184337617, 0.628478654164, 0.765798826876 };
    const double wA[3] = { 1.75079757453, 5.33119329554, 8.53631570907 };
    const double eA[3] = { 1.46512849354, 0.388793645451, 0.237553094818 };
    const double mrB[3] = { 0.775499017646, 0.0456509863944, 0.0119700747107 };
    const double zB[3] = { 0.0078463621865, 0.0560695169334, 0.113998830216 };
    const double wB[3] = { 0.558487205522, 3.16316913392, 5.5721797119 };
    const double eB[3] = { 0.701042367771, 0.346025736804, 0.232598332915 };
    const double mrC[3] = { 0.113627198195, 0.00341971556882, 0.000815003981621 };
    const double zC[3] = { 4.37267562541, 5.43730053493, 5.64856035423 };
    const double wC[3] = { 1.22745486997, 3.55429518235, 5.69087757756 };
    const double eC[3] = { 1.54076558938, 0.3888118394, 0.237553113109 };
    SloshRow r;
    SloshCylinderRow(1.0, 1.0, 3, &r);
    Check("(1) R = 1, h = 1: max rel. deviation from the reference",
          RowDev(r, mrA, zA, wA, eA, 0.550865075333, 0.721835273154, 0.052898982162), 1.0e-6);
    double dm = r.m0r + r.mr[0] + r.mr[1] + r.mr[2] - 1.0;
    double dz = r.m0r * r.z0 + r.mr[0] * r.z[0] + r.mr[1] * r.z[1] + r.mr[2] * r.z[2] - r.zLiq;
    Check("(2) R = 1, h = 1: sum of masses - liquid mass (fraction)", dm, 1.0e-14);
    Check("(2) R = 1, h = 1: CoM of the analog - CoM of the liquid, m", dz, 1.0e-14);
    SloshCylinderRow(1.5, 0.4, 3, &r);
    Check("(1) R = 1.5, h = 0.4: max rel. deviation from the reference",
          RowDev(r, mrB, zB, wB, eB, 0.166879921249, 1.13848897543, 0.350516000885), 1.0e-6);
    SloshCylinderRow(1.5, 6.0, 3, &r);
    Check("(1) R = 1.5, h = 6: max rel. deviation from the reference",
          RowDev(r, mrC, zC, wC, eC, 0.882138082255, 2.8112917348, 1.64595869141), 1.0e-6);
}

// (3) присоединение тонов сверх nModes; (4) зеркальные положения при a < 0;
// (5) пределы применимости.
void TestParams()
{
    SloshTank k = TankExample();
    k.nModes = 1;
    double mL = 0.4 * MassFull(k);
    double h = 0.4 * k.length;
    SloshRow r3;
    SloshRow r1;
    SloshCylinderRow(k.radius, h, 3, &r3);
    SloshCylinderRow(k.radius, h, 1, &r1);
    SloshTankParams pa;
    SloshTankParams pb;
    SloshRowToParams(k, r3, mL, 20.0, 0.1, 100.0, &pa);
    SloshRowToParams(k, r1, mL, 20.0, 0.1, 100.0, &pb);
    double d = fabs(pa.m0 / pb.m0 - 1.0) + fabs(pa.arm0 - pb.arm0) + fabs(pa.i0 / pb.i0 - 1.0)
               + fabs(pa.mode[0].m / pb.mode[0].m - 1.0) + fabs(pa.mode[0].w / pb.mode[0].w - 1.0);
    Check("(3) 3-mode row lumped to 1 mode - direct 1-mode row", d, 1.0e-12);
    k.nModes = 3;
    SloshTankParams pu;
    SloshTankParams pInv;
    SloshRowToParams(k, r3, mL, 20.0, 0.1, 100.0, &pu);
    SloshRowToParams(k, r3, mL, -20.0, 0.1, 100.0, &pInv);
    double ends = 2.0 * k.armAft + k.length;
    double dm = fabs(pu.arm0 + pInv.arm0 - ends) + fabs(pu.armSurf + pInv.armSurf - ends);
    for (int n = 0; n < 3; ++n)
    {
        dm += fabs(pu.mode[n].arm + pInv.mode[n].arm - ends);
        dm += fabs(pu.mode[n].w - pInv.mode[n].w);
    }
    Check("(4) a < 0: positions mirrored about the tank middle, m", dm, 1.0e-12);
    SloshTankParams pv;
    SloshRowToParams(k, r3, mL, 0.05, 0.1, 100.0, &pv);
    Check("(5) |a| = 0.05 < aMin: valid flag (expected 0)", pv.valid, 0.0);
    Check("(5) |a| = 0.05 < aMin: |a| clamped to aMin - 0.1", pv.a - 0.1, 1.0e-15);
    printf("    Bond number of the example tank at a = 20 m/s^2: %.3e (valid %d)\n", pu.bond,
           pu.valid);
    k.sigma = 1.0e5;
    SloshRowToParams(k, r3, mL, 20.0, 0.1, 100.0, &pv);
    Check("(5) Bond number below bondMin: valid flag (expected 0)", pv.valid, 0.0);
    printf("    smooth-wall zeta of the example tank, q = 0.4, a = 20 m/s^2: %.3e\n",
           pu.mode[0].zeta0);
    printf("    Miles zeta: ring w = 0.15 m at d = 0.5 m, R = 1.5 m, wave 0.05 m: %.4f\n",
           SloshZetaBaffle(0.5, 0.15, 1.5, 0.05));
}

// (6) интерполяция таблицы против прямого расчёта в серединах интервалов;
// (7) запись и чтение CSV, отказ на испорченном файле.
//   dir -- каталог вывода
void TestTable(const char* dir)
{
    SloshTank k = TankExample();
    static SloshTable tab;
    static SloshTable back;
    SloshTableCylinder(k.radius, k.length, 3, 48, &tab);
    double ew = 0.0;
    double em = 0.0;
    double ei = 0.0;
    double ez = 0.0;
    for (int i = 0; i + 1 < tab.n; ++i)
    {
        double q = 0.5 * (tab.q[i] + tab.q[i + 1]);
        SloshRow a;
        SloshRow b;
        SloshTableRow(tab, q, &a);
        SloshCylinderRow(k.radius, q * k.length, 3, &b);
        double vw = fabs(sqrt(a.w2a[0] / b.w2a[0]) - 1.0);
        double vm = fabs(a.mr[0] / b.mr[0] - 1.0);
        double vi = fabs(a.i0r / b.i0r - 1.0);
        double vz = fabs(a.z[0] - b.z[0]);
        ew = (vw > ew) ? vw : ew;
        em = (vm > em) ? vm : em;
        ei = (vi > ei) ? vi : ei;
        ez = (vz > ez) ? vz : ez;
    }
    Check("(6) 48-row table, mid points: max rel. error of f1", ew, 1.0e-3);
    Check("(6) 48-row table, mid points: max rel. error of m1", em, 5.0e-3);
    Check("(6) 48-row table, mid points: max rel. error of I0", ei, 1.0e-2);
    Check("(6) 48-row table, mid points: max error of the arm of m1, m", ez, 5.0e-3);
    char path[512];
    char err[256];
    PathJoin(path, sizeof(path), dir, "slosh_table.csv");
    bool okSave = SloshTableSaveCsv(path, tab, err, sizeof(err));
    bool okLoad = SloshTableLoadCsv(path, &back, err, sizeof(err));
    double d = (okSave && okLoad && back.n == tab.n) ? 0.0 : 1.0;
    for (int i = 0; i < back.n && d < 1.0; ++i)
    {
        const SloshRow& a = tab.row[i];
        const SloshRow& b = back.row[i];
        double v = fabs(tab.q[i] - back.q[i]) + fabs(a.zSurf - b.zSurf) + fabs(a.i0r / b.i0r - 1.0)
                   + fabs(a.z0 - b.z0) + fabs(a.m0r - b.m0r);
        for (int m = 0; m < 3; ++m)
        {
            v += fabs(a.mr[m] / b.mr[m] - 1.0) + fabs(a.w2a[m] / b.w2a[m] - 1.0);
            v += fabs(a.z[m] - b.z[m]) + fabs(a.eta[m] - b.eta[m]);
        }
        d = (v > d) ? v : d;
    }
    printf("    slosh_table.csv: %d rows\n", back.n);
    Check("(7) table CSV save/load round trip: max deviation", d, 1.0e-8);
    PathJoin(path, sizeof(path), dir, "slosh_table_bad.csv");
    FILE* f = fopen(path, "w");
    if (f != 0)
    {
        fprintf(f, "q,zSurf\n0.1,0.5,0.25\n");
        fclose(f);
    }
    SloshTable keep = back;
    bool okBad = SloshTableLoadCsv(path, &back, err, sizeof(err));
    printf("    bad table refused: %s\n", err);
    double changed = (back.n != keep.n || back.q[0] != keep.q[0]) ? 1.0 : 0.0;
    Check("(7) short row refused (1 = accepted)", okBad ? 1.0 : 0.0, 0.0);
    Check("(7) output kept on refusal (1 = changed)", changed, 0.0);
}

} // namespace

int main(int argc, char** argv)
{
    const char* dir = (argc > 1) ? argv[1] : ".";
    TestReference();
    TestParams();
    TestTable(dir);
    RunConservation();
    RunFrequency(dir);
    RunParamsSweep(dir);
    RunTimeHistory(dir);
    printf(BadCount() == 0 ? "TEST_SLOSH PASS\n" : "TEST_SLOSH FAIL\n");
    return BadCount();
}

// test_point.cpp -- gnclib 1.4.0 / 1.5.0: the cheap paths of a model point
// give the numbers of the full ones -- EcpToAlt / GEO::AltOf and (1.5.0)
// EcpToGeoSc, GeoRadius against EcpToGeo and GeoToEcp,
// EnvAtmosphere::State against Temperature, Pressure and Density, the
// strip loop of CylinderStrip against the general products of 1.3.0, and
// CylinderStripPair against two CylinderStrip calls.
// Build: g++ -std=c++17 -O2 -Iinclude tests/test_point.cpp -o build/test_point
#include "lib.h"
#include "lib_random.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

using namespace gnc;

static int g_bad = 0;

// Check of a condition with a print.
static void Check(const char* name, bool ok)
{
    printf("%-64s %s\n", name, ok ? "OK" : "FAIL");
    g_bad += ok ? 0 : 1;
}

// Two doubles equal to the bit.
static bool Same(double a, double b)
{
    return memcmp(&a, &b, sizeof(double)) == 0;
}

// The strip loop of 1.3.0 with the general products: the reference.
static void StripRef(const CfBodyVector& vRel, const CfBodyOmega& omg, double rho, double length,
                     double radius, double cmFromTail, double cdLateral, int nStrips,
                     CfBodyVector* force, CfBodyMoment* moment)
{
    CfBodyVector f0 { };
    CfBodyMoment m0 { };
    *force = f0;
    *moment = m0;
    double dx = length / (double)nStrips;
    double dia = 2.0 * radius;
    int i = 0;
    for (i = 0; i < nStrips; i = i + 1)
    {
        double xi = ((double)i + 0.5) * dx - cmFromTail;
        CfBodyVector r { };
        r.bu = xi;
        CfBodyVector vloc = Add(vRel, OmegaCrossR(omg, r));
        CfBodyVector vlat { };
        vlat.be = vloc.be;
        vlat.bn = vloc.bn;
        double vm = Norm(vlat);
        if (vm < 1.0e-9)
        {
            continue;
        }
        CfBodyVector df = Scale(vlat, -0.5 * rho * cdLateral * dia * dx * vm);
        *force = Add(*force, df);
        *moment = Add(*moment, Torque(r, df));
    }
}

// Force and moment equal to the bit, the zero components by value.
static bool SameFm(const CfBodyVector& f, const CfBodyMoment& m, const CfBodyVector& g, const CfBodyMoment& n)
{
    return Same(f.be, g.be) && Same(f.bn, g.bn) && f.bu == g.bu && Same(m.be, n.be) && Same(m.bn, n.bn) &&
           m.bu == n.bu;
}

// Height: EcpToAlt against EcpToGeo on a grid.
static void TestHeight()
{
    CfGeoid el = ELL_PZ9011;
    double worst = 0.0;
    int i = 0;
    for (i = 0; i <= 900; i = i + 1)
    {
        int j = 0;
        for (j = 0; j <= 400; j = j + 1)
        {
            CfGeodetic g { };
            g.lat = (double)i * 0.1 * PI / 180.0;
            g.lon = 0.7;
            g.alt = -2000.0 + (double)j * 1000.0;
            CfEcefPoint p = GeoToEcp(g, el);
            double d = fabs(EcpToAlt(p, el) - EcpToGeo(p, el).alt);
            worst = d > worst ? d : worst;
        }
    }
    printf("  EcpToAlt - EcpToGeo().alt, worst %.3e m\n", worst);
    Check("EcpToAlt: the height of EcpToGeo within 2e-8 m", worst < 2.0e-8);
    CfEcefPoint pole { };
    pole.z = el.a * (1.0 - el.f) + 500.0;
    Check("EcpToAlt: the pole as EcpToGeo", EcpToAlt(pole, el) == EcpToGeo(pole, el).alt);
    GEO geo;
    CfGeodetic org { };
    org.lat = 62.925 * PI / 180.0;
    org.lon = 40.578 * PI / 180.0;
    org.alt = 100.0;
    geo.SetOrigin(org);
    CfLocalPoint l { };
    l.n = 30000.0;
    l.e = 60000.0;
    l.u = 130000.0;
    Check("GEO::AltOf: the height of GeoOf within 2e-8 m", fabs(geo.AltOf(l) - geo.GeoOf(l).alt) < 2.0e-8);
}

// Latitude and radius (1.5.0): EcpToGeoSc against EcpToGeo, GeoRadius
// against the length of GeoToEcp, on a grid of both hemispheres.
static void TestLatitude()
{
    CfGeoid el = ELL_PZ9011;
    double worstSc = 0.0;
    double worstR = 0.0;
    bool altSame = true;
    int i = 0;
    for (i = -900; i <= 900; i = i + 1)
    {
        int j = 0;
        for (j = 0; j <= 80; j = j + 1)
        {
            CfGeodetic g { };
            g.lat = (double)i * 0.1 * PI / 180.0;
            g.lon = 0.3 + (double)i * 0.001;
            g.alt = -2000.0 + (double)j * 5000.0;
            CfEcefPoint p = GeoToEcp(g, el);
            CfGeoSc sc = EcpToGeoSc(p, el);
            CfGeodetic ref = EcpToGeo(p, el);
            double ds = fabs(sc.sinLat - sin(ref.lat));
            double dc = fabs(sc.cosLat - cos(ref.lat));
            worstSc = ds > worstSc ? ds : worstSc;
            worstSc = dc > worstSc ? dc : worstSc;
            altSame = altSame && Same(sc.alt, EcpToAlt(p, el));
            CfGeodetic t = ref;
            t.alt = 300000.0 + (double)j * 100.0;
            CfEcefPoint e = GeoToEcp(t, el);
            double r = sqrt(e.x * e.x + e.y * e.y + e.z * e.z);
            double dr = fabs(GeoRadius(sc.sinLat, sc.cosLat, t.alt, el) - r);
            worstR = dr > worstR ? dr : worstR;
        }
    }
    printf("  EcpToGeoSc - sin, cos of EcpToGeo().lat, worst %.3e; GeoRadius - |GeoToEcp|, worst %.3e m\n",
           worstSc, worstR);
    Check("EcpToGeoSc: the latitude of EcpToGeo within 1e-14", worstSc < 1.0e-14);
    Check("EcpToGeoSc: the height is EcpToAlt to the bit", altSame);
    Check("GeoRadius: the length of GeoToEcp within 2e-8 m", worstR < 2.0e-8);
    CfEcefPoint pole { };
    pole.z = -(el.a * (1.0 - el.f) + 700.0);
    CfGeoSc sp = EcpToGeoSc(pole, el);
    Check("EcpToGeoSc: the south pole -- sine -1, cosine 0, height of EcpToGeo",
          sp.sinLat == -1.0 && sp.cosLat == 0.0 && sp.alt == EcpToGeo(pole, el).alt);
}

// Atmosphere: State against the three calls, every 0.1 m, three offsets.
static void TestAtmosphere()
{
    double dts[3] = { 0.0, 15.0, -20.0 };
    long bad = 0;
    long n = 0;
    int k = 0;
    for (k = 0; k < 3; k = k + 1)
    {
        EnvAtmosphere a(dts[k]);
        long i = 0;
        for (i = -2000; i <= 1600000; i = i + 1)
        {
            double h = (double)i * 0.1 + (double)(i % 7) * 1.0e-5;
            AtmState s = a.State(h);
            bool ok = Same(s.t, a.Temperature(h)) && Same(s.p, a.Pressure(h)) && Same(s.rho, a.Density(h));
            bad = bad + (ok ? 0 : 1);
            n = n + 1;
        }
    }
    printf("  State against Temperature, Pressure, Density: %ld heights, %ld differ\n", n, bad);
    Check("EnvAtmosphere::State: the three calls to the bit", bad == 0);
}

// Strips: the 1.4.0 loop and the pair against the reference.
static void TestStrips()
{
    EnvAeroforces aero(ATM_STD);
    RNG rng(12345);
    long bad = 0;
    long badPair = 0;
    int n = 0;
    for (n = 0; n < 20000; n = n + 1)
    {
        CfBodyVector v { };
        v.bu = rng.UniformRange(-1500.0, 1500.0);
        v.be = rng.UniformRange(-300.0, 300.0);
        v.bn = rng.UniformRange(-300.0, 300.0);
        CfBodyOmega w { };
        w.bu = rng.UniformRange(-0.5, 0.5);
        w.be = rng.UniformRange(-0.5, 0.5);
        w.bn = rng.UniformRange(-0.5, 0.5);
        if (n % 5 == 0)
        {
            v.be = 0.0;
            v.bn = 0.0;
        }
        double rho = rng.UniformRange(1.0e-6, 1.3);
        double len = rng.UniformRange(10.0, 60.0);
        double cm = rng.UniformRange(5.0, 30.0);
        int strips = 1 + n % 32;
        CfBodyVector f { };
        CfBodyMoment m { };
        aero.CylinderStrip(v, w, rho, 2.0, 0.5, len, 1.85, cm, 1.2, 0.0, 0.0, strips, &f, &m);
        CfBodyVector fr { };
        CfBodyMoment mr { };
        StripRef(v, w, rho, len, 1.85, cm, 1.2, strips, &fr, &mr);
        bad = bad + (SameFm(f, m, fr, mr) ? 0 : 1);
        CfBodyVector f1 { };
        CfBodyMoment m1 { };
        CfBodyVector f0 { };
        CfBodyMoment m0 { };
        aero.CylinderStripPair(v, w, rho, len, 1.85, cm, 1.2, strips, &f1, &m1, &f0, &m0);
        CfBodyOmega still { };
        CfBodyVector fs { };
        CfBodyMoment ms { };
        aero.CylinderStrip(v, still, rho, 2.0, 0.5, len, 1.85, cm, 1.2, 0.0, 0.0, strips, &fs, &ms);
        badPair = badPair + (SameFm(f1, m1, f, m) && SameFm(f0, m0, fs, ms) ? 0 : 1);
    }
    printf("  20000 random bodies: strip loop %ld differ, pair %ld differ\n", bad, badPair);
    Check("CylinderStrip: the general products of 1.3.0 to the bit", bad == 0);
    Check("CylinderStripPair: two CylinderStrip calls to the bit", badPair == 0);
}

// Entry point.
int main()
{
    printf("TEST_POINT (gnclib %s)\n", GNCLIB_VERSION);
    TestHeight();
    TestLatitude();
    TestAtmosphere();
    TestStrips();
    printf("TEST_POINT %s\n", g_bad == 0 ? "PASS" : "FAIL");
    return g_bad == 0 ? 0 : 1;
}

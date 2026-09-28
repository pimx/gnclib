// gnclib 1.1.0: navigation ESKF (moved from FLY90 / HOP90 V90.43, with
// the Earth rate, IMU-lever noise and gate recovery of HOP90).
//============================================================================
// eskf_predict.cpp -- ESKF: настройка, инициализация, такт
// комплексирования (компенсация смещений, рычаг ИНС, механизация БИНС,
// адаптивный процессный шум), распространение ковариации.
//============================================================================
#include "lib_eskf.h"
#include "lib_earth.h"

namespace gnc
{

// Конструктор: фильтр не настроен.
Eskf::Eskf()
    : m_frame(0), m_g(gnc::G0_STD), m_qScale(1.0), m_gnssAge(0.0),
      m_time(0.0), m_wLast { }, m_lever { }, m_imuLever { }, m_wPrevLever { },
      m_leverInit(false), m_v { }, m_p { }, m_ba { }, m_bg { }, m_wDotF { }, m_wDotVar(0.0), m_wDotV { }
{
    int ir = 0;
    for (ir = 0; ir < N; ++ir)
    {
        m_rejectRun[ir] = 0;
    }
    int i = 0;
    int j = 0;
    for (i = 0; i < N; ++i)
    {
        for (j = 0; j < N; ++j)
        {
            m_P[i][j] = 0.0;
        }
    }
}

// Настройка: шумы, площадка, тяжесть, рычаг высотомера.
void Eskf::setup(const Config& c, const gnc::GEO* frame, double gravity,
                 gnc::CfBodyVector leverBody)
{
    m_cfg = c;
    m_frame = frame;
    m_g = gravity;
    m_lever = leverBody;
    m_imuLever = mkBody(0.0, 0.0, 0.0); // задаётся setImuLever
    m_wPrevLever = mkBody(0.0, 0.0, 0.0);
    m_leverInit = false;
}

// Инициализация точно известным начальным состоянием; начальная
// ковариация -- малые диагональные блоки.
void Eskf::init(gnc::CfLocalVector p0, gnc::CfLocalVector v0,
                const gnc::CfRotation& att0, Time t0)
{
    m_time = t0;
    m_p = p0;
    m_v = v0;
    m_R = mat3OfRot(att0);
    m_ba = mkBody(0.0, 0.0, 0.0);
    m_bg = mkBody(0.0, 0.0, 0.0);
    m_wLast = mkBody(0.0, 0.0, 0.0);
    m_qScale = 1.0;
    m_gnssAge = 0.0;
    int i = 0;
    int j = 0;
    for (i = 0; i < N; ++i)
    {
        for (j = 0; j < N; ++j)
        {
            m_P[i][j] = 0.0;
        }
    }
    setDiag(0, 0.005 * 0.005); // dp
    setDiag(3, 0.005 * 0.005); // dv
    setDiag(6, 0.002 * 0.002); // dTheta
    setDiag(9, 0.01 * 0.01);   // dba
    setDiag(12, 2e-4 * 2e-4);  // dbg
}

// Установка диагонального блока 3x3.
void Eskf::setDiag(int i0, double v)
{
    m_P[i0][i0] = v;
    m_P[i0 + 1][i0 + 1] = v;
    m_P[i0 + 2][i0 + 2] = v;
}

// Заполнение навигационного решения.
void Eskf::fillOut(NavData& out) const
{
    out.posLop = lopOf(m_p);
    out.velLoc = m_v;
    out.att = rotOfMat3(m_R);
    out.omgBdw = omgOf(m_wLast);
    out.gnssAge = m_gnssAge;
}

// Рычаг ИНС: omgDot -- разностью гироскопа по такту (шум давится тем,
// что добавка мала при малых dt); центростремительный член -- точно.
void Eskf::leverCompensate(const gnc::CfBodyVector& w, gnc::CfBodyVector* a,
                           double dt)
{
    if (norm(m_imuLever) <= 1e-9)
    {
        return;
    }
    gnc::CfBodyVector wDot = mkBody(0.0, 0.0, 0.0);
    if (m_leverInit && dt > 1e-9)
    {
        // V90.35: производная гироскопа сглаживается (20 мс) -- импульсы
        // РСУ и изгиб корпуса иначе дают ложные ускорения на рычаге ИНС
        gnc::CfBodyVector raw = scale(sub(w, m_wPrevLever), 1.0 / dt);
        double k = (m_cfg.wDotTau > 1e-9) ? 1.0 - exp(-dt / m_cfg.wDotTau) : 1.0;
        gnc::CfBodyVector dw = sub(raw, m_wDotF);
        m_wDotF = add(m_wDotF, scale(dw, k));
        wDot = m_wDotF;
        // the lever-noise variance always from the 20 ms smoothing (independent
        // of the compensation smoothing wDotTau; identical when it is 0.02)
        double kv = 1.0 - exp(-dt / 0.02);
        gnc::CfBodyVector dv = sub(raw, m_wDotV);
        m_wDotV = add(m_wDotV, scale(dv, kv));
        // V90.36 шум рычага: сглаживание не успевает за импульсами РСУ и
        // изгибом -- остаток (dw) на рычаге r даёт ошибку компенсации
        // ускорения ~ |dw| |r|; его дисперсия (0.2 с) идёт в шум процесса
        double ks = 1.0 - exp(-dt / 0.2);
        double d2 = dot(dv, dv) * (1.0 - kv) * (1.0 - kv);
        m_wDotVar += (d2 - m_wDotVar) * ks;
    }
    m_wPrevLever = w;
    m_leverInit = true;
    gnc::CfBodyVector aLev =
        add(cross(wDot, m_imuLever), cross(w, cross(w, m_imuLever)));
    *a = sub(*a, aLev);
}

// Один такт: механизация, ковариация, коррекции, выдача.
void Eskf::update(NavData& outNav, const ImuData& imu, const GnssData& gnss,
                  const AltimeterData& alt, const AltimeterData& ralt, Time t)
{
    double dt = (double)(t - m_time);
    if (dt < 0.0)
    {
        dt = 0.0;
    }
    m_time = t;
    if (dt <= 0.0)
    {
        fillOut(outNav);
        return;
    }
    // компенсация смещений и рычага ИНС
    gnc::CfBodyVector w = sub(omgVec(imu.gyro), m_bg);
    gnc::CfBodyVector a = sub(imu.accel, m_ba);
    leverCompensate(w, &a, dt);

    // механизация БИНС (V90.34: вращающаяся Земля -- поле тяжести в
    // оценке положения, кориолисово ускорение, угловая скорость Земли)
    Mat3 Reb = m_R;
    gnc::CfLocalVector gL = mkLoc(0.0, 0.0, -m_g);
    gnc::CfLocalVector omL = mkLoc(0.0, 0.0, 0.0);
    if (m_frame != 0)
    {
        gL = EarthGravityLoc(*m_frame, mkLop(m_p.e, m_p.n, m_p.u));
        omL = EarthRateLoc(*m_frame);
    }
    gnc::CfBodyVector omB = mkBody(
        Reb.m[0][0] * omL.e + Reb.m[1][0] * omL.n + Reb.m[2][0] * omL.u,
        Reb.m[0][1] * omL.e + Reb.m[1][1] * omL.n + Reb.m[2][1] * omL.u,
        Reb.m[0][2] * omL.e + Reb.m[1][2] * omL.n + Reb.m[2][2] * omL.u);
    w = sub(w, omB);
    gnc::CfLocalVector fLoc = mulv(Reb, a);
    gnc::CfLocalVector acc = sub(add(fLoc, gL), scale(cross(omL, m_v), 2.0));
    m_p = add(add(m_p, scale(m_v, dt)), scale(acc, 0.5 * dt * dt));
    m_v = add(m_v, scale(acc, dt));
    // кинематика на матрице: R <- R * expm([w dt]x) (Родриг)
    double wn = norm(w);
    if (wn > 1e-12)
    {
        m_R = m_R.mul(Mat3::fromAxisAngle(vec3Of(scale(w, 1.0 / wn)),
                                          wn * dt));
    }
    m_R = m_R.orthonormalize();

    // распространение ковариации с манёвренным масштабом Q
    double qInst = 1.0 + 4.0 * norm(acc) + 10.0 * norm(w);
    if (qInst > 30.0)
    {
        qInst = 30.0;
    }
    m_qScale *= exp(-dt / 2.0);
    if (qInst > m_qScale)
    {
        m_qScale = qInst;
    }
    propagate(Reb, a, w, dt, m_qScale);

    m_gnssAge += dt;
    gnssCorrect(gnss);
    altCorrect(alt, m_cfg.altNoiseStd, 0.0);
    altCorrect(ralt, m_cfg.raltNoiseStd, m_cfg.raltNoiseRel);

    m_wLast = w;
    fillOut(outNav);
}

// P = Fd P Fd^T + Qd, Fd = I + F*dt. Ненулевые блоки F:
//   F[dp][dv] = I; F[dv][dTheta] = -Reb*[a]x; F[dv][dba] = -Reb;
//   F[dTheta][dTheta] = -[w]x; F[dTheta][dbg] = -I.
void Eskf::propagate(const Mat3& Reb, const gnc::CfBodyVector& a,
                     const gnc::CfBodyVector& w, double dt, double qScale)
{
    double Fd[N][N];
    int i = 0;
    int j = 0;
    int k = 0;
    for (i = 0; i < N; ++i)
    {
        for (j = 0; j < N; ++j)
        {
            Fd[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    for (i = 0; i < 3; ++i)
    {
        Fd[i][3 + i] += dt; // dp <- dv
    }
    Mat3 RA = Reb.mul(Mat3::skew(vec3Of(a)));
    for (i = 0; i < 3; ++i)
    {
        for (j = 0; j < 3; ++j)
        {
            Fd[3 + i][6 + j] -= RA.m[i][j] * dt; // dv <- dTheta
        }
    }
    for (i = 0; i < 3; ++i)
    {
        for (j = 0; j < 3; ++j)
        {
            Fd[3 + i][9 + j] -= Reb.m[i][j] * dt; // dv <- dba
        }
    }
    Mat3 Wx = Mat3::skew(vec3Of(w));
    for (i = 0; i < 3; ++i)
    {
        for (j = 0; j < 3; ++j)
        {
            Fd[6 + i][6 + j] -= Wx.m[i][j] * dt; // dTheta <- dTheta
        }
        Fd[6 + i][12 + i] -= dt; // dTheta <- dbg
    }

    // T = Fd * P; P = T * Fd^T
    double T[N][N];
    for (i = 0; i < N; ++i)
    {
        for (j = 0; j < N; ++j)
        {
            double s = 0.0;
            for (k = 0; k < N; ++k)
            {
                s += Fd[i][k] * m_P[k][j];
            }
            T[i][j] = s;
        }
    }
    for (i = 0; i < N; ++i)
    {
        for (j = 0; j < N; ++j)
        {
            double s = 0.0;
            for (k = 0; k < N; ++k)
            {
                s += T[i][k] * Fd[j][k];
            }
            m_P[i][j] = s;
        }
    }

    // Qd (диагональ); dv и dTheta -- с манёвренным масштабом
    double qv = m_cfg.accelNoiseStd * dt;
    qv = qv * qv * qScale;
    // V90.36 шум рычага ИНС: sigma_a = sigma_wdot |r| (адаптивно)
    if (m_cfg.leverNoise != 0)
    {
        double sigL = m_cfg.leverNoiseK * sqrt(fmax(m_wDotVar, 0.0)) * norm(m_imuLever) * dt;
        qv += sigL * sigL;
    }
    double qt = m_cfg.gyroNoiseStd * dt;
    qt = qt * qt * qScale;
    double qa = m_cfg.accelBiasRw * m_cfg.accelBiasRw * dt;
    double qg = m_cfg.gyroBiasRw * m_cfg.gyroBiasRw * dt;
    for (i = 0; i < 3; ++i)
    {
        m_P[3 + i][3 + i] += qv;
        m_P[6 + i][6 + i] += qt;
        m_P[9 + i][9 + i] += qa;
        m_P[12 + i][12 + i] += qg;
    }
}

}   // namespace gnc

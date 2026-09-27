// gnclib 1.1.0: navigation ESKF (moved from FLY90 / HOP90 V90.43, with
// the Earth rate, IMU-lever noise and gate recovery of HOP90).
//============================================================================
// eskf_update.cpp -- ESKF: коррекции СНС и высотомера, скалярные
// коррекции (полная и ограниченная вертикальным каналом), симметризация.
//============================================================================
#include "lib_eskf.h"

namespace gnc
{

// Коррекция СНС: позиция и скорость ENU (из ECEF), 6 скаляров.
void Eskf::gnssCorrect(const GnssData& gnss)
{
    if (gnss.valid == 0 || m_frame == 0)
    {
        return;
    }
    m_gnssAge = 0.0;
    gnc::CfLocalVector pm = vecOf(m_frame->E2L(gnss.posEcp));
    gnc::CfLocalVector vm = m_frame->E2L(gnss.velEcf);
    double Rp = m_cfg.gnssPosStd * m_cfg.gnssPosStd;
    double Rv = m_cfg.gnssVelStd * m_cfg.gnssVelStd;
    scalarUpdate(0, 1.0, pm.e - m_p.e, Rp);
    scalarUpdate(1, 1.0, pm.n - m_p.n, Rp);
    scalarUpdate(2, 1.0, pm.u - m_p.u, Rp);
    scalarUpdate(3, 1.0, vm.e - m_v.e, Rv);
    scalarUpdate(4, 1.0, vm.n - m_v.n, Rv);
    scalarUpdate(5, 1.0, vm.u - m_v.u, Rv);
    symmetrize();
}

// Коррекция высотомера (1 скаляр, только вертикаль). Вблизи земли
// (дальность < 2.5 см по оценке опор) измерение не используется:
// показание ограничено нулём и шум несимметричен.
void Eskf::altCorrect(const AltimeterData& alt, double std0, double stdRel)
{
    if (alt.valid == 0)
    {
        return;
    }
    // луч высотомера в ENU (ось -bu связанной системы)
    gnc::CfLocalVector beam = mulv(m_R, mkBody(-1.0, 0.0, 0.0));
    double gearZ = add(m_p, mulv(m_R, m_lever)).u;
    if (beam.u < -0.2 && gearZ > 0.025)
    {
        double c = 1.0 / (-beam.u); // d(range)/d(p_u)
        double rangePred = gearZ * c;
        double sd = std0 + stdRel * rangePred;
        scalarUpdateMasked(2, c, alt.range - rangePred, sd * sd);
    }
}

// Скалярная коррекция: измерение z = c * x[idx] + шум(R); строб 6
// сигм; ошибка немедленно вводится в номинал (позиция/скорость -- ENU,
// ориентация и смещения -- BODY, инъекция R <- R * expm([dth]x)).
void Eskf::scalarUpdate(int idx, double c, double innov, double R)
{
    double S = c * c * m_P[idx][idx] + R;
    if (S < 1e-12)
    {
        return;
    }
    if (innov * innov > 36.0 * S)
    {
        // строб от грубых выбросов; V90.35: долгая серия отказов -- фильтр
        // разошёлся (P мала относительно фактической ошибки): дисперсия
        // измеряемой компоненты поднимается до квадрата невязки, и
        // измерение принимается -- иначе расхождение необратимо
        m_rejectRun++;
        if (m_rejectRun < 60)
        {
            return;
        }
        m_P[idx][idx] += innov * innov / (c * c + 1e-12);
        S = c * c * m_P[idx][idx] + R;
    }
    m_rejectRun = 0;
    double K[N];
    int i = 0;
    int j = 0;
    for (i = 0; i < N; ++i)
    {
        K[i] = m_P[i][idx] * c / S;
    }
    // P = (I - K h) P
    double row[N];
    for (j = 0; j < N; ++j)
    {
        row[j] = c * m_P[idx][j];
    }
    for (i = 0; i < N; ++i)
    {
        for (j = 0; j < N; ++j)
        {
            m_P[i][j] -= K[i] * row[j];
        }
    }
    gnc::CfLocalVector dp = mkLoc(K[0] * innov, K[1] * innov, K[2] * innov);
    gnc::CfLocalVector dv = mkLoc(K[3] * innov, K[4] * innov, K[5] * innov);
    Vec3 dth(K[6] * innov, K[7] * innov, K[8] * innov);
    gnc::CfBodyVector dba = mkBody(K[9] * innov, K[10] * innov,
                                   K[11] * innov);
    gnc::CfBodyVector dbg = mkBody(K[12] * innov, K[13] * innov,
                                   K[14] * innov);
    m_p = add(m_p, dp);
    m_v = add(m_v, dv);
    double an = norm(dth);
    if (an > 1e-15)
    {
        m_R = m_R.mul(Mat3::fromAxisAngle(scale(dth, 1.0 / an), an));
    }
    m_R = m_R.orthonormalize();
    m_ba = add(m_ba, dba);
    m_bg = add(m_bg, dbg);
}

// Скалярная коррекция, ограниченная вертикальным каналом
// {dp_u, dv_u, dba_bu}; ковариация -- форма Джозефа.
void Eskf::scalarUpdateMasked(int idx, double c, double innov, double R)
{
    double S = c * c * m_P[idx][idx] + R;
    if (S < 1e-12)
    {
        return;
    }
    if (innov * innov > 36.0 * S)
    {
        return; // строб 6 сигм
    }
    double K[N];
    int i = 0;
    int j = 0;
    for (i = 0; i < N; ++i)
    {
        K[i] = 0.0;
    }
    K[2] = m_P[2][idx] * c / S; // dp (вертикаль ENU)
    K[5] = m_P[5][idx] * c / S; // dv (вертикаль ENU)
    K[9] = m_P[9][idx] * c / S; // dba_bu (продольная BODY)
    // A P, где A = I - K h (h: столбец idx со значением c)
    double AP[N][N];
    for (i = 0; i < N; ++i)
    {
        for (j = 0; j < N; ++j)
        {
            AP[i][j] = m_P[i][j] - K[i] * c * m_P[idx][j];
        }
    }
    // (A P) A^T + K R K^T
    for (i = 0; i < N; ++i)
    {
        for (j = 0; j < N; ++j)
        {
            m_P[i][j] = AP[i][j] - AP[i][idx] * c * K[j] + K[i] * R * K[j];
        }
    }
    m_p.u += K[2] * innov;
    m_v.u += K[5] * innov;
    m_ba.bu += K[9] * innov;
}

// Симметризация ковариации после обновления.
void Eskf::symmetrize()
{
    int i = 0;
    int j = 0;
    for (i = 0; i < N; ++i)
    {
        for (j = i + 1; j < N; ++j)
        {
            double s = 0.5 * (m_P[i][j] + m_P[j][i]);
            m_P[i][j] = s;
            m_P[j][i] = s;
        }
    }
}

}   // namespace gnc

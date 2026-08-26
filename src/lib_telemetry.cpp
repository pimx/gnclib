#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

// lib_telemetry.cpp -- отправка UDP-телеметрии v2 через UdpSender.
#include "lib_telemetry.h"

UdpTelemetry::UdpTelemetry() : m_counter(0), m_maskAccum(0)
{
}

UdpTelemetry::~UdpTelemetry()
{
    close();
}

bool UdpTelemetry::setup(const Config& c)
{
    m_cfg = c;
    if (m_cfg.host == 0)
    {
        m_cfg.host = "127.0.0.1";
    }
    if (m_cfg.port <= 0)
    {
        m_cfg.port = 9000;
    }
    if (m_cfg.decimation < 1)
    {
        m_cfg.decimation = 1;
    }
    m_counter = 0;
    m_maskAccum = 0;
    if (!m_cfg.enable)
    {
        return true;
    }
    return m_sender.open();
}

void UdpTelemetry::close()
{
    m_sender.close();
}

// Сбор и отправка пакета v2 (см. lib_telemetry.h и TELEMETRY.md)
void UdpTelemetry::send(double t, VecEnu posEnu, Quat qBodyToEnu,
                        VecEnu thrustEnu, double mass, VecEnu velEnu, int mode,
                        int rcsMask)
{
    if (!m_sender.isOpen())
    {
        return;
    }
    m_maskAccum |= rcsMask; // защёлка активности на интервале
    if (m_counter++ % m_cfg.decimation != 0)
    {
        return;
    }

    double pkt[PKT_DOUBLES];
    pkt[0] = t;
    pkt[1] = posEnu.e;
    pkt[2] = posEnu.n;
    pkt[3] = posEnu.u;
    pkt[4] = qBodyToEnu.w;
    pkt[5] = qBodyToEnu.x;
    pkt[6] = qBodyToEnu.y;
    pkt[7] = qBodyToEnu.z;
    pkt[8] = thrustEnu.e;
    pkt[9] = thrustEnu.n;
    pkt[10] = thrustEnu.u;
    pkt[11] = mass;
    pkt[12] = velEnu.e;
    pkt[13] = velEnu.n;
    pkt[14] = velEnu.u;
    pkt[15] = (double)mode;
    pkt[16] = (double)m_maskAccum;
    m_maskAccum = 0;

    m_sender.send(m_cfg.host, m_cfg.port, pkt, (int)sizeof(pkt));
}

void UdpTelemetry::sendV3(double t, VecEnu posEnu, Quat qBodyToEnu,
                          VecEnu thrustEnu, double mass, VecEnu velEnu,
                          int mode, int rcsMask, double dFeet)
{
    if (!m_sender.isOpen())
    {
        return;
    }
    m_maskAccum |= rcsMask;
    if (m_counter++ % m_cfg.decimation != 0)
    {
        return;
    }
    double pkt[PKT3_DOUBLES];
    pkt[0] = t;
    pkt[1] = posEnu.e;
    pkt[2] = posEnu.n;
    pkt[3] = posEnu.u;
    pkt[4] = qBodyToEnu.w;
    pkt[5] = qBodyToEnu.x;
    pkt[6] = qBodyToEnu.y;
    pkt[7] = qBodyToEnu.z;
    pkt[8] = thrustEnu.e;
    pkt[9] = thrustEnu.n;
    pkt[10] = thrustEnu.u;
    pkt[11] = mass;
    pkt[12] = velEnu.e;
    pkt[13] = velEnu.n;
    pkt[14] = velEnu.u;
    pkt[15] = (double)mode;
    pkt[16] = (double)m_maskAccum;
    pkt[17] = dFeet;
    m_maskAccum = 0;
    m_sender.send(m_cfg.host, m_cfg.port, pkt, (int)sizeof(pkt));
}

void UdpTelemetry::sendV4(double t, VecEnu posEnu, Quat qBodyToEnu,
                          VecEnu thrustEnu, double mass, VecEnu velEnu,
                          int mode, int rcsMask, double dFeet,
                          double n2Left, double propLeft)
{
    if (!m_sender.isOpen())
    {
        return;
    }
    m_maskAccum |= rcsMask;
    if (m_counter++ % m_cfg.decimation != 0)
    {
        return;
    }
    double pkt[PKT4_DOUBLES];
    pkt[0] = t;
    pkt[1] = posEnu.e;
    pkt[2] = posEnu.n;
    pkt[3] = posEnu.u;
    pkt[4] = qBodyToEnu.w;
    pkt[5] = qBodyToEnu.x;
    pkt[6] = qBodyToEnu.y;
    pkt[7] = qBodyToEnu.z;
    pkt[8] = thrustEnu.e;
    pkt[9] = thrustEnu.n;
    pkt[10] = thrustEnu.u;
    pkt[11] = mass;
    pkt[12] = velEnu.e;
    pkt[13] = velEnu.n;
    pkt[14] = velEnu.u;
    pkt[15] = (double)mode;
    pkt[16] = (double)m_maskAccum;
    pkt[17] = dFeet;
    pkt[18] = n2Left;
    pkt[19] = propLeft;
    m_maskAccum = 0;
    m_sender.send(m_cfg.host, m_cfg.port, pkt, (int)sizeof(pkt));
}

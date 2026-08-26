#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

// lib_udp_posix.cpp -- реализация UdpReader / UdpSender для POSIX
// (BSD-сокеты).
#ifndef _WIN32

#include "lib_udp.h"

#include <string.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>

UdpReader::UdpReader() : m_sock(-1), m_open(false)
{
}

UdpReader::~UdpReader()
{
    close();
}

bool UdpReader::open(int port, bool loopbackOnly, int timeoutMs)
{
    int s = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    if (s < 0)
    {
        return false;
    }
    struct sockaddr_in sa;
    memset(&sa, 0, sizeof(sa));
    sa.sin_family = AF_INET;
    sa.sin_port = htons((unsigned short)port);
    sa.sin_addr.s_addr = htonl(loopbackOnly ? INADDR_LOOPBACK : INADDR_ANY);
    if (bind(s, (struct sockaddr*)&sa, sizeof(sa)) < 0)
    {
        ::close(s);
        return false;
    }
    if (timeoutMs > 0)
    {
        struct timeval tv;
        tv.tv_sec = timeoutMs / 1000;
        tv.tv_usec = (timeoutMs % 1000) * 1000;
        setsockopt(s, SOL_SOCKET, SO_RCVTIMEO, &tv, sizeof(tv));
    }
    m_sock = (long long)s;
    m_open = true;
    return true;
}

void UdpReader::close()
{
    if (m_open)
    {
        ::close((int)m_sock);
        m_open = false;
    }
}

int UdpReader::read(void* buf, int maxLen)
{
    if (!m_open)
    {
        return -1;
    }
    return (int)::recv((int)m_sock, buf, (size_t)maxLen, 0);
}

UdpSender::UdpSender() : m_sock(-1), m_open(false)
{
}

UdpSender::~UdpSender()
{
    close();
}

bool UdpSender::open()
{
    int s = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    if (s < 0)
    {
        return false;
    }
    m_sock = (long long)s;
    m_open = true;
    return true;
}

void UdpSender::close()
{
    if (m_open)
    {
        ::close((int)m_sock);
        m_open = false;
    }
}

bool UdpSender::send(const char* host, int port, const void* buf, int len)
{
    if (!m_open)
    {
        return false;
    }
    struct sockaddr_in sa;
    memset(&sa, 0, sizeof(sa));
    sa.sin_family = AF_INET;
    sa.sin_port = htons((unsigned short)port);
    if (host == 0)
    {
        host = "127.0.0.1";
    }
    if (inet_pton(AF_INET, host, &sa.sin_addr) != 1)
    {
        return false;
    }
    return ::sendto((int)m_sock, buf, (size_t)len, 0, (struct sockaddr*)&sa,
                    sizeof(sa)) == len;
}

#endif // !_WIN32

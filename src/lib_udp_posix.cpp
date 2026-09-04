//======================================================================
//  lib_udp_posix.cpp - UdpReader / UdpSender for POSIX (BSD sockets).
//======================================================================

#ifndef _WIN32

#include "lib_udp.h"

#include <string.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>

namespace gnc
{

UdpReader::UdpReader()
    : m_sock(-1)
    , m_open(false)
{
}

UdpReader::~UdpReader()
{
    Close();
}

bool UdpReader::Open(int port, bool loopbackOnly, int timeoutMs)
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
    unsigned long addr = INADDR_ANY;
    if (loopbackOnly)
    {
        addr = INADDR_LOOPBACK;
    }
    sa.sin_addr.s_addr = htonl(addr);
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

void UdpReader::Close()
{
    if (m_open)
    {
        ::close((int)m_sock);
        m_open = false;
    }
}

int UdpReader::Read(void* buf, int maxLen)
{
    if (!m_open)
    {
        return -1;
    }
    return (int)::recv((int)m_sock, buf, (size_t)maxLen, 0);
}

UdpSender::UdpSender()
    : m_sock(-1)
    , m_open(false)
{
}

UdpSender::~UdpSender()
{
    Close();
}

bool UdpSender::Open()
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

void UdpSender::Close()
{
    if (m_open)
    {
        ::close((int)m_sock);
        m_open = false;
    }
}

bool UdpSender::Send(const char* host, int port, const void* buf, int len)
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
    long long n = (long long)::sendto((int)m_sock, buf, (size_t)len, 0,
                                      (struct sockaddr*)&sa, sizeof(sa));
    return n == (long long)len;
}

}   // namespace gnc

#endif  // !_WIN32

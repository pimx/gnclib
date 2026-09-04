//======================================================================
//  lib_udp_windows.cpp - UdpReader / UdpSender for Windows (WinSock2;
//  WSAStartup/WSACleanup per object - Windows keeps a reference count).
//======================================================================

#ifdef _WIN32

#include "lib_udp.h"

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <winsock2.h>
#include <ws2tcpip.h>
#include <string.h>
#pragma comment(lib, "ws2_32.lib")

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
    WSADATA wsa;
    if (WSAStartup(MAKEWORD(2, 2), &wsa) != 0)
    {
        return false;
    }
    SOCKET s = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    if (s == INVALID_SOCKET)
    {
        WSACleanup();
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
    if (bind(s, (struct sockaddr*)&sa, sizeof(sa)) != 0)
    {
        closesocket(s);
        WSACleanup();
        return false;
    }
    if (timeoutMs > 0)
    {
        DWORD tmo = (DWORD)timeoutMs;
        setsockopt(s, SOL_SOCKET, SO_RCVTIMEO, (const char*)&tmo, sizeof(tmo));
    }
    m_sock = (long long)s;
    m_open = true;
    return true;
}

void UdpReader::Close()
{
    if (m_open)
    {
        closesocket((SOCKET)m_sock);
        WSACleanup();
        m_open = false;
    }
}

int UdpReader::Read(void* buf, int maxLen)
{
    if (!m_open)
    {
        return -1;
    }
    return recv((SOCKET)m_sock, (char*)buf, maxLen, 0);
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
    WSADATA wsa;
    if (WSAStartup(MAKEWORD(2, 2), &wsa) != 0)
    {
        return false;
    }
    SOCKET s = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    if (s == INVALID_SOCKET)
    {
        WSACleanup();
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
        closesocket((SOCKET)m_sock);
        WSACleanup();
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
    int n = sendto((SOCKET)m_sock, (const char*)buf, len, 0,
                   (struct sockaddr*)&sa, sizeof(sa));
    return n == len;
}

}   // namespace gnc

#endif  // _WIN32

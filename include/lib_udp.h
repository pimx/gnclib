//======================================================================
//  lib_udp.h
//  UDP exchange classes. The interface is OS-invariant and pulls no OS
//  headers; the implementations live in udp_windows.cpp (WinSock2,
//  WSAStartup per object - Windows reference-counts it) and
//  udp_posix.cpp (BSD sockets). Both files are guarded by #ifdef and
//  can be fed to the build system unconditionally.
//
//  UdpReader - datagram reception: a socket bound to a port
//  (optionally loopback only) with a configurable receive timeout
//  (for loops that poll a stop flag).
//  UdpSender - datagram transmission to host:port (host 0 -> 127.0.0.1).
//
//  The classes own a socket, so they are non-copyable (a copy would
//  double-close the descriptor). C++17, no STL, no heap, no exceptions.
//======================================================================

#ifndef LIB_UDP_H
#define LIB_UDP_H

namespace gnc
{

class UdpReader
{
public:

    UdpReader();
    ~UdpReader();

    //  Open reception on a port; loopbackOnly - bind to 127.0.0.1;
    //  timeoutMs - receive timeout (0 - blocking reception).
    bool Open(int port, bool loopbackOnly, int timeoutMs);

    void Close();

    bool IsOpen() const
    {
        return m_open;
    }

    //  Receive one datagram; returns its length, -1 - timeout or error.
    int Read(void* buf, int maxLen);

private:

    UdpReader(const UdpReader&) = delete;
    UdpReader& operator=(const UdpReader&) = delete;

    long long m_sock;   // socket descriptor (OS-dependent)
    bool m_open;
};

class UdpSender
{
public:

    UdpSender();
    ~UdpSender();

    bool Open();

    void Close();

    bool IsOpen() const
    {
        return m_open;
    }

    //  Send a datagram to host:port (host 0 -> "127.0.0.1").
    bool Send(const char* host, int port, const void* buf, int len);

private:

    UdpSender(const UdpSender&) = delete;
    UdpSender& operator=(const UdpSender&) = delete;

    long long m_sock;   // socket descriptor (OS-dependent)
    bool m_open;
};

}   // namespace gnc

#endif  // LIB_UDP_H

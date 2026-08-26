//============================================================================
// lib_udp.h -- библиотечные классы обмена по UDP.
//
// UdpReader -- приём датаграмм: сокет, привязанный к порту (опционально
// только loopback), с настраиваемым таймаутом приёма (для циклов с
// проверкой флага останова).
// UdpSender -- отправка датаграмм на host:port (host 0 -> 127.0.0.1).
//
// Интерфейс инвариантен к операционной системе; реализации -- в файлах
// с суффиксами ОС: lib_udp_windows.cpp (WinSock2, WSAStartup на объект),
// lib_udp_posix.cpp (BSD-сокеты). Без STL, исключений и динамических
// аллокаций.
//============================================================================
#ifndef LIB_UDP_H
#define LIB_UDP_H

class UdpReader
{
public:
    UdpReader();
    ~UdpReader();

    // Открыть приём на порту; loopbackOnly -- привязка к 127.0.0.1;
    // timeoutMs -- таймаут приёма (0 -- блокирующий приём)
    bool open(int port, bool loopbackOnly, int timeoutMs);
    void close();
    bool isOpen() const
    {
        return m_open;
    }

    // Приём одной датаграммы; возврат длины, -1 -- таймаут или ошибка
    int read(void* buf, int maxLen);

private:
    // --- состояние сокета ---
    long long m_sock; // описатель сокета (ОС-зависимый)
    bool m_open;      // сокет открыт
};

class UdpSender
{
public:
    UdpSender();
    ~UdpSender();

    // Открыть сокет отправки
    bool open();
    void close();
    bool isOpen() const
    {
        return m_open;
    }

    // Отправить датаграмму на host:port (host 0 -> "127.0.0.1")
    bool send(const char* host, int port, const void* buf, int len);

private:
    // --- состояние сокета ---
    long long m_sock; // описатель сокета (ОС-зависимый)
    bool m_open;      // сокет открыт
};

#endif // LIB_UDP_H

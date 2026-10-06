// test_ini.cpp -- gnclib 1.2.0: разбор файла параметров Ini и значений
// ключей (Cli::parseEndpoint, listItem, splitPair) -- lib_cli.
// Build: g++ -std=c++17 -O2 -Iinclude tests/test_ini.cpp build/libgnc.a -o build/test_ini
#include "lib_cli.h"
#include <stdio.h>
#include <string.h>

static int g_bad = 0;

// Проверка условия с печатью.
static void Check(const char* name, bool ok)
{
    printf("%-50s %s\n", name, ok ? "OK" : "FAIL");
    g_bad += ok ? 0 : 1;
}

// Точка входа.
int main()
{
    static Ini ini;
    const char* text =
        "; comment line\n"
        "top = 1.5   ; trailing comment\n"
        "[vehicle]\r\n"
        "\tm0\t=  23000.0 # hash comment\r\n"
        "name = HOP90\n"
        "no equals sign here\n"
        "[ sim ]\n"
        "m0 = 21000\n"
        "empty =\n"
        "  = 5\n";
    Check("loadText", ini.loadText(text));
    Check("count == 5", ini.count() == 5);
    Check("top = 1.5 before any section", ini.getDouble("top", 0.0) == 1.5
                                           && strcmp(ini.section(0), "") == 0);
    Check("last occurrence wins (m0 = 21000)", ini.getDouble("m0", 0.0) == 21000.0);
    Check("section lookup (vehicle m0 = 23000)", ini.getDouble("m0", 0.0, "vehicle") == 23000.0);
    Check("section name trimmed ([ sim ])", ini.getInt("m0", 0, "sim") == 21000);
    char v[CLI_VALUE_MAX];
    Check("string value", ini.getValue("name", v) && strcmp(v, "HOP90") == 0);
    Check("empty value present", ini.hasKey("empty") && ini.value(4)[0] == 0);
    Check("missing key -> default", ini.getDouble("nokey", 7.0) == 7.0);
    Check("line numbers", ini.line(1) == 4);
    const char* keys[] = { "top", "m0", "name" };
    Check("verifyKeys finds unknown 'empty'", ini.verifyKeys(keys, 3, false) == false);
    Check("load of a missing file fails", ini.load("no_such_file.ini") == false);
    char a[CLI_VALUE_MAX];
    int port = 0;
    Check("endpoint empty -> 127.0.0.1:9000", Cli::parseEndpoint("", a, 64, &port)
                                               && strcmp(a, "127.0.0.1") == 0 && port == 9000);
    Check("endpoint 10.0.0.2:9100", Cli::parseEndpoint("10.0.0.2:9100", a, 64, &port)
                                     && strcmp(a, "10.0.0.2") == 0 && port == 9100);
    Check("endpoint port only 9200", Cli::parseEndpoint("9200", a, 64, &port) && port == 9200);
    Check("endpoint addr only", Cli::parseEndpoint("192.168.1.5", a, 64, &port)
                                 && strcmp(a, "192.168.1.5") == 0 && port == 9000);
    Check("endpoint :9300", Cli::parseEndpoint(":9300", a, 64, &port) && port == 9300);
    Check("endpoint bad address rejected", Cli::parseEndpoint("host:9000", a, 64, &port) == false);
    Check("endpoint bad port rejected", Cli::parseEndpoint("1.2.3.4:70000", a, 64, &port) == false);
    char it[CLI_VALUE_MAX];
    Check("listItem 0 / 2", Cli::listItem("dyn,sens,act", 0, it) && strcmp(it, "dyn") == 0
                            && Cli::listItem("dyn,sens,act", 2, it) && strcmp(it, "act") == 0);
    Check("listItem past end", Cli::listItem("dyn,sens,act", 3, it) == false);
    char nm[CLI_VALUE_MAX];
    Check("splitPair ctl=10.0.0.3", Cli::splitPair("ctl=10.0.0.3", nm, it)
                                    && strcmp(nm, "ctl") == 0 && strcmp(it, "10.0.0.3") == 0);
    Check("splitPair without '='", Cli::splitPair("ctl", nm, it) == false);
    printf(g_bad == 0 ? "TEST_INI PASS\n" : "TEST_INI FAIL\n");
    return g_bad == 0 ? 0 : 1;
}

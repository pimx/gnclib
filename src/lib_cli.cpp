// lib_cli.cpp -- реализация разбора командной строки (см. lib_cli.h).
#include "lib_cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Cli::Cli() : m_defs(0), m_nDefs(0), m_nKeys(0), m_nPos(0)
{
}

Cli::Cli(int argc, char** argv) : m_defs(0), m_nDefs(0), m_nKeys(0), m_nPos(0)
{
    init(argc, argv);
}

// Токен -- ключ: '-' или '/' и далее ИМЯ ТОЛЬКО ИЗ БУКВ (до '=' или
// конца токена, не короче одной буквы). Следствия: отрицательные числа
// ("-2.0") -- значения; пути ("/tmp/app", "/usr/bin/x") -- позиционные
// ('/' -- не буква); токены вида "/data.2=x", "/a-b=c" -- позиционные
// (в имени не-буквы). Значение после '=' может содержать любые символы
bool Cli::isKeyToken(const char* s)
{
    if (s == 0 || (s[0] != '-' && s[0] != '/'))
    {
        return false;
    }
    if (s[1] == 0 || s[1] == '=')
    {
        return false; // пустое имя
    }
    {
        char c0 = s[1];
        if (!((c0 >= 'a' && c0 <= 'z') || (c0 >= 'A' && c0 <= 'Z')))
        {
            return false; // имя обязано НАЧИНАТЬСЯ с буквы
        }
    }
    for (int i = 2; s[i] != 0 && s[i] != '='; ++i)
    {
        char c = s[i];
        bool letter = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
        bool digit = (c >= '0' && c <= '9');
        if (!letter && !digit)
        {
            return false; // в имени посторонний символ: не ключ
        }
    }
    return true;
}

void Cli::init(int argc, char** argv)
{
    m_defs = 0;
    m_nDefs = 0;
    m_nKeys = 0;
    m_nPos = 0;
    int i = 0;
    while (i < argc)
    {
        const char* a = argv[i];
        if (isKeyToken(a))
        {
            const char* name = a + 1; // после '-' или '/'
            const char* eq = strchr(name, '=');
            const char* value = 0;
            int nameLen;
            if (eq != 0)
            {
                nameLen = (int)(eq - name); // /key=value
                value = eq + 1;
            }
            else
            {
                nameLen = (int)strlen(name);
                // значение -- следующий токен, если он не ключ
                if (i + 1 < argc && !isKeyToken(argv[i + 1]))
                {
                    value = argv[i + 1];
                    i = i + 1; // токен поглощён как значение
                }
            }
            if (m_nKeys < CLI_MAX_KEYS)
            {
                m_keys[m_nKeys].name = name;
                m_keys[m_nKeys].nameLen = nameLen;
                m_keys[m_nKeys].value = value;
                m_nKeys = m_nKeys + 1;
            }
        }
        else
        {
            if (m_nPos < CLI_MAX_POS)
            {
                m_pos[m_nPos] = a;
                m_nPos = m_nPos + 1;
            }
        }
        i = i + 1;
    }
}

// Имя из строки (len символов) -- НЕПУСТОЙ ПРЕФИКС запрошенного key,
// регистр не различается. Пропуски букв середины ("-rt" к realtime)
// НЕ допускаются -- при необходимости коротких имён вводить alias
// (прикладной слой), а не ослаблять правило
bool Cli::nameMatches(const char* name, int len, const char* key)
{
    if (len <= 0)
    {
        return false;
    }
    // Ключ, объявленный с ведущим '=' -- ТОЧНОЕ совпадение (V70):
    // не участвует в префиксных сокращениях и не создаёт
    // неоднозначности прежним однобуквенным флагам.
    bool exactOnly = false;
    if (key[0] == '=')
    {
        exactOnly = true;
        key = key + 1;
    }
    for (int i = 0; i < len; ++i)
    {
        char a = name[i];
        char b = key[i];
        if (b == 0)
        {
            return false; // имя в строке длиннее запрошенного
        }
        if (a >= 'A' && a <= 'Z')
        {
            a = (char)(a - 'A' + 'a');
        }
        if (b >= 'A' && b <= 'Z')
        {
            b = (char)(b - 'A' + 'a');
        }
        if (a != b)
        {
            return false;
        }
    }
    if (exactOnly && key[len] != 0)
    {
        return false; // сокращение точного ключа запрещено
    }
    return true;
}

void Cli::copyValue(char dst[], const char* src)
{
    int i = 0;
    if (src != 0)
    {
        while (src[i] != 0 && i < CLI_VALUE_MAX - 1)
        {
            dst[i] = src[i];
            i = i + 1;
        }
    }
    dst[i] = 0;
}

const Cli::Key* Cli::findKey(const char* key) const
{
    // Точный ключ ("=key" в объявленном списке, V70): токен строки
    // должен совпасть с именем ПОЛНОСТЬЮ -- сокращения не принимаются
    // и не перехватывают чужие однобуквенные флаги.
    bool exactOnly = false;
    if (m_defs != 0)
    {
        int klen = (int)strlen(key);
        for (int i = 0; i < m_nDefs; ++i)
        {
            const char* d = m_defs[i];
            if (d != 0 && d[0] == '=' &&
                nameMatches(key, klen, d))
            {
                exactOnly = true;
                break;
            }
        }
    }
    for (int i = 0; i < m_nKeys; ++i)
    {
        if (exactOnly && m_keys[i].nameLen != (int)strlen(key))
        {
            continue;
        }
        if (nameMatches(m_keys[i].name, m_keys[i].nameLen, key))
        {
            return &m_keys[i];
        }
    }
    return 0;
}

bool Cli::getKeyName(int index, char name[]) const
{
    if (index < 0 || index >= m_nKeys)
    {
        name[0] = 0;
        return false;
    }
    int n = m_keys[index].nameLen;
    if (n > CLI_VALUE_MAX - 1)
    {
        n = CLI_VALUE_MAX - 1;
    }
    for (int i = 0; i < n; ++i)
    {
        name[i] = m_keys[index].name[i];
    }
    name[n] = 0;
    return true;
}

// Публичная обёртка правила соответствия имён
bool Cli::nameMatchesKey(const char* name, const char* key)
{
    return nameMatches(name, (int)strlen(name), key);
}

bool Cli::hasKey(const char* key) const
{
    requireUnambiguous(key, "hasKey");
    return findKey(key) != 0;
}

bool Cli::getKeyValue(const char* key, char value[], const char* def) const
{
    requireUnambiguous(key, "getKeyValue");
    const Key* k = findKey(key);
    if (k != 0 && k->value != 0)
    {
        copyValue(value, k->value);
        return true;
    }
    copyValue(value, def);
    return false;
}

bool Cli::getPosValue(int pos, char value[]) const
{
    if (pos < 0 || pos >= m_nPos)
    {
        copyValue(value, "");
        return false;
    }
    copyValue(value, m_pos[pos]);
    return true;
}

// Аварийное завершение с диагностикой: молчаливое срабатывание не того
// ключа обходится дороже остановки задачи
static void cliFatal(const char* what, const char* a, const char* b)
{
    if (b != 0)
    {
        fprintf(stderr, "CLI FATAL: %s: '%s' vs '%s'\n", what, a, b);
    }
    else
    {
        fprintf(stderr, "CLI FATAL: %s: '%s'\n", what, a);
    }
    exit(2);
}

int Cli::matchCount(const char* name) const
{
    if (m_defs == 0 || name == 0)
    {
        return 1; // список не объявлен -- проверок нет
    }
    int len = (int)strlen(name);
    int n = 0;
    for (int i = 0; i < m_nDefs; ++i)
    {
        if (nameMatches(name, len, m_defs[i]))
        {
            n = n + 1;
        }
    }
    return n;
}

// Проверка командной строки БЕЗ аварийного завершения (требование владельца,
// этап 119). Разбор здесь независимый: метод пригоден и до init, и для
// проверки чужого argv. Список объявленных ключей обязателен -- именно он
// задаёт множество допустимых имён.
bool Cli::verifyArgs(int argc, const char** argv, bool verbose) const
{
    if (m_defs == 0 || m_nDefs <= 0)
    {
        if (verbose)
        {
            fprintf(stderr,
                    "CLI verify: no declared key list (init with keys)\n");
        }
        return false;
    }
    int seen[CLI_MAX_DEFS];
    for (int i = 0; i < CLI_MAX_DEFS; ++i)
    {
        seen[i] = 0;
    }
    bool ok = true;
    for (int a = 1; a < argc; ++a)
    {
        const char* s = argv[a];
        if (!isKeyToken(s))
        {
            continue; // значение или позиционный аргумент
        }
        int len = 0;
        while (s[1 + len] != 0 && s[1 + len] != '=')
        {
            ++len;
        }
        char nm[CLI_VALUE_MAX];
        int n = (len > CLI_VALUE_MAX - 1) ? CLI_VALUE_MAX - 1 : len;
        for (int c2 = 0; c2 < n; ++c2)
        {
            nm[c2] = s[1 + c2];
        }
        nm[n] = 0;
        // допустимость: имя -- префикс РОВНО ОДНОГО объявленного ключа
        int hits = 0, which = -1;
        for (int i = 0; i < m_nDefs; ++i)
        {
            if (nameMatches(nm, n, m_defs[i]))
            {
                ++hits;
                which = i;
            }
        }
        if (hits == 0)
        {
            if (verbose)
            {
                fprintf(stderr, "CLI verify: unknown key '%s'\n", nm);
            }
            ok = false;
            continue;
        }
        if (hits > 1)
        {
            if (verbose)
            {
                fprintf(
                    stderr,
                    "CLI verify: ambiguous key '%s' (matches %d declared)\n",
                    nm, hits);
            }
            ok = false;
            continue;
        }
        // повтор: то же имя или другое сокращение того же ключа
        if (which >= 0 && which < CLI_MAX_DEFS)
        {
            if (seen[which] != 0)
            {
                if (verbose)
                {
                    fprintf(stderr,
                            "CLI verify: duplicate key '%s' (same as declared "
                            "'%s')\n",
                            nm, m_defs[which]);
                }
                ok = false;
            }
            seen[which] = 1;
        }
    }
    return ok;
}

void Cli::requireUnambiguous(const char* name, const char* where) const
{
    if (m_defs == 0)
    {
        return;
    }
    int n = matchCount(name);
    if (n > 1)
    {
        fprintf(stderr,
                "CLI FATAL: ambiguous key name in %s: '%s' matches %d declared "
                "keys:",
                where, name, n);
        int len = (int)strlen(name);
        for (int i = 0; i < m_nDefs; ++i)
        {
            if (nameMatches(name, len, m_defs[i]))
            {
                fprintf(stderr, " %s", m_defs[i]);
            }
        }
        fprintf(stderr, "\n");
        exit(2);
    }
}

void Cli::init(int argc, char** argv, const char* const* keys, int nkeys,
               bool fatalChecks)
{
    init(argc, argv);
    m_defs = keys;
    m_nDefs = (nkeys < CLI_MAX_DEFS) ? nkeys : CLI_MAX_DEFS;
    if (keys == 0 || nkeys <= 0)
    {
        m_defs = 0;
        m_nDefs = 0;
        return;
    }
    if (nkeys > CLI_MAX_DEFS)
    {
        cliFatal("declared key list exceeds CLI_MAX_DEFS", "", 0);
    }
    if (!fatalChecks)
    {
        return; // этап 119: контроль строки берёт на себя verifyArgs
    }
    // (1) список не должен содержать ключ, являющийся префиксом другого
    for (int i = 0; i < m_nDefs; ++i)
    {
        if (m_defs[i] == 0 || m_defs[i][0] == 0)
        {
            cliFatal("empty declared key name", "", 0);
        }
        for (int j = 0; j < m_nDefs; ++j)
        {
            if (i == j)
            {
                continue;
            }
            int li = (int)strlen(m_defs[i]);
            if (nameMatches(m_defs[i], li, m_defs[j]))
            {
                cliFatal("declared key is a prefix of another", m_defs[i],
                         m_defs[j]);
            }
        }
    }
    // (2) каждое имя из командной строки должно быть однозначным
    for (int i = 0; i < m_nKeys; ++i)
    {
        char nm[CLI_VALUE_MAX];
        int n = m_keys[i].nameLen;
        if (n > CLI_VALUE_MAX - 1)
        {
            n = CLI_VALUE_MAX - 1;
        }
        for (int c = 0; c < n; ++c)
        {
            nm[c] = m_keys[i].name[c];
        }
        nm[n] = 0;
        if (matchCount(nm) == 0)
        {
            cliFatal("unknown key on the command line", nm, 0);
        }
        requireUnambiguous(nm, "command line");
    }
}

//============================================================================
// Значения ключей (1.2.0): адрес, списки, пары.
//============================================================================

// Адрес IPv4 a.b.c.d (каждая часть 0..255), n символов s; 1 -- верно.
static int cliIpv4(const char* s, int n)
{
    int parts = 0;
    int val = -1;
    int i = 0;
    for (i = 0; i <= n; i++)
    {
        char c = (i < n) ? s[i] : 0;
        if (c >= '0' && c <= '9')
        {
            val = (val < 0) ? (c - '0') : val * 10 + (c - '0');
            if (val > 255)
            {
                return 0;
            }
        }
        else if (c == '.' || c == 0)
        {
            if (val < 0)
            {
                return 0;
            }
            parts++;
            val = -1;
        }
        else
        {
            return 0;
        }
    }
    return (parts == 4) ? 1 : 0;
}

// Сетевой адрес из значения ключа (см. lib_cli.h; перенесено из
// tlm90Endpoint проектов без изменения логики).
bool Cli::parseEndpoint(const char* s, char addr[], int addrSize, int* port,
                        const char* defAddr, int defPort)
{
    int a = 0;
    int la = 0;
    int p = 0;
    int lp = 0;
    int n = 0;
    int colon = -1;
    int dot = 0;
    int i = 0;
    if (s != 0)
    {
        while (s[n] != 0)
        {
            if (s[n] == ':' && colon < 0)
            {
                colon = n;
            }
            if (s[n] == '.')
            {
                dot = 1;
            }
            n++;
        }
    }
    if (colon >= 0)
    {
        la = colon;
        p = colon + 1;
        lp = n - p;
    }
    else if (dot != 0)
    {
        la = n;
    }
    else
    {
        p = 0;
        lp = n;
    }
    *port = defPort;
    if (lp > 0)
    {
        long v = 0;
        for (i = 0; i < lp; i++)
        {
            char c = s[p + i];
            if (c < '0' || c > '9' || v > 65535)
            {
                return false;
            }
            v = v * 10 + (c - '0');
        }
        if (v < 1 || v > 65535)
        {
            return false;
        }
        *port = (int)v;
    }
    if (la > 0)
    {
        if (cliIpv4(s + a, la) == 0 || la >= addrSize)
        {
            return false;
        }
        for (i = 0; i < la; i++)
        {
            addr[i] = s[a + i];
        }
        addr[la] = 0;
    }
    else
    {
        for (i = 0; defAddr[i] != 0 && i < addrSize - 1; i++)
        {
            addr[i] = defAddr[i];
        }
        addr[i] = 0;
    }
    return true;
}

// Элемент index списка через запятую (без пробелов по краям не
// обрезается -- как прежние разборщики FLY90).
bool Cli::listItem(const char* s, int index, char item[])
{
    item[0] = 0;
    if (s == 0 || s[0] == 0 || index < 0)
    {
        return false;
    }
    int k = 0;
    const char* p = s;
    while (k < index)
    {
        const char* q = strchr(p, ',');
        if (q == 0)
        {
            return false;
        }
        p = q + 1;
        k++;
    }
    if (*p == 0)
    {
        return false;
    }
    const char* e = strchr(p, ',');
    int len = (e != 0) ? (int)(e - p) : (int)strlen(p);
    len = (len < CLI_VALUE_MAX - 1) ? len : CLI_VALUE_MAX - 1;
    memcpy(item, p, (size_t)len);
    item[len] = 0;
    return true;
}

// Пара "имя=значение": разделение по первому '='; false -- '=' нет.
bool Cli::splitPair(const char* s, char name[], char value[])
{
    name[0] = 0;
    value[0] = 0;
    const char* eq = (s != 0) ? strchr(s, '=') : 0;
    if (eq == 0)
    {
        return false;
    }
    int len = (int)(eq - s);
    len = (len < CLI_VALUE_MAX - 1) ? len : CLI_VALUE_MAX - 1;
    memcpy(name, s, (size_t)len);
    name[len] = 0;
    strncpy(value, eq + 1, CLI_VALUE_MAX);
    value[CLI_VALUE_MAX - 1] = 0;
    return true;
}

//============================================================================
// Ini (1.2.0) -- разбор файла параметров (см. lib_cli.h).
//============================================================================

// Конструктор: пусто.
Ini::Ini() : m_used(0), m_n(0), m_overflow(0)
{
    m_buf[0] = 0;
}

// Отбрасывание пробелов и табуляций (и \r, \n) по краям; возврат -- начало.
static char* iniTrim(char* s)
{
    while (*s == ' ' || *s == '\t')
    {
        s++;
    }
    int n = (int)strlen(s);
    while (n > 0 && (s[n - 1] == ' ' || s[n - 1] == '\t' || s[n - 1] == '\r' || s[n - 1] == '\n'))
    {
        s[--n] = 0;
    }
    return s;
}

// Копия len символов s в буфер с нулём; смещение или -1 при переполнении.
int Ini::store(const char* s, int len)
{
    if (m_used + len + 1 > INI_BUF)
    {
        m_overflow = 1;
        return -1;
    }
    int at = m_used;
    memcpy(m_buf + at, s, (size_t)len);
    m_buf[at + len] = 0;
    m_used += len + 1;
    return at;
}

// Разбор одной строки (строка портится); sec -- смещение текущей секции.
bool Ini::parseLine(char* s, int lineNo, int* sec)
{
    char* cm = strpbrk(s, ";#");
    if (cm != 0)
    {
        *cm = 0;
    }
    char* t = iniTrim(s);
    if (t[0] == '[')
    {
        char* e = strchr(t, ']');
        int len = (e != 0) ? (int)(e - t - 1) : (int)strlen(t + 1);
        char name[INI_LINE_MAX];
        memcpy(name, t + 1, (size_t)len);
        name[len] = 0;
        char* nm = iniTrim(name);
        int at = store(nm, (int)strlen(nm));
        if (at < 0)
        {
            return false;
        }
        *sec = at;
        return true;
    }
    char* eq = strchr(t, '=');
    if (eq == 0)
    {
        return true;
    }
    *eq = 0;
    char* key = iniTrim(t);
    char* val = iniTrim(eq + 1);
    if (key[0] == 0 || key[0] == '[')
    {
        return true;
    }
    if (m_n >= INI_MAX_ENTRIES)
    {
        m_overflow = 1;
        return false;
    }
    int k = store(key, (int)strlen(key));
    int v = store(val, (int)strlen(val));
    if (k < 0 || v < 0)
    {
        return false;
    }
    m_e[m_n].sec = *sec;
    m_e[m_n].key = k;
    m_e[m_n].val = v;
    m_e[m_n].line = lineNo;
    m_n++;
    return true;
}

// Разбор текста в памяти.
bool Ini::loadText(const char* text)
{
    m_used = 0;
    m_n = 0;
    m_overflow = 0;
    int sec = store("", 0);
    int lineNo = 0;
    const char* p = text;
    while (p != 0 && *p != 0)
    {
        const char* e = strchr(p, '\n');
        int len = (e != 0) ? (int)(e - p) : (int)strlen(p);
        char line[INI_LINE_MAX];
        int n = (len < INI_LINE_MAX - 1) ? len : INI_LINE_MAX - 1;
        memcpy(line, p, (size_t)n);
        line[n] = 0;
        lineNo++;
        if (parseLine(line, lineNo, &sec) == false)
        {
            fprintf(stderr, "INI: limits exceeded at line %d (%d keys, %d bytes)\n", lineNo,
                    INI_MAX_ENTRIES, INI_BUF);
            return false;
        }
        p = (e != 0) ? e + 1 : 0;
    }
    return true;
}

// Загрузка файла.
bool Ini::load(const char* path)
{
    FILE* f = (path != 0) ? fopen(path, "rb") : 0;
    if (f == 0)
    {
        m_n = 0;
        return false;
    }
    m_used = 0;
    m_n = 0;
    m_overflow = 0;
    int sec = store("", 0);
    int lineNo = 0;
    char line[INI_LINE_MAX];
    bool ok = true;
    while (ok && fgets(line, sizeof(line), f) != 0)
    {
        lineNo++;
        ok = parseLine(line, lineNo, &sec);
    }
    fclose(f);
    if (!ok)
    {
        fprintf(stderr, "INI: %s: limits exceeded at line %d (%d keys, %d bytes)\n", path, lineNo,
                INI_MAX_ENTRIES, INI_BUF);
    }
    return ok;
}

// Ключ по индексу.
const char* Ini::key(int i) const
{
    return (i >= 0 && i < m_n) ? m_buf + m_e[i].key : "";
}

// Значение по индексу.
const char* Ini::value(int i) const
{
    return (i >= 0 && i < m_n) ? m_buf + m_e[i].val : "";
}

// Секция по индексу.
const char* Ini::section(int i) const
{
    return (i >= 0 && i < m_n) ? m_buf + m_e[i].sec : "";
}

// Номер строки по индексу.
int Ini::line(int i) const
{
    return (i >= 0 && i < m_n) ? m_e[i].line : 0;
}

// Значение числом (atof).
double Ini::number(int i) const
{
    return atof(value(i));
}

// Значение целым (atoi).
int Ini::integer(int i) const
{
    return atoi(value(i));
}

// Индекс последнего вхождения ключа или -1.
int Ini::find(const char* k, const char* sec) const
{
    int i = 0;
    for (i = m_n - 1; i >= 0; i--)
    {
        if (strcmp(m_buf + m_e[i].key, k) == 0
            && (sec == 0 || strcmp(m_buf + m_e[i].sec, sec) == 0))
        {
            return i;
        }
    }
    return -1;
}

// Наличие ключа.
bool Ini::hasKey(const char* k, const char* sec) const
{
    return find(k, sec) >= 0;
}

// Значение ключа строкой.
bool Ini::getValue(const char* k, char out[], const char* def, const char* sec) const
{
    int i = find(k, sec);
    const char* src = (i >= 0) ? value(i) : ((def != 0) ? def : "");
    strncpy(out, src, CLI_VALUE_MAX);
    out[CLI_VALUE_MAX - 1] = 0;
    return i >= 0;
}

// Значение ключа числом.
double Ini::getDouble(const char* k, double def, const char* sec) const
{
    int i = find(k, sec);
    return (i >= 0) ? number(i) : def;
}

// Значение ключа целым.
int Ini::getInt(const char* k, int def, const char* sec) const
{
    int i = find(k, sec);
    return (i >= 0) ? integer(i) : def;
}

// Проверка ключей файла по объявленному списку.
bool Ini::verifyKeys(const char* const* keys, int nkeys, bool verbose) const
{
    bool ok = true;
    int i = 0;
    for (i = 0; i < m_n; i++)
    {
        int found = 0;
        int j = 0;
        for (j = 0; j < nkeys && found == 0; j++)
        {
            found = (strcmp(keys[j], key(i)) == 0) ? 1 : 0;
        }
        if (found == 0)
        {
            ok = false;
            if (verbose)
            {
                fprintf(stderr, "INI: unknown key '%s' (line %d)\n", key(i), line(i));
            }
        }
    }
    return ok;
}

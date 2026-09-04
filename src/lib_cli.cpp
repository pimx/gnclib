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

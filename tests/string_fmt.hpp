#pragma once

#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <stdarg.h>  // For va_start, etc.


// from https://stackoverflow.com/a/8098080
inline std::string vstrfmt(const std::string &fmt_str, va_list ap) {
    int final_n;
    size_t n = fmt_str.size() * 2; /* Reserve two times as much as the length of the fmt_str */
    if (n < 16) {
        n = 16;
    }

    std::string buffer;
    while(true) {
        buffer.resize(n, '\0');
        strcpy(const_cast<char *>(buffer.c_str()), fmt_str.c_str());
        final_n = vsnprintf(const_cast<char *>(buffer.c_str()), n, fmt_str.c_str(), ap);
        if (final_n < 0 || final_n >= (int)n) {
            n += abs(final_n - (int)n + 1);
        } else {
            break;
        }
    }
    return std::string(buffer.c_str(), (size_t)final_n);
}


inline std::string strfmt(const std::string fmt_str, ...)
{
    va_list ap;
    va_start(ap, fmt_str);
    std::string ret = vstrfmt(fmt_str, ap);
    va_end(ap);
    return ret;
}

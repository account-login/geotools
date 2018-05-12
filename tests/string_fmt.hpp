#pragma once

#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <stdarg.h>  // For va_start, etc.


namespace tz {

    // from https://stackoverflow.com/a/8098080
    inline std::string vstrfmt(const std::string &fmt_str, va_list ap) {
        int final_n;
        size_t n = fmt_str.size() * 2; /* Reserve two times as much as the length of the fmt_str */
        if (n < 16) {
            n = 16;
        }

        std::vector<char> buffer;
        while (true) {
            // init buffer
            buffer.resize(n, '\0');

            // copy va_list for reuse
            va_list ap_copy;
            va_copy(ap_copy, ap);
            final_n = vsnprintf(buffer.data(), n, fmt_str.c_str(), ap_copy);
            va_end(ap_copy);

            if (final_n >= (int)n) {
                n += abs(final_n - (int)n + 1);
            } else if (final_n < 0) {
                // vsnprintf failed!
                return fmt_str;
            } else {
                // finished
                break;
            }
        }
        return std::string(buffer.data(), (size_t)final_n);
    }

    // FIXME: fmt_str should be const char *
    inline std::string strfmt(const std::string fmt_str, ...) {
        va_list ap;
        va_start(ap, fmt_str);
        std::string ret = vstrfmt(fmt_str, ap);
        va_end(ap);
        return ret;
    }

    template <class T>
    inline std::string str(const T &value) {
        std::stringstream ss;
        ss << value;
        return ss.str();
    }

    // NOTE: uint8_t and int8_t is considered char type by stringstream
    template <>
    inline std::string str(const uint8_t &value) {
        return str(uint32_t(value));
    }

    template <>
    inline std::string str(const int8_t &value) {
        return str(int32_t(value));
    }

    template <class C>
    inline std::string repr_set(const C &c, size_t limit = 5) {
        std::string ans = "{";
        const char *comma = "";
        size_t count = 0;
        for (typename C::const_iterator it = c.begin(); it != c.end(); ++it) {
            if (count++ < limit) {
                ans += comma + str(*it);
                comma = ", ";
            } else {
                ans += comma + std::string("...");
                break;
            }
        }
        ans += "}";
        return ans;
    }

    template <class M>
    inline std::string repr_map(const M &mapping) {
        std::string ans;
        for (typename M::const_iterator it = mapping.begin(); it != mapping.end(); ++it) {
            ans += "[" + str(it->first) + ":" + str(it->second) + "]";
        }
        return ans;
    }

}   // ::tz

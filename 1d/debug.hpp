#pragma once

#ifndef __NPDEBUG__
#define __NPDEBUG__

#include <iostream>
#include <sstream>

#ifndef NDEBUG
    #pragma message "Compiling in debug mode"
    #define __FILENAME__ (\
        __builtin_strrchr(__FILE__, '/') ? \
        __builtin_strrchr(__FILE__, '/') + 1 : __FILE__)

    #define npdebug_prep(); { \
            std::cerr << "[" << __FILENAME__ \
                      << ":" << __LINE__ \
                      << ", " << __func__ \
                      << "] " ; \
    }

    #define npdebug(...); { \
        npdebug_prep(); \
        np::va_debug(__VA_ARGS__); \
    }

    namespace np {
        template<typename... Args>
        inline void va_debug(Args&&... args) {
            (std::cerr << ... << args) << std::endl;
        }

        template<typename T>
        void range_debug(const T& t) {
            range_debug("", t);
        }

        template<typename T>
        void range_debug(const std::string& msg, const T& t) {
            std::string out;
            for (auto elem : t)
                out += elem += ", ";

            npdebug(msg, out);
        }

        template<typename T>
        T inspect(const T& t) {
            npdebug(t);
            return t;
        }

        template<typename T>
        T inspect(const std::string& msg, const T& t) {
            npdebug(msg, t);
            return t;
        }
    }
#else
    #define npdebug(...) {}

    namespace np {
        template<typename... Args>
        inline void va_debug(Args&... args) {}

        template<typename T>
        inline void range_debug(const T& t) {}

        template<typename T>
        inline void range_debug(const std::string& msg, const T& t) {}

        template<typename T>
        T inspect(const T& t) { return t; }

        template<typename T>
        T inspect(const std::string& msg, const T& t) { return t; }
    }
#endif // NDEBUG
#endif // __NPDEBUG__

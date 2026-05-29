#ifndef RANDOM_H
#define RANDOM_H

#include <ctime>
#include <random>

inline double randDouble(double a, double b) {
    using rt = std::mt19937_64::result_type;

    static const rt seed{
        static_cast<rt>(::time(nullptr)) + static_cast<rt>(::clock())
    };
    static std::mt19937_64 mersenne{seed};
    static std::uniform_real_distribution<double> distribution{a, b};

    return distribution(mersenne);
}

#endif // RANDOM_H

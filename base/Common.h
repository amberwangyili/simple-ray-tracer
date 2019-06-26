//
// Created by YILI WANG on 6/24/19.
//

#ifndef TRACER_FINAL_COMMON_H
#define TRACER_FINAL_COMMON_H

#include <bits/stdc++.h>

const double pi = 3.14159265358979323;
const double epsilon = 1e-6;
const double eps = 1e-6;
const double inf = 1e20;
const double inv_pi = 0.318309886;
const double Gamma = 2.2;
const int prime[] = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
        31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
        73, 79, 83, 89, 97, 101, 103, 107, 109, 113
};

typedef unsigned char uint8;
typedef unsigned int  uint32;
typedef float  float32;
typedef double float64;
typedef std::pair<double, double> Pair;

inline double square(const double x) {
    return x * x;
}


template <typename T>
inline T abs(const T &a) {
    return a < 0 ? -a : a;
}

inline double clamp(double x) {
    if (x < 0) return 0;
    else if (x > 1) return 1;
    else return x;
}

inline double toInt(double x) {
    //return int(clamp(x) * 255. + 0.5);
    return int(pow(clamp(x), 1 / 2.2) * 255 + .5);
}
#endif //TRACER_FINAL_COMMON_H

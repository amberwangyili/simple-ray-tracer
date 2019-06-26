//
// Created by YILI WANG on 6/24/19.
//

#ifndef TRACER_FINAL_RAY_H
#define TRACER_FINAL_RAY_H

#include "Vec.h"

class Ray {

    //成员：origin 出射点，direction 出射方向（为单位向量）
    //操作：给定参数t（在ray上离出射点的距离度量），求得绝对坐标
public:
    Vec origin;
    Vec direction;
    Ray() { }
    Ray(Vec s, Vec d) {
        this->origin = s;
        this->direction = d.normalize();
    }
    inline Vec get(double t) const {
        return origin + direction * t;
    }
};
#endif //TRACER_FINAL_RAY_H

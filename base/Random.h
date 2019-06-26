//
// Created by YILI WANG on 6/24/19.
//

#ifndef TRACER_FINAL_RANDOM_H
#define TRACER_FINAL_RANDOM_H

#include "Common.h"
#include "Ray.h"
using namespace std;
class Random{

    //成员：随机种子
    //操作:
    // 1.get_rand() 线性同余法得到随机数，更新种子
    // 2.norm_coordinate() 得到以一个法向量为z轴的相对坐标系
    // 3.sample_hemisphere() 输入一个面的法向量，得到以该法向量为轴的半球上随机射出的一条光线的方向
    // 4.sample_disk() 输入一个圆盘的中心点，半径，法向量，得到在这个圆盘上随机射出的一条光线（包括它的顶点和方向）
    // 5.sample_triangle() 输入一个三角形三个顶点坐标和法向量，得到这个三角形上随机射出的一条光线（包括它的顶点和方向）
    // 6.sample_sphere() 输入一个球体的半径和圆心，得到从该球体射出的一条随机光线（包括顶点和方向）

public:
    unsigned seed;
public:

    Random(unsigned s = 1){
        seed = s;
    }
    ~Random(){};

    inline double getrand(){
        seed = 1664525u * seed + 1013904223u;
        return seed * (1.0 / 4294967296.0);
    }

    inline void norm_coordinate(const Vec &norm, Vec &u , Vec &v, Vec &w){
        w = norm;
        Vec w_orthon = abs(w.x)>0.1? Vec::YAxis : Vec::XAxis;
        u = cross(w_orthon,w).normalize();
        v = cross(w,u);
    }


    inline Vec sample_hemisphere(const Vec &norm) {
        double phi = 2 * pi * getrand();
        double r2 = getrand();
        double sin_theta = sqrt(r2);
        Vec u,v,w;
        norm_coordinate(norm,u,v,w);
        return (u * cos(phi) *sin_theta+ v * sin(phi) * sin_theta + w * sqrt(1 - r2)).normalize();
    }

    inline Ray sample_disk(const Vec &center, double r, const Vec &norm) {
        double r0 = sqrt(getrand());
        double theta = getrand() * (2.0 * pi);
        double rx = r * r0 * cos(theta);
        double ry = r * r0 * sin(theta);
        Vec u,v,w;
        norm_coordinate(norm,u,v,w);
        Vec origin = center + u * rx + v * ry;
        Vec direct = norm;
        return Ray(origin,direct);
    }

    inline Ray sample_triangle(const Vec &a, const Vec &b, const Vec &c, const Vec &norm) {
        double r1 = getrand(), r2 = getrand();
        if (r1 + r2 >= 1.0) {
            r1 = 1.0 - r1;
            r2 = 1.0 - r2;
        }

        Vec origin = a +  r1 * (b - a) + r2 * (c - a);
        Vec direct = norm;
        return Ray(origin,direct);
    }
    inline Ray sample_sphere(const Vec &center, double r) {
        double cos_theta = 2.0 * getrand() - 1.0;
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        double phi = 2.0 * pi * getrand();
        double x = sin_theta * cos(phi);
        double y = sin_theta * sin(phi);
        double z = cos_theta;
        Vec direct = Vec(x, y, z);
        Vec origin = r * direct + center;
        return Ray(origin,direct);

    }


};
#endif //TRACER_FINAL_RANDOM_H

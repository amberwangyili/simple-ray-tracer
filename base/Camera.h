//
// Created by YILI WANG on 6/24/19.
//

#ifndef TRACER_FINAL_CAMERA_H
#define TRACER_FINAL_CAMERA_H


#include "Common.h"
#include "Vec.h"
#include "Ray.h"
#include "Random.h"




class Camera {
public:
    double fov1, fov_scale1;
    double fov2, fov_scale2;
    Vec eye;
    Vec front, up, right;
    Camera(void) { }

    Camera(const Vec &e, const Vec &f, const Vec &u, double fov1_, double fov2_) {
        eye = e, front = f, up = u, right = cross(f, u);
        fov1 = fov1_, fov_scale1 = tan(fov1_ * 0.5 * pi / 180) * 2;
        fov2 = fov2_, fov_scale2 = tan(fov2_ * 0.5 * pi / 180) * 2;
    }

    inline virtual Ray generate(double x, double y, Random *rng) {
        Vec r = right * ((x - 0.5) * fov_scale1);
        Vec u = up * ((y - 0.5) * fov_scale2);
        return Ray(eye, (front + r + u).normalize());
    }
};

class NaiveCamera : public Camera {
public:
    double fov, fov_scale;

    NaiveCamera(void) { }

    NaiveCamera(const Vec &e, const Vec &f, const Vec &u, double fov_) {
        eye = e, front = f, up = u, right = cross(f, u);
        fov = fov_, fov_scale = tan(fov_ * 0.5 * pi / 180) * 2;
    }

    inline virtual Ray generate(double x, double y, Random *rng) {
        Vec r = right * ((x - 0.5) * fov_scale);
        Vec u = up * ((y - 0.5) * fov_scale);
        return Ray(eye, (front + r + u).normalize());
    }
};


class DofCamera: public Camera {
    //带景深效果的相机
    //成员：视点eye
    //     视线方向front，up一般是(0,1,0)，right是up和front对应的坐标系向量
    //     fov为以front为轴上下视野的夹角
    //     lens_r是起抗锯齿作用的采样圆盘的半径
    //     focus_distance 是焦距
    //操作：生成方射到图像的x,y位置pixel的一条随机光线
public:
    Vec eye;
    Vec front, up, right;
    double fov, fov_scale, lens_r, focus_distance;

    DofCamera(void) { }
    DofCamera(const Vec &e, const Vec &f, const Vec &u, double fov_, double lens_r_, double focus_distance_) {
        eye = e, front = f, up = u, right = cross(f, u);
        fov = fov_, fov_scale = tan(fov_ * 0.5 * pi / 180) * 2;
        lens_r = lens_r_, focus_distance = focus_distance_;
    }

    inline virtual Ray generate(double x, double y, Random *rng) {
        Vec r = right * ((x - 0.5) * fov_scale);
        Vec u = up * ((y - 0.5) * fov_scale);
        Ray center_ray = Ray(eye, (front + r + u).normalize());
        Ray disk_ray = rng->sample_disk(eye, lens_r, front);
        Vec look_at = center_ray.origin + center_ray.direction * focus_distance / dot(center_ray.direction, front);
        disk_ray.direction = (look_at - disk_ray.origin).normalize();
        return disk_ray;
    }

};

class Dof{
    //带景深效果的相机
    //成员：视点eye
    //     视线方向front，up一般是(0,1,0)，right是up和front对应的坐标系向量
    //     fov为以front为轴上下视野的夹角
    //     lens_r是起抗锯齿作用的采样圆盘的半径
    //     focus_distance 是焦距
    //操作：生成方射到图像的x,y位置pixel的一条随机光线
public:
    Vec eye;
    Vec front, up, right;
    double fov, fov_scale, lens_r, focus_distance;

    Dof(void) { }
    Dof(const Vec &e, const Vec &f, const Vec &u, double fov_, double lens_r_, double focus_distance_) {
        eye = e, front = f, up = u, right = cross(f, u);
        fov = fov_, fov_scale = tan(fov_ * 0.5 * pi / 180) * 2;
        lens_r = lens_r_, focus_distance = focus_distance_;
    }

    inline virtual Ray generate(double x, double y, Random *rng) {
        Vec r = right * ((x - 0.5) * fov_scale);
        Vec u = up * ((y - 0.5) * fov_scale);
        Ray center_ray = Ray(eye, (front + r + u).normalize());
        Ray disk_ray = rng->sample_disk(eye, lens_r, front);
        Vec look_at = center_ray.origin + center_ray.direction * focus_distance / dot(center_ray.direction, front);
        disk_ray.direction = (look_at - disk_ray.origin).normalize();
        return disk_ray;
    }

};
#endif //TRACER_FINAL_CAMERA_H

//
// Created by YILI WANG on 6/24/19.
//

#ifndef TRACER_FINAL_IMPLICIT_H
#define TRACER_FINAL_IMPLICIT_H

#include "../object/Object.h"

//third order implicit curve


#define NEWTON_ITER 100
#define NEWTON_DELTA 0.1
#define LAMBDA 2
class Teardrop: public Object{

public:
    AABB box;
    Vec pos;
    double scale;
    Vec vmin, vmax;
public:
    Teardrop(){
        vmax.x = 0.5;
        vmin.x = -0.5;
        vmax.z = 0.5;
        vmin.z =-0.5;
        vmin.y = 0;
        vmax.y = 1;
        box = AABB(vmin,vmax);
    }


    double func(double t,const Ray &ray){
        double x, y, z;
        x = ray.get(t).x;
        y = ray.get(t).y;
        z = ray.get(t).z;
        return x*x+z*z-(y)*(1-y)*(1-y);
    }

    double func_deri(double t, const Ray &r){
        Vec val;
        double x, y, z;
        x = r.get(t).x;
        y = r.get(t).y;
        z = r.get(t).z;
        double dx, dy, dz;
        dx = r.direction.x;
        dy = r.direction.y;
        dz = r.direction.z;
        return 2*x*dx+2*z*dz - (3*y*y-4*y+1)*dy;
    }


    virtual Hitpoint intersect(const Ray &ray);
    virtual void sample(Random *rng, Ray &ray, double &pdf){
        Hitpoint pt = this->intersect(ray);
        Vec dir = rng->sample_hemisphere(pt.norm);
        ray = Ray(pt.position,dir);
        pdf = 1;
    }
};


Hitpoint Teardrop::intersect(const Ray &ray) {

    if (!(box.intersect(ray))) return Hitpoint::null;
    double t;

    int idx=-1;
    float min_dir_val = -1.0f;
    for(int i=0;i<3;i++){
        if(fabs(ray.direction[i])>min_dir_val){
            idx = i;
            min_dir_val = fabs(ray.direction[i]);
        }
    }
    if(ray.direction[idx]>0)
        t = (vmin[idx] - ray.origin[idx]) / ray.direction[idx];
    else
        t = (vmax[idx] - ray.origin[idx]) / ray.direction[idx];

    for (int i = 0; i < NEWTON_ITER; ++i) {
        t = t - NEWTON_DELTA * func(t, ray) / (func_deri(t, ray));
//        cout<<t<<endl;
    }

    if (abs(func(t, ray)) > 1e-2) return Hitpoint::null;
    Hitpoint res;
    res.distance = t;
    res.object = static_cast<Object *>(this);
    res.position = ray.get(t);
    res.norm = Vec(func_deri(t, ray) / (ray.direction.x + eps), func_deri(t, ray) / (ray.direction.y + eps),
                   func_deri(t, ray) / (ray.direction.z + eps));
    res.norm.normalize();
    return res;
//}
}

#endif //TRACER_FINAL_IMPLICIT_H

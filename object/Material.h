//
// Created by YILI WANG on 6/24/19.
//

#ifndef TRACER_FINAL_MATERIAL_H
#define TRACER_FINAL_MATERIAL_H

#include "Texture.h"
#include "../base/Random.h"
#include "../base/Common.h"
#include "../base/Vec.h"




class Light;
//成员: emission (r,g,b)的光亮度
//操作：初始化及get emission

class Material;
//成员：Texture类 reflectance
//操作：sampling虚函数，输入入射ray，入射点位置，法向量，输出出射ray
//     get_reflectance & get_type
class Diffuse;
class Specular;
class Refractive;
class Isotropic;


class Light{
public:
    Light(const Vec &emission_ = Vec(0,0,0)) : emission(emission_){}
    inline Vec get_emission() const {return emission;}
private:
    Vec emission;
};



class Material {
public:
    Material(Texture *const reflectance) : _reflectance(reflectance) {}
    virtual void sampling(const Ray &in, const Vec &pos, const Vec &norm, Random *rng,
                          Ray &out, double &pdf) const = 0;
    virtual Vec get_reflectance(const Vec &pos, const Vec &norm) const {
        return _reflectance->get(pos,norm);
    }
    virtual int get_type()const = 0;


private:
    Texture *_reflectance;
};



class Diffuse : public Material {
public:
    Diffuse(Texture *const reflectance) : Material(reflectance) { }
    virtual void sampling(const Ray &in, const Vec &pos, const Vec &norm, Random *rng,
                          Ray &out, double &pdf) const override {

        Vec abs_norm = (dot(norm, in.direction) < 0) ? norm : -1*norm;
        out.origin = pos;
        out.direction = rng->sample_hemisphere(abs_norm);
        pdf = 1;
    }
    virtual int get_type() const override{
        return 1;
    }


};


class Specular : public Material {
public:
    Specular(Texture *const reflectance) : Material(reflectance) { }

    virtual void sampling(const Ray &in, const Vec &pos, const Vec &norm, Random *rng,
                          Ray &out, double &pdf) const override {

        out.origin = pos;
        out.direction = reflect(in.direction, norm);
        pdf = 1;
    }
    virtual int  get_type() const override{
        return 2;
    }

};

class Refractive : public Material {
public:
    Refractive(Texture *const reflectance, double nt = 1.5) : Material(reflectance), _nt(nt) {

    }


    virtual void sampling(const Ray &in, const Vec &pos, const Vec &norm, Random *rng,
                          Ray &out, double &pdf) const override {

        out.origin = pos;

        Vec refl_dirction = in.direction - norm * 2 * dot(norm, in.direction);
        Vec abs_norm = (dot(norm, in.direction) < 0) ? norm : norm * -1;
        bool   into = dot(norm, abs_norm) > 0;
        double nc = 1, nt = _nt;
        double nnt = into ? nc / nt : nt / nc;
        double ddn = dot(in.direction, abs_norm);
        double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);


        if (cos2t < 0) {
            out.direction = refl_dirction;
            pdf = 1;
        }
        else {
            Vec refr_direction = (in.direction * nnt - norm * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normalize();

            double a = nt - nc, b = nt + nc;
            double R0 = square(a) / square(b);
            double c = 1 - (into ? -ddn : dot(norm, refr_direction));
            double Re = R0 + (1 - R0) * (c * c * c * c * c);
            double Tr = 1 - Re;
            double P = 0.25 + 0.5 * Re;

            if (rng->getrand() < P) {
                out.direction = refl_dirction;
                pdf = P / Re;
            } else {
                out.direction = refr_direction;
                pdf = (1 - P) / Tr;
            }

        }
    }
    virtual int get_type() const override{
        return 3;
    }
private:
    double _nt;
};

class Isotropic: public Material{
public:
    Isotropic(Texture *const reflectance):Material(reflectance){}
    virtual void sampling(const Ray &in, const Vec &pos, const Vec &norm, Random *rng,
                            Ray &out, double &pdf) const override{
        out = rng->sample_sphere(pos,1);
        pdf = 1;
    }
    virtual int get_type() const override{
        return 4;
    }
};
#endif //TRACER_FINAL_MATERIAL_H

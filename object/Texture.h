//
// Created by YILI WANG on 6/24/19.
//

#ifndef TRACER_FINAL_TEXTURE_H
#define TRACER_FINAL_TEXTURE_H
#include "../base/Common.h"
#include "../base/Vec.h"
#include <opencv2/opencv.hpp>



class Texture{
public:
    virtual Vec get(const Vec &pos,const Vec &norm) const = 0;
};

class ConstantTexture: public Texture{
public:
    ConstantTexture(const Vec &c): color(c){}
    virtual Vec get(const Vec &pos,const Vec &norm) const {
        return color;
    }
private:
    Vec color;

};

class CheckerTexture : public Texture {
public:
    CheckerTexture() { }
    CheckerTexture(Texture *t0, Texture *t1): even(t0), odd(t1) { }
    virtual Vec get(const Vec& p,const Vec &norm) const {
        float sines = sin(10*p.x)*sin(10*p.y)*sin(10*p.z);
        if (sines < 0)
            return odd->get(p,norm);
        else
            return even->get(p,norm);
    }
    Texture *odd;
    Texture *even;
};


class ImageTexture : public Texture {
public:
    ImageTexture(std::string filename, const Vec &origin, const Vec &u, const Vec &v, double du, double dv, double gamma = 1.)
            : _origin(origin), _u(u), _v(v), _du(du), _dv(dv) {

        cv::Mat raw = cv::imread(filename, 1);
        cv::Mat rawf;
        raw.convertTo(rawf, CV_32FC3, 1. / 255, 0);
        if (abs(gamma - 1) > eps) {
            cv::pow(rawf, gamma, _picture);
        } else {
            _picture = rawf;
        }
    }

    virtual Vec get(const Vec &pos, const Vec &norm) const override {
        int uu = static_cast<int>(dot((pos - _origin), _u) * _du);
        int vv = static_cast<int>(dot((pos - _origin), _v) * _dv);
        uu %= _picture.cols, vv %= _picture.rows;
        if (uu < 0) uu += _picture.cols;
        if (vv < 0) vv += _picture.rows;
        cv::Vec3f res = _picture.at<cv::Vec3f>(vv, uu);
        return Vec(res[2], res[1], res[0]);
    }

protected:
    cv::Mat _picture;
    Vec _origin, _u, _v;
    double _du, _dv;
};



class BmpTexture : public Texture {
public:
    BmpTexture(std::string filename, const Vec &origin, const Vec &u, const Vec &v, double du, double dv, double gamma = 1.)
            : _origin(origin), _u(u), _v(v), _du(du), _dv(dv) {

        cv::Mat raw = cv::imread(filename, 1);
        cv::Mat rawf;
        raw.convertTo(rawf, CV_32FC3, 1. / 255, 0);
        if (abs(gamma - 1) > eps) {
            cv::pow(rawf, gamma, _picture);
        } else {
            _picture = rawf;
        }
    }

    virtual Vec get(const Vec &pos, const Vec &norm) const override {
        int uu = static_cast<int>(dot((pos - _origin), _u) * _du);
        int vv = static_cast<int>(dot((pos - _origin), _v) * _dv);
        uu %= _picture.cols, vv %= _picture.rows;
        if (uu < 0) uu += _picture.cols;
        if (vv < 0) vv += _picture.rows;
        cv::Vec3f res = _picture.at<cv::Vec3f>(vv, uu);
        return Vec(res[2], res[1], res[0]);
    }

protected:
    cv::Mat _picture;
    Vec _origin, _u, _v;
    double _du, _dv;
};





class BumpTexture : public BmpTexture {
public:
    BumpTexture(std::string filename, const Vec &origin, const Vec &u, const Vec &v, double du, double dv, double gamma = 1., double coeff = 1.)
            : BmpTexture(filename, origin, u, v, du, dv, gamma), _coeff(coeff) {

    }

    virtual Vec get(const Vec &pos, const Vec &norm) const override {
        int uu = static_cast<int>(dot((pos - _origin), _u) * _du);
        int vv = static_cast<int>(dot((pos - _origin), _v) * _dv);

        uu %= _picture.cols, vv %= _picture.rows;
        if (uu < 0) uu += _picture.cols;
        if (vv < 0) vv += _picture.rows;
        cv::Vec3f du, dv;
        if (uu == 0) {
            du = (_picture.at<cv::Vec3f>(vv, uu+1) - _picture.at<cv::Vec3f>(vv, uu)) * 2;
        } else if (uu == _picture.cols - 1) {
            du = (_picture.at<cv::Vec3f>(vv, uu) - _picture.at<cv::Vec3f>(vv, uu-1)) * 2;
        } else {
            du = _picture.at<cv::Vec3f>(vv, uu+1) - _picture.at<cv::Vec3f>(vv, uu-1);
        }
        if (vv == 0) {
            dv = (_picture.at<cv::Vec3f>(vv+1, uu) - _picture.at<cv::Vec3f>(vv, uu)) * 2;
        } else if (vv == _picture.rows - 1) {
            dv = (_picture.at<cv::Vec3f>(vv, uu) - _picture.at<cv::Vec3f>(vv-1, uu)) * 2;
        } else {
            dv = _picture.at<cv::Vec3f>(vv+1, uu) - _picture.at<cv::Vec3f>(vv-1, uu);
        }

        return (norm + du[0] * _u * _coeff + dv[0] * _v * _coeff).normalize();
    }

private:
    double _coeff;
};

#endif //TRACER_FINAL_TEXTURE_H

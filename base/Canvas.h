//
// Created by YILI WANG on 6/24/19.
//

#ifndef TRACER_FINAL_CANVAS_H
#define TRACER_FINAL_CANVAS_H

#include "Common.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>



class Canvas{
public:
    uint8 *img;
    int height, width, c;
    int st0,st1,size;

    Canvas(){}
    Canvas(int h_, int w_ , int c_):height(h_),width(w_),c(c_){
        img = new uint8[h_*w_*c_];
        st1 = c, st0 = c*width, size = c*width*height;
        memset(img,0,size* sizeof(uint8));
    }

    inline uint8 &get(int x, int y, int c){
        return img[x*st0 + y*st1 + c];
    }

    inline void set(int x, int y, int c, uint8 v){
        img[x*st0 + y*st1 + c] = v;
    }

    void show(const std::string &, bool wait = true);
    void write(const std::string &filename);

protected:
    cv::Mat toMat(void);

};

void Canvas::show(const std::string &windowname, bool wait) {
    cv::Mat res = toMat();
    cv::imshow(windowname, res);
    if (wait)
        cv::waitKey(0);
}

void Canvas::write(const std::string &filename) {
    cv::Mat res = toMat();
    cv::imwrite(filename, res);

}

cv::Mat Canvas::toMat() {
    assert(c == 3);
    cv::Mat res(height, width, CV_8UC3);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            res.at<cv::Vec3b>(i, j)[2] = get(i, j, 0);
            res.at<cv::Vec3b>(i, j)[1] = get(i, j, 1);
            res.at<cv::Vec3b>(i, j)[0] = get(i, j, 2);
        }
    }
    return res;
}

#endif //TRACER_FINAL_CANVAS_H

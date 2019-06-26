//
// Created by YILI WANG on 6/24/19.
//

#ifndef TRACER_FINAL_PATHTRACER_H
#define TRACER_FINAL_PATHTRACER_H


#include "../base/Camera.h"
#include "../base/Common.h"
#include "../base/Ray.h"
#include "../base/Canvas.h"
#include "../object/Object.h"
#include "../Scene.h"
#include "../Config.h"




class Pathtracer{
public:
    Pathtracer(Scene &_scene) : scene(_scene){}
    Vec trace(const Ray &ray, int depth, Random *rng);
    void render(Camera &camera, Canvas *canvas, int l=-1, int r=-1);
//    void render(Dof &camera, Canvas *canvas);
protected:
    Scene &scene;

};
//void Pathtracer::render(Dof &camera, Canvas *canvas) {
void Pathtracer::render(Camera &camera, Canvas *canvas, int l, int r) {
    int h = canvas->height, w = canvas->width;
    if (l == -1) l = 0, r = h;
    double t0 = clock();
#pragma omp parallel for schedule(dynamic, 1)

    for (int y = l; y < r ; ++y) {
//    for (int y = 250; y < 400 ; ++y) {
        Random *rng = new Random(19971226+rand());
        double rate = (y - l) / (r - l + 1e-6);
        int t = int((clock() - t0) / rate * (1 - rate) / CLOCKS_PER_SEC);
        fprintf(stderr, "\rRendering (%d spp) %5.2f%% res %02d:%02d:%02d        ", SPP_NUM, 100. * rate, t / 3600, (t / 60) % 60, t % 60);
        for (int x = 0 ; x < w; ++x) {
//        for (int x = 450 ; x < 550; ++x) {
            Vec color = Vec(0,0,0);
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    double sy = 1.0 - (double(y) / h + 1.0 / h * i);
                    double sx = double(x) / w + 1.0 / h * j;
                    Ray ray = camera.generate(sx, sy, rng);
                    Vec sub_color = Vec(0,0,0);
                    for (int k = 0; k < SAMPLE; ++k) {
                        Vec tmp = trace(ray, 0, rng);
                        sub_color = sub_color + tmp / double(SAMPLE);
                    }
                    color = color + Vec(clamp(sub_color.x), clamp(sub_color.y), clamp(sub_color.z))/4;
                }
            }
            canvas->set(y, x, 0, toInt(color.x));
            canvas->set(y, x, 1, toInt(color.y));
            canvas->set(y, x, 2, toInt(color.z));
        }
        delete rng;
    }
}

Vec Pathtracer::trace(const Ray &ray, int depth, Random *rng) {
    Hitpoint pt = scene.intersect(ray);

    if (!pt.object) {
        return Vec(0,0,0);
    }


    //get the intersection object
    Object *object = pt.object;
    Vec pos = pt.position;
    Vec norm = pt.norm;

    //check if it light
    if (object->is_light) {
        return object->light->get_emission();
    }

    //get its material type
    class Material *material = object->material;
    //get its albedo
    Vec reflectance = material->get_reflectance(pos,norm);


    //ref: smallpt
    int new_depth = depth + 1;
    bool is_max_depth = new_depth >= MAX_DEPTH;
    bool use_roulette = new_depth > 5;
    bool roulette = use_roulette && rng->getrand() < reflectance.max();
    if (is_max_depth || (use_roulette && !roulette)) return Vec(0,0,0);

    Vec flux = (use_roulette && roulette) ? reflectance / reflectance.max() : reflectance;

    Ray new_ray; double pdf;
    material->sampling(ray, pos, norm, rng, new_ray, pdf);//get new_rat and pdf

    return flux * trace(new_ray, new_depth, rng) / pdf;
}
#endif //TRACER_FINAL_PATHTRACER_H

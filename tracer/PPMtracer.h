//
// Created by YILI WANG on 6/24/19.
//

#ifndef TRACER_FINAL_PPMTRACER_H
#define TRACER_FINAL_PPMTRACER_H


#include "../Scene.h"
#include "Photon.h"
#include "../base/Canvas.h"
class PPMtracer {
public:
    PPMtracer(Scene *scene, Scene *light, int max_depth)
            : _scene(scene), _light(light), _max_depth(max_depth) {

    }

    ~PPMtracer() {
        delete _photon_map;
    }

    Vec trace(const Ray &ray, int depth, Random *rng, int global_n, double global_r, int caustic_n, double caustic_r);
    virtual void render(Camera *camera, Canvas *canvas);

protected:
    Scene  *_scene, *_light;
    int    _max_depth;

    PhotonMap *_photon_map;
};

void PPMtracer::render(Camera *camera, Canvas *canvas) {
//    for (Object *object : _light->objects)
//        _scene->objects.push_back(object);

    int h = canvas->height, w = canvas->width;
    double alpha = 0.7;
    Vec *colors = new Vec[h * w];
    memset(colors, 0, sizeof(Vec)*h*w);

    for (int k = 0; k < SAMPLE; ++k) {
        _photon_map = new PhotonMap(_scene, _light);
        _photon_map->initialize();

#pragma omp parallel for schedule(dynamic, 1)
        for (int y = 0; y < h; ++y) {
            Random *rng = new Random(19971226+ y + rand());
            if(k%ITER==0)
            fprintf(stderr, "\rRendering (iteration %d) %5.2f%%", k, 100. * y / (h - 1));
            for (int x = 0; x < w; ++x) {
                Vec color = Vec::Zero;
                for (int i = 0; i < 2; ++i) {
                    for (int j = 0; j < 2; ++j) {
                        double sy = 1.0 - (double(y) / h + 1.0 / h * i);
                        double sx = double(x) / w + 1.0 / h * j;
                        Ray ray = camera->generate(sx, sy, rng);

                        Vec sub_color = Vec::Zero;
                        int global_n = PHOTON_GLOBAL_N;
                        int caustic_n = PHOTON_CAUSTIC_N;
                        double global_r = PHOTON_GLOBAL_R;
                        double caustic_r = PHOTON_CAUSTIC_R;
                        global_r = ((i + 1) + alpha) / ((i + 1) + 1.0) * global_r;

                        Vec tmp = trace(ray, 0, rng, global_n, global_r, caustic_n, caustic_r);
                        color += tmp / double(4);
                    }
                }
                colors[y * w + x] += color / double(SAMPLE);
            }
            delete rng;
        }
        double scale = double(SAMPLE) / double(k+1);
        {
            if (k%ITER==0){
                fprintf(stderr, "\rRendering (iteration %d) finished\n", k);
                std::cerr << "rescale factor " << scale << std::endl;
            }

        }


        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                Vec color = colors[y * w + x];
                canvas->set(y, x, 0, toInt(color.x * scale));
                canvas->set(y, x, 1, toInt(color.y * scale));
                canvas->set(y, x, 2, toInt(color.z * scale));
            }
        }
        if(k%ITER==0){
            char buffer[128];
            sprintf(buffer, "result_pm_%d.bmp", k);
            canvas->write(buffer);
            sprintf(buffer, "Photon Mapping Iteration %d", k);
            canvas->show(buffer, false);
        }

    }
}

Vec PPMtracer::trace(const Ray &ray, int depth, Random *rng,
                           int global_n, double global_r, int caustic_n, double caustic_r) {

    Hitpoint inter = _scene->intersect(ray);

    if (!inter.object) {
        return Vec::Zero;
    }

    int new_depth = depth + 1;
    bool is_max_depth = new_depth >= _max_depth;
    Object *object = inter.object;
    Vec pos = inter.position;
    Vec norm = inter.norm;

    if (object->is_light) {
        return object->light->get_emission();
    }
    Material *material = object->material;
    Vec reflectance = material->get_reflectance(pos, norm);

    bool use_rr = new_depth > 5;
    bool roulette = use_rr && rng->getrand() < reflectance.max();
    if (is_max_depth || (use_rr && !roulette)) {
        return Vec::Zero;
    }

    Vec radiance = Vec::Zero;
    if (material->get_type()==1) {
        radiance = _photon_map->sample(pos, norm, global_n, global_r, caustic_n, caustic_r);
    }
    
        if (material->get_type()==3 && new_depth < 5) {

            Ray new_ray1, new_ray2;
            double pdf1, pdf2;
            material->sampling(ray, pos, norm, rng, new_ray1, pdf1);
            material->sampling(ray, pos, norm, rng, new_ray2, pdf2);
            radiance += trace(new_ray1, new_depth, rng, global_n, global_r, caustic_n, caustic_r) * pdf1;
            if (pdf2 > 1e-4) radiance += trace(new_ray2, new_depth, rng, global_n, global_r, caustic_n, caustic_r) * pdf2;
        }

        else {
            Ray new_ray; double pdf = 1.0;
           material->sampling(ray, pos, norm, rng, new_ray, pdf);
            radiance = trace(new_ray, new_depth, rng, global_n, global_r, caustic_n, caustic_r) / pdf;
        }

    return reflectance * radiance;
}



#endif //TRACER_FINAL_PPMTRACER_H

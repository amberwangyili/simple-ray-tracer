//
// Created by YILI WANG on 6/24/19.
//

#ifndef TRACER_FINAL_PHOTON_H
#define TRACER_FINAL_PHOTON_H
#include "../object/Object.h"
#include "../kdtree/PhotonKdtree.h"
#include "../Config.h"

struct Photon : Vec {
    Photon(const Vec &pos, const Vec &flux, const Vec &direction)
            : Vec(pos), flux(flux), direction(direction) {

    }

    Vec flux;
    Vec direction;
};

typedef Photon *ptr_photon_t;
typedef VecKDTree<Photon> PhotonKDTree;
typedef std::vector<ptr_photon_t> photon_vec_t;

enum class PhotonState : int {
    Direct    = 0x0100,
    Indirect  = 0x0101,
    Caustic   = 0x0200,
};

class PhotonMap {
public:
    struct PhotonDistanceCompare {
        PhotonDistanceCompare(const Vec &center) : center(center) {

        }

        inline bool operator ()(const ptr_photon_t &lhs, const ptr_photon_t &rhs) {
            return (*lhs - center).len2() < (*rhs - center).len2();
        }

        Vec center;
    };

    PhotonMap(Scene *scene = NULL, Scene *light = NULL) : _scene(scene), _light(light), _global(), _caustic(), _volume() {

    }

    void initialize();
    Vec sample(const Vec &position, const Vec &norm, int global_n, double global_r, int caustic_n, double caustic_r);

private:
    void find_knn(const Photon &center, int global_n, double global_r, int caustic_n, double caustic_r, photon_vec_t &global_result, photon_vec_t &caustic_result);
    void trace(const Ray &ray, const Vec &flux, int depth, PhotonState state, Random *rng);

    PhotonKDTree _global, _caustic, _volume;
    Scene *_scene, *_light;
};


void PhotonMap::initialize() {
    bool finished_global = false;
    bool finished_caustic = false;

    Random *rng = new Random(19971226 + rand());
    for (int i = 0; i < PHOTON_SAMPLE; ++i) {
        for (Object *object : _light->objects) {

            Ray light;
            double pdf;
            object->sample(rng, light, pdf);

            Vec flux = object->light->get_emission() * pi / pdf / PHOTON_SAMPLE;
            light.direction = rng->sample_hemisphere(light.direction);

            trace(light, flux, 0, PhotonState::Direct, rng);
        }
    }
    delete rng;

//    std::cerr << "global: " << _global.wrapper.size() << std::endl;
//    std::cerr << "caustic: " << _caustic.wrapper.size() << std::endl;
    _global.initialize();
    _caustic.initialize();
}

Vec PhotonMap::sample(const Vec &position, const Vec &norm, int global_n, double global_r, int caustic_n, double caustic_r) {
    Photon query(position, Vec::Zero, Vec::Zero);
    photon_vec_t global_result, caustic_result;
    std::vector<double> distances;
    find_knn(query, global_n, global_r, caustic_n, caustic_r, global_result, caustic_result);

    Vec flux = Vec::Zero;

    photon_vec_t valid;

    double max_dist = 0.0;
    for (Photon *photon : global_result) {
        Vec diff = position - *photon;
        double dist = diff.len();
        double dt   = dot(norm, diff) / dist;
        if (abs(dt) < global_r * global_r * 0.01) {
            valid.push_back(photon);
            max_dist = max(max_dist, dist);
        }
    }

    if (max_dist > eps) {
        double k = 1.1;

        for (size_t i = 0; i < valid.size(); ++i) {
            double dist = (position - *(valid[i])).len();
            double w = 1.0 - (dist / (k * max_dist));
            Vec v = valid[i]->flux * inv_pi;
            flux += w * v;
        }
        flux = flux / (1.0 - 2.0 / (3.0 * k));

        return flux / (pi * square(max_dist));
    }

    return Vec::Zero;
}

void PhotonMap::find_knn(const Photon &center, int global_n, double global_r, int caustic_n, double caustic_r,
                         photon_vec_t &global_result, photon_vec_t &caustic_result) {

    photon_vec_t result1 = _global.find_knn(center, global_r, global_n);

    for (int i = 0; i < min(global_n, static_cast<int>(result1.size())); ++i) {
        result1[i]->flux = result1[i]->flux * PHOTON_GLOBAL_MUL;
        global_result.push_back(result1[i]);
    }
}

void PhotonMap::trace(const Ray &ray, const Vec &flux, int depth, PhotonState state, Random *rng) {
    if (depth > MAX_DEPTH) return ;
    if (flux.max() <= 0) return ;

    Hitpoint inter = _scene->intersect(ray);
    if (inter.object == NULL) return ;
    if (inter.object->is_light) return ;

    int new_depth = depth + 1;

    Object *object = inter.object;
    Material *material = object->material;
    Vec pos = inter.position;
    Vec norm = inter.norm;


    Vec reflectance = material->get_reflectance(pos, norm);

    double photon_pdf = 1.0;
    if (material->get_type()==1) {
        if (_global.wrapper.size() < PHOTON_GLOBAL)
            _global.wrapper.push_back(new Photon(pos, flux, ray.direction));

        const double prob = reflectance.mean();
        if (rng->getrand() < prob) {
            photon_pdf *= prob;
        } else {
            return ;
        }
    }

    PhotonState new_state = state;

    Ray new_ray; double sample_pdf;
    material->sampling(ray, pos, norm, rng, new_ray, sample_pdf);
    Vec new_flux = flux * reflectance / (sample_pdf * photon_pdf);

    trace(new_ray, new_flux, new_depth, new_state, rng);
}
#endif //TRACER_FINAL_PHOTON_H

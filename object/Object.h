//
// Created by YILI WANG on 6/24/19.
//

#ifndef TRACER_FINAL_OBJECT_H
#define TRACER_FINAL_OBJECT_H

#include "../base/Common.h"
#include "../base/Ray.h"
#include "Material.h"
class KdTree;



class Hitpoint;
//成员：交点所在物体object，
//     距离ray出射点的相对距离distance
//     绝对坐标position
//     交点在对应面上的法向量

class Object;
//成员: 材质类material或光源体 bool值（是否为光源）
//操作：intersect() 虚函数：求交,输入Ray，返回Hitpoint类
//     set_material() 初始化：设置材质material，且设置光源bool值为false
//     set_light()初始化：设置光源体light，且设置光源bool值为true

class Scene;
//成员：vector<Object*> 场景中存在的物体
//操作：add() 加入物体
//     intersect() 光线与场景中相交的第一个物体

class AABB;
//成员：vmin vmax（axis-aligned）
//操作：intersect() 继承：求交,输入Ray，返回Hitpoint类


class Plane;
//成员：中心点point，法向量norm，上下绝对长度d（两倍）
//     intersect()求交

class Sphere;
//成员：半径，圆心，半径平方
//操作：intersect() 求交
//     sample() 射入一条射线，和对应的pdf，返回反射的随机一条射线

class Tri;
//成员：三个顶点，法向量，反法向量，重心
//操作：area（）算面积
//     intersect() 求交
//     sample() 射入一条射线，和对应的pdf，返回反射的随机一条射线 (sample in triangle)

class Quad;
//成员：两个三角形拼接而成的四边形tri1，tri2, tri1,tri2的面积
//操作：intersect()求交
//     sample() 在四边形内sample出一条射线

class TriMesh;
//成员：友元类 Kd树
//     vector(tri)
//     AABB bounding box
//操作：intersect()求交
//     sample() 在四边形内sample出一条射线

typedef std::vector<Tri *> TriVec;
typedef std::vector<Vec> VertexVec;

using namespace std;

class Hitpoint {
public:
    Object *object;
    double distance;
    Vec position;
    Vec norm;

    Hitpoint() : object(NULL), distance(0) { }
    const static Hitpoint null;
};


class Scene {
public:
    vector<Object*> objects;
    Scene(void) { }

public:
    virtual Scene *add(Object *obj){objects.push_back(obj);return this;}
    Hitpoint intersect(const Ray &ray);
};




class Object {
public:
    union {
        class Material *material;
        Light *light;
    };
    bool is_light;

    Object() : material(NULL), is_light(false) {}
    virtual Hitpoint intersect(const Ray &ray) { return Hitpoint::null; }
    virtual void sample(Random *rng, Ray &ray, double &pdf) = 0;
    inline Object *set_material(class Material *const m) {
        material = m, is_light = false;
        return this;
    }
    inline Object *set_light(Light *const l) {
        light = l, is_light = true;
        return this;
    }
};


class AABB {
public:
    Vec vmin, vmax;

    AABB(void) { }

    AABB(Vec vmin, Vec vmax) : vmin(vmin), vmax(vmax) {
        initialize();
    }

    inline void initialize(void) {

    }

    inline Pair get_near_far(const Ray &ray) {
        Vec result_min = vmin - ray.origin;
        Vec result_max = vmax - ray.origin;
        result_min = result_min / (ray.direction + eps);
        result_max = result_max / (ray.direction + eps);
        for (int i = 0; i < 3; ++i)
            if (result_min[i] > result_max[i])
                std::swap(result_min[i], result_max[i]);
        return Pair(result_min.max(), result_max.min());
    }

    inline bool intersect(const Ray &ray) {
        Pair near_far = get_near_far(ray);
        return !(near_far.first > near_far.second || near_far.second < 0);
    }
};



class Sphere : public Object {
public:
    Vec center;
    double radius;
    double sqrrad;

    Sphere(void) { }

    Sphere(const Vec &o, double r) : center(o), radius(r) {
        sqrrad = r * r;
    }
    virtual Hitpoint intersect(const Ray &ray);
    virtual void sample(Random *rng, Ray &ray, double &pdf){
        ray = rng->sample_sphere(center, radius);
        pdf = 1;
    }
    Hitpoint intersect(const Ray &ray, double t_min, double t_max){
        Vec oc = center - ray.origin;
        double b = dot(oc, ray.direction);
        double det = b * b - oc.len2() + sqrrad;
        if (det > 0) {
            double sdet = sqrt(det);
            double temp  = b-sdet;
            if(temp<t_max&&temp>t_min){
                Hitpoint res;
                res.distance  = temp;
                res.object = static_cast<Object*>(this);
                res.position = ray.get(temp);
                res.norm = (res.position-center).normalize();
                return  res;
            }
            temp = b+sdet;
            if(temp<t_max&&temp>t_min){
                Hitpoint res;
                res.distance  = temp;
                res.object = static_cast<Object*>(this);
                res.position = ray.get(temp);
                res.norm = (res.position-center).normalize();
                return  res;
            }
        }
        return Hitpoint::null;
    }
};



class Plane : public Object {
public:
    Vec point, norm;
    double d;
    Texture *texture;

    Plane(void) : texture(NULL) { }

    Plane(const Vec &p, const Vec &n) : point(p), norm(n), texture(NULL) {
        initialize();
    }

    inline void initialize(void) {
        d = -dot(point, norm);
    }


    virtual Hitpoint intersect(const Ray &ray);
    virtual void sample(Random *rng, Ray &ray, double &pdf) {}
};

class Tri : public Object {
public:
    Vec a, b, c, norm, nnorm, center;
    double d;

    Tri(void) { }

    Tri(const Vec &a_, const Vec &b_, const Vec &c_) : a(a_), b(b_), c(c_) {
        initialize();
    }


    inline void initialize(void) {
        norm = cross(b - a, c - a).normalize();
        nnorm = -norm;
        center = (a + b + c) / 3;
        d = -dot(norm, a);
    }

    inline double area() {
        return abs(cross(c - a, b - a).len()) / 2;
    }

    virtual Hitpoint intersect(const Ray &ray);
    virtual void sample(Random *rng, Ray &ray, double &pdf){
        ray = rng->sample_triangle(a, b, c, norm);
        pdf = 1;
    }
};


class Quad : public Object {
public:
    Tri tri1, tri2;
    double area1, area2;
    Quad(void) { }

    Quad(const Vec &a, const Vec &b, const Vec &c, const Vec &d) {
        tri1 = Tri(a, b, c);
        tri2 = Tri(a, c, d);
        initialize();
    }
    inline void initialize(void) {
        area1 = tri1.area();
        area2 = tri2.area();
    }


    virtual Hitpoint intersect(const Ray &ray) {
        Hitpoint inter = tri1.intersect(ray);
        if (inter.object != NULL) {
            inter.object = static_cast<Object *>(this);
            return inter;
        }
        inter = tri2.intersect(ray);
        if (inter.object != NULL)
            inter.object = static_cast<Object *>(this);
        return inter;
    }
    virtual void sample(Random *rng, Ray &ray, double &pdf) {
        if (rng->getrand() < area1 / (area1 + area2)) tri1.sample(rng, ray, pdf);
        else tri2.sample(rng, ray, pdf);
        pdf = 1 / (area1 + area2);
    }
};



class TriMesh: public Object{
    friend class KdTree;
public:
    TriMesh(void) { }

    TriMesh(const TriVec &triangles, const AABB &bbox)
            : _triangles(triangles), _bbox(bbox) {

    }

    TriMesh(const TriVec &triangles, const VertexVec &vertexes)
            : _triangles(triangles) {

        initialize(vertexes);
    }
    TriMesh(const string &infile_name, Material *material, const Vec &resize,const Vec &delta,bool flip_normal, double rotate,Vec axis);

    inline void initialize(VertexVec vertexes) {
        Vec vmin(inf, inf, inf), vmax(-inf, -inf, -inf);
        for (auto v : vertexes) {
            if (v.x < vmin.x) vmin.x = v.x;
            if (v.y < vmin.y) vmin.y = v.y;
            if (v.z < vmin.z) vmin.z = v.z;

            if (v.x > vmax.x) vmax.x = v.x;
            if (v.y > vmax.y) vmax.y = v.y;
            if (v.z > vmax.z) vmax.z = v.z;
        }
        _bbox = AABB(vmin, vmax);
    }

    virtual Hitpoint intersect(const Ray &ray);
    virtual void sample(Random *rng, Ray &ray, double &pdf){
        Hitpoint pt = this->intersect(ray);
        Vec dir = rng->sample_hemisphere(pt.norm);
        ray = Ray(pt.position,dir);
        pdf = 1;
    }


public:
    TriVec _triangles;
    AABB    _bbox;

};


TriMesh::TriMesh(const string &infile_name, Material *material, const Vec &resize,
                 const Vec &delta, bool flip_normal,double rotate,Vec axis) {
    FILE *inobj = fopen(infile_name.c_str(), "r");
    char *buffer = new char[1024];
    VertexVec vertexes;
    std::vector<int> fx, fy, fz;
    fx.clear(), fy.clear(), fz.clear();
    while (!feof(inobj)) {
        fgets(buffer, 1024, inobj);
        if (buffer[0] == '#') {
            continue;
        } else if (buffer[0] == 'v') {
            char op;
            double x, y, z;
            sscanf(buffer, "%c %lf %lf %lf", &op, &x, &y, &z);
            Vec vertex = (Vec(x, y, z) * resize) + delta;
            vertex = vertex.rotate(axis,2*pi*rotate);
            vertexes.push_back(vertex);
        } else if (buffer[0] == 'f') {
            char op;
            int x, y, z;
            sscanf(buffer, "%c %d %d %d", &op, &x, &y, &z);
            if (x == y || y == z || x == z)
                continue;
            fx.push_back(x), fy.push_back(y), fz.push_back(z);
        }
    }

    for (size_t i = 0; i < fx.size(); ++i) {
        int x = fx[i], y = fy[i], z = fz[i];
        if (flip_normal) std::swap(x, y);
        Tri *triangle = new Tri(vertexes[x - 1], vertexes[y - 1], vertexes[z - 1]);
        triangle->set_material(material);
        _triangles.push_back(triangle);
    }

    initialize(vertexes);
}





class Constant_medium: public Object{
public:
    Sphere *boundary;
    float density;
    Constant_medium(Sphere *b, float d):boundary(b),density(d){}
    virtual Hitpoint intersect(const Ray &ray);
    virtual void sample(Random *rng, Ray &ray, double &pdf) {}

};


const Hitpoint Hitpoint::null = Hitpoint();

Hitpoint Constant_medium::intersect(const Ray &ray) {
    double db = (drand48()<1e-5);
    db = false;
    Hitpoint rec1, rec2;
    rec1=boundary->intersect(ray,-FLT_MAX,FLT_MAX);
    rec2 = boundary->intersect(ray,rec1.distance+1e-4,FLT_MAX);

    if(rec1.object!=NULL&&rec2.object!=NULL){
        if(db)std::cerr<<"\n t0 t1 " << rec1.distance<<rec2.distance<<endl;
            if(rec1.distance>=rec2.distance) return Hitpoint::null;
            if(rec1.distance<0) rec1.distance = 0;
            double distance_inside_boundary = rec2.distance-rec1.distance;
            double hit_distance = -(1/density)*log(drand48());
            if (hit_distance<distance_inside_boundary){
                Hitpoint res;
                res.distance = rec1.distance+hit_distance;
                res.position = ray.get(res.distance);
                res.norm = Vec(1,0,0);
                res.object = static_cast<Object*>(this);
                if(db)std::cerr<<res.distance<<endl;
                if(db)std::cerr<<res.position<<endl;
                return res;
            }
    }
}

Hitpoint Sphere::intersect(const Ray &ray) {
    Vec oc = center - ray.origin;
    double b = dot(oc, ray.direction);
    double det = b * b - oc.len2() + sqrrad;
    if (det > 0) {
        double sdet = sqrt(det);
        double distance = 0;
        if (b - sdet > eps)
            distance = b - sdet;
        else if (b + sdet > eps)
            distance = b + sdet;
        if (distance > eps) {
            Hitpoint res;
            res.object = static_cast<Object *>(this);
            res.distance = distance;
            res.position = ray.get(distance);
            res.norm = (res.position - center).normalize();
            return res;
        }
    }

    return Hitpoint::null;
}


Hitpoint Plane::intersect(const Ray &ray) {
    //check parallel with plane
    double a = - dot(norm, ray.direction);
    if (a <= eps)
        return Hitpoint::null;
    //compute hitpoint
    double b = dot(norm, ray.origin - point);
    Hitpoint res;
    res.object = static_cast<Object *>(this);
    res.distance = b / a;
    res.position = ray.get(res.distance);
    res.norm = norm;
    return res;
}


Hitpoint Scene::intersect(const Ray &ray) {
    Hitpoint res = Hitpoint::null;
    for (Object *obj : objects) {
        Hitpoint tmp = obj->intersect(ray);
        if (tmp.object != NULL) {
            if (res.object == NULL || res.distance > tmp.distance) {
                res = tmp;
            }
        }
    }
    return res;
}

inline static bool checkside(const Vec &a, const Vec &b, const Vec &c, const Vec &p) {
    Vec ab = b - a, ac = c - a, ap = p - a;
    Vec v1 = cross(ab, ac);
    Vec v2 = cross(ab, ap);
    return dot(v1, v2) >= 0;
}

inline static bool in_tri(const Vec &p, const Tri &t) {
    return checkside(t.a, t.b, t.c, p) && checkside(t.b, t.c, t.a, p) && checkside(t.c, t.a, t.b, p);
}

Hitpoint Tri::intersect(const Ray &ray) {
    double t = -(d + dot(norm, ray.origin)) / (dot(norm, ray.direction) + eps);
    if (t > eps) {
        Vec point = ray.get(t);
        bool inter = in_tri(point, *this);
        if (inter) {
            Hitpoint res;
            res.object = static_cast<Object *>(this);
            res.distance = t;
            res.position = point;
            res.norm = norm;
            return res;
        }
    }
    return Hitpoint::null;
}

Hitpoint TriMesh::intersect(const Ray &ray) {
    if (!_bbox.intersect(ray)) {
        return Hitpoint::null;
    }

    Hitpoint res = Hitpoint::null;
    for (Tri *triangle : _triangles) {
        Hitpoint tmp = triangle->intersect(ray);
        if (tmp.object != NULL) {
            if (res.object == NULL || res.distance > tmp.distance)
                res = tmp;
        }
    }

    return res;
}



#endif //TRACER_FINAL_OBJECT_H

//
// Created by YILI WANG on 6/24/19.
//

#ifndef TRACER_FINAL_KDTREE_H
#define TRACER_FINAL_KDTREE_H

#include "../object/Object.h"
using namespace std;



template <typename NT, typename VT>

struct Node{
    typedef NT node_t;
    typedef VT value_t;

    Node(void) : axis(-1), split(0), lson(NULL), rson(NULL) { }
    inline bool is_leaf(void) const { return lson == NULL && rson == NULL; }

    int     axis;
    value_t split;
    node_t  *lson, *rson;
};

template <typename PT>
struct Compare{
    typedef PT ptr_t;
    virtual bool operator ()(const ptr_t &lhs, const ptr_t &rhs) const = 0;
};



class Kdtree: public Object{
public:
    struct KdNode:Node<KdNode,double>{
        KdNode():Node<KdNode,double>(){}
        AABB aabb;
        TriMesh mesh;
    };

    Kdtree(TriMesh *wrap, int max_leaf = 64):wrap(wrap),root(NULL),vmax_leaf(max_leaf){
        initialize();
    }
    void initialize();
    virtual Hitpoint intersect(const Ray &ray);
    virtual void sample(Random *rng, Ray &ray, double &pdf) {}
    void build(KdNode *&root, const std::vector<Tri *> &a, AABB aabb, int current);
    Hitpoint traverse(KdNode *root,const Ray &ray);
public:
    TriMesh *wrap;
    KdNode *root;
    int vmax_leaf = 64;

};

inline static int get_split_side(const Tri *triangle, int axis, double split) {
    int count[2] = {0, 0};
    split -= eps;
    ++count[triangle->a[axis] < split];
    ++count[triangle->b[axis] < split];
    ++count[triangle->c[axis] < split];
    if (count[1] == 3)
        return -1;
    if (count[0] == 3)
        return 1;
    return 0;
}

void Kdtree::build(KdNode *&root, const std::vector<Tri *> &a, AABB aabb, int current) {
    if (a.empty()) {
        return;
    }

    root = new KdNode();
    root->aabb = aabb;

    size_t n = a.size();
    if (n <= vmax_leaf) {
        root->mesh = TriMesh(a, aabb);
    } else {
        root->axis = current;
        root->split = (aabb.vmax[current] + aabb.vmin[current]) / 2;

        std::vector<Tri *> la, ra;
        for (Tri *triangle : a) {
            int res = get_split_side(triangle, root->axis, root->split);
            if (res <= 0) la.push_back(triangle);
            if (res >= 0) ra.push_back(triangle);
        }

        Vec lvmin = aabb.vmin, lvmax = aabb.vmax;
        Vec rvmin = aabb.vmin, rvmax = aabb.vmax;
        lvmax[current] = root->split, rvmin[current] = root->split;

        int next = (current + 1) % 3;
        build(root->lson, la, AABB(lvmin, lvmax), next);
        build(root->rson, ra, AABB(rvmin, rvmax), next);
    }
}


Hitpoint Kdtree::traverse(KdNode *root, const Ray &ray) {
    if (root == NULL)
        return Hitpoint::null;

    if (root->is_leaf()) {
        return root->mesh.intersect(ray);
    } else {
        KdNode *lson = root->lson;
        KdNode *rson = root->rson;

        int axis = root->axis;
        Pair p = root->aabb.get_near_far(ray);
        double near = p.first, far = p.second;

        if (near > far || far < -1)
            return Hitpoint::null;

        double split = root->split;
        double near_x = ray.get(near)[axis];
        double far_x = ray.get(far)[axis];
        bool ignored = false, swaped = false;

        if (near_x > split) {
            std::swap(lson, rson), swaped = true;
            if (far_x > split) ignored = true;
        } else if (far_x < split) {
            ignored = true;
        }

        Hitpoint interl = traverse(lson, ray);
        if (interl.object != NULL) {
            if (ignored)
                return interl;
            if (!swaped && interl.position[axis] < split + 1e-9)
                return interl;
            if (swaped && interl.position[axis] > split - 1e-9)
                return interl;
        }

        if (!ignored) {
            Hitpoint interr = traverse(rson, ray);
            if (interr.distance < far + 1e-9)
                return interr;
        }
    }

    return Hitpoint::null;
}
void Kdtree::initialize() {
    build(root, wrap->_triangles, wrap->_bbox, 0);
}

Hitpoint Kdtree::intersect(const Ray &ray) {
    Hitpoint res = traverse(root, ray);
    return res;
}




#endif //TRACER_FINAL_KDTREE_H

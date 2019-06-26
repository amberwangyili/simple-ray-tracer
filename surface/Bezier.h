//
// Created by YILI WANG on 6/24/19.
//

#ifndef TRACER_FINAL_BEZIER_H
#define TRACER_FINAL_BEZIER_H

#include "../object/Object.h"
#include <cmath>
#include <ctime>

#include "../object/Object.h"
#include <cmath>
#include <ctime>
#define tol 0.1
using namespace std;
//
//
//class Vec2{
//public:
//    union {
//        struct {
//            double x, y;
//        };
//        struct {
//            double p[2];
//        };
//    };
//    Vec2(double x_=0, double y_=0){ x=x_; y=y_;}
//
//    inline double &operator[](int index) {
//        return p[index];
//    }
//
//    inline const double &operator[](int index) const {
//        return p[index];
//    }
//};
//Vec2 operator*(const Vec2 &a, const double &b) {
//    return Vec2(a.x * b, a.y * b);
//}
//
//Vec2 operator*( const double &b,const Vec2 &a){
//    return Vec2(a.x * b, a.y * b);
//}
//Vec2 operator+(const Vec2 &a, const Vec2 &b) {
//    return Vec2(a.x + b.x, a.y + b.y);
//}
//
//Vec2 operator-(const Vec2 &a, const Vec2 &b) {
//    return Vec2(a.x - b.x, a.y - b.y);
//}
//
//double sphere_hit(Ray r, Vec o, double radius) {
//    Vec ro = o - r.origin;
//    double b = r.direction.dot(ro);
//    double d = b*b- ro.dot(ro) + radius*radius;
//    if (d < 0) return -1;
//    else d = sqrt(d);
//    double t = b - d > eps ? b - d : b + d > eps? b + d : -1;
//    if (t < 0)
//        return -1;
//    return t;
//}
//class BezierCurve {
//public:
//    vector<Vec2> control_points;
//    double height;
//    int n;
//    double *dx,*dy;
//
//    BezierCurve(vector<Vec2> c_pts): control_points(c_pts){
//        n = control_points.size();
//        --n;
//        dx = new double[n];
//        dy = new double[n];
//
//        for (int i = 0; i <= n; ++i) {
//            dx[i] = c_pts[0].x;
//            dy[i] = c_pts[0].y;
//            for (int j = 0; j <= n - i; ++j) {
//                c_pts[j].x = c_pts[j+1].x-c_pts[j].x;
//                c_pts[j].y = c_pts[j+1].y-c_pts[j].y;
//            }
//        }
//        double n_down = 1, fac=1,nxt = n;
//        for (int i = 0; i <=n; ++i,--nxt) {
//            fac = fac*(i==0?1:i);
//            dx[i] = dx[i]*n_down/fac;
//            dy[i] = dy[i]*n_down/fac;
//            n_down*=nxt;
//        }
//        height = eval(1)[1];
//    }
//
//    double compute_width();
////    static Vec2 eval(const BezierCurve &curve, double t);
//    Vec2 eval(double t){
//        double ans_x = 0, ans_y = 0, t_pow = 1;
//        for (int i = 0; i <=n ; ++i) {
//            ans_x += dx[i]*t_pow;
//            ans_y += dy[i]*t_pow;
//            t_pow *= t;
//        }
//        return Vec2(ans_x,ans_y);
//    }
//
////    static Vec2 deri(const BezierCurve &curve, double t);
//    Vec2 deri(double t){
//        double ans_x = 0, ans_y = 0, t_pow = 1;
//        for (int i = 1; i <=n; ++i) {
//            ans_x += dx[i]*i*t_pow;
//            ans_y += dy[i]*i*t_pow;
//            t_pow *= t;
//        }
//        return Vec2(ans_x,ans_y);
//    }
//};
//double BezierCurve::compute_width() {
//    double temp = -inf;
//    for (double i = 0; i < 1; i+=0.001) {
//        if (abs(eval(i)[0])>temp) temp = abs(eval(i)[0]);
//    }
//    return temp+tol;
//}
////Vec2 BezierCurve::eval(const BezierCurve &curve, double t) {
////    // bezier curve evaluation using de Casteljau's algo
////    if (curve.control_points.size() == 1)
////        return curve.control_points[0];
////    int size = curve.control_points.size();
////    vector<Vec2> new_points(size - 1);
////    for (int i = 0; i < size - 1; i++)
////        new_points[i] = curve.control_points[i] * (1 - t) + curve.control_points[i + 1] * t;
////    BezierCurve reduced = BezierCurve(new_points);
////    return eval(reduced, t);
////}
////
////
////
////Vec2 BezierCurve::deri(const BezierCurve &curve, double t) {
////    int size = curve.control_points.size();
////    vector<Vec2> new_points(size - 1);
////    for (int i = 0; i < size - 1 ; i++)
////        new_points[i] = (curve.control_points[i + 1] - curve.control_points[i]) * (size - 1);
////    BezierCurve reduced = BezierCurve(new_points);
////    return eval(reduced, t);
////}
//
//class BezierSurface: public Object{
//public:
//    BezierCurve curve;
//    Vec pos;
//    AABB box;
//    BezierSurface(const Vec pos_, const BezierCurve c_):pos(pos_),curve(c_){
//        Vec vmin  = Vec(pos.x-curve.compute_width()-tol,pos.y-tol,pos.z-curve.compute_width()-tol);
//        Vec vmax =  Vec(pos.x+curve.compute_width()+tol,pos.y+curve.height+tol,pos.z+curve.compute_width()+tol);
//        box = AABB(vmin,vmax);
//    }
//
//    double func(double t, const Ray &ray){
//        double xt, yt;
//        xt = curve.eval(t).x;
//        yt = curve.eval(t).y;
//
//        double x,z;
//        x = ray.origin.x+ray.direction.x/ray.direction.y*(yt-ray.origin.y);
//        z = ray.origin.z+ray.direction.z/ray.direction.y*(yt-ray.origin.y);
//        double res = xt*xt-x*x-z*z;
//        return res;
//    }
//
//    double func_deri(double t, const Ray &ray){
//
//        double xt, yt;
//        xt = curve.eval(t).x;
//        yt = curve.eval(t).y;
//        double x,z;
//        x = ray.origin.x+ray.direction.x/ray.direction.y*(yt-ray.origin.y);
//        z = ray.origin.z+ray.direction.z/ray.direction.y*(yt-ray.origin.y);
//
//        double deri_xt, deri_yt;
//        deri_xt = curve.deri(t).x;
//        deri_yt = curve.deri(t).y;
//        double res = 2*xt*deri_xt-2*x*deri_yt*ray.direction.x/ray.direction.y-2*deri_yt*z*ray.direction.z/ray.direction.y;
//        return res;
//    }
//
//
//
//    //ref: trickle for dealing with degenerate case;
//    Vec get_norm(Vec p);
//    Vec2 get_theta(Vec hit);
//    double get_curve_t(double ypos);
//    Hitpoint degenerate_case(const Ray &ray);
//    double newton(double t, const Ray &ray);
//    virtual Hitpoint intersect(const Ray &ray);
//
//
//};
//
//double BezierSurface::newton(double t, const Ray &ray){
//    for (int i = 0; i < 10; ++i) {
//        double ft =func(t,ray);
//        double dt = func_deri(t,ray);
//        t -= ft/dt;
//        if(abs(ft)<5e-3){return t;}
//    }
//    return -1;
//}
//
//
//Hitpoint BezierSurface::intersect(const Ray &ray) {
////    static int cnt = 0;
////    static double t = 0;
////    if (cnt == 0) t = clock();
////    ++cnt;
////    if (cnt == 10000)
////    {
////        printf("%lf\n", (clock() - t) / CLOCKS_PER_SEC);
////        exit(0);
////    }
//    //inline compute aabbox:
//    if(!box.intersect(ray)) {return  Hitpoint::null;}
//    if(abs(ray.direction.y)>eps){
//        double ray_t =inf;
//        double curve_t = -1;
//        double temp;
//        for (double t = 0; t <=1 ; t+=0.1) {
//            temp = newton(t,ray);
//            if(temp<1-eps&&temp>0+eps){
//                double ray_temp = (curve.eval(temp).y-ray.origin.y)/ray.direction.y;
//                if(ray_temp>0&&ray_temp<ray_t){
//                    ray_t = ray_temp;
//                    curve_t = temp;
//                }
//            }
//        }
//        if(curve_t<1-eps&&curve_t>0+eps&&func(curve_t,ray)<5e-4){
//            Hitpoint res;
//            res.object = static_cast<Object*>(this);
//            res.distance = ray_t;
//            res.position = ray.get(ray_t);
////            cout<<"intersect!!!!!"<<endl;
//
//            res.norm = get_norm(res.position).normalize();
//            return res;
//        }
//
//
//    }
//    else{
//        Hitpoint hit;
//        hit = degenerate_case(ray);
//        return hit;
//    }
//    return  Hitpoint::null;
//}
//
//Vec BezierSurface::get_norm(Vec p){
//    Vec2 tmp = get_theta(p);
//    Vec2 dir = curve.deri(tmp.y);
//    Vec d_surface  = Vec(cos(tmp.x), dir.y / dir.x, sin(tmp.x));
//    Vec d_circ = Vec(-sin(tmp.x), 0, cos(tmp.x));
//    return cross(d_surface,d_circ).normalize();
//}
//
//Vec2 BezierSurface::get_theta(Vec hit) {
//    double t  = get_curve_t(hit.y-pos.y);
//    double u = atan2(hit.z-pos.z,hit.x-pos.x);
//    if(u<0) u += 2*pi;
//    return Vec2(u,t);
//
//}
//
//Hitpoint BezierSurface::degenerate_case(const Ray &ray)
//{
//    double final_dis = inf;
//    double dis_to_axis = (Vec(pos.x, ray.origin.y, pos.z) - ray.origin).len();
//    double hit = ray.get(dis_to_axis).y;
//    if (hit < pos.y + eps || hit > pos.y + curve.height - eps)
//        return Hitpoint::null;
//    // solve function pos.y+y(t)=ray.origin.y to get x(t)
//    double t = get_curve_t(hit - pos.y);
//    if (t < 0+eps || t > 1-eps)
//        return Hitpoint::null;
//    Vec2 loc = curve.eval(t);
//    double ft = pos.y + loc.y - hit;
//    if (std::abs(ft) > eps)
//        return Hitpoint::null;
//    // assume sphere (pos.x, pos.y + loc.y, pos.z) - loc.x
//    final_dis = sphere_hit(ray, Vec(pos.x, pos.y + loc.y, pos.z), loc.x);
//    if (final_dis < 0)
//        return Hitpoint::null;
//    Vec inter_p = ray.get(final_dis);
//    // printf("y %f small!!!",std::abs((inter_p - Vec(pos.x, inter_p.y, pos.z)).len2() - sqr(loc.x)));
//    if (std::abs((inter_p - Vec(pos.x, inter_p.y, pos.z)).len2() - (loc.x)*(loc.x)) > 5e-4)
//        return Hitpoint::null;
//    // second iteration, more accuracy
//    hit = inter_p.y;
//    if (hit < pos.y + eps || hit > pos.y + curve.height - eps)
//        return Hitpoint::null;
//    t = get_curve_t(hit - pos.y);
//    loc = curve.eval(t);
//    ft = pos.y + loc.y - hit;
//    if (std::abs(ft) > eps)
//        return Hitpoint::null;
//    final_dis = sphere_hit(ray, Vec(pos.x, pos.y + loc.y, pos.z), loc.x);
//    if (final_dis < 0)
//        return Hitpoint::null;
//    inter_p = ray.get(final_dis);
//    if (std::abs((inter_p - Vec(pos.x, hit, pos.z)).len2() - (loc.x)*(loc.x)) > 5e-4)
//        return Hitpoint::null;
//    // printf("---y %f small!!!",std::abs((inter_p - Vec(pos.x, inter_p.y, pos.z)).len2() - sqr(loc.x)));
//    Hitpoint res;
//    res.distance = final_dis;
//    res.position = inter_p;
//    res.norm = get_norm(inter_p);
//    res.object = static_cast<Object*>(this);
//    return res;
//}
//
//double BezierSurface::get_curve_t(double ypos){
//    double deriv,func;
//    double t = .5;
//    for (int i = 0; i < 20; ++i) {
//        if (t<0) t = 0 ;
//        else if (t>1) t = 1;
//        func = curve.eval(t).y - ypos;
//        deriv = curve.deri(t).y;
//        if(abs(func)<eps) return t;
//        t = t-func/deriv;
//    }
//    return -1;
//}
//



class BezierCurve
{// f(y)=x, y goes up?
public:
    double *dx, *dy, max, height, max2, r, num;
    int n;
    struct D{
        double t0, t1, width, y0, y1, width2;
    }data[20];
    // x(t) = \sum_{i=0}^n dx_i * t^i
    // y(t) = \sum_{i=0}^n dy_i * t^i
    BezierCurve(double* px, double* py, int n_, int num_, double r_): num(num_), n(n_), r(r_) {
        dx = new double[n];
        dy = new double[n];
        --n;
        // preproces
        for(int i = 0; i <= n; ++i)
        {
            dx[i] = px[0];
            dy[i] = py[0];
            for (int j = 0; j <= n - i; ++j)
            {
                px[j] = px[j + 1] - px[j];
                py[j] = py[j + 1] - py[j];
            }
        }
        double n_down = 1, fac = 1, nxt = n;
        for (int i = 0; i <= n; ++i, --nxt)
        {
            fac = fac * (i == 0 ? 1 : i);
            dx[i] = dx[i] * n_down / fac;
            dy[i] = dy[i] * n_down / fac;
            n_down *= nxt;
        }
        max = 0;
        double interval = 1. / (num - 1), c = 0;
        for (int cnt = 0; cnt <= num; c += interval, ++cnt)
        {
            data[cnt].width = 0;
            data[cnt].t0 = std::max(0., c - r);
            data[cnt].t1 = std::min(1., c + r);
            data[cnt].y0 = getpos(data[cnt].t0).y;
            data[cnt].y1 = getpos(data[cnt].t1).y;
            for (double t = data[cnt].t0; t <= data[cnt].t1; t += 0.00001)
            {
                Vec pos = getpos(t);
                if (data[cnt].width < pos.x)
                    data[cnt].width = pos.x;
            }
            if (max < data[cnt].width)
                max = data[cnt].width;
            data[cnt].width += eps;
            data[cnt].width2 = square(data[cnt].width);
        }
        max += eps;
        max2 = max * max;
        height = getpos(1).y;
    }
    Vec getpos(double t)
    {
        double ans_x = 0, ans_y = 0, t_pow = 1;
        for (int i = 0; i <= n; ++i)
        {
            ans_x += dx[i] * t_pow;
            ans_y += dy[i] * t_pow;
            t_pow *= t;
        }
        return Vec(ans_x, ans_y);
    }
    Vec getdir(double t)
    {
        double ans_x = 0, ans_y = 0, t_pow = 1;
        for(int i = 1; i <= n; ++i)
        {
            ans_x += dx[i] * i * t_pow;
            ans_y += dy[i] * i * t_pow;
            t_pow *= t;
        }
        return Vec(ans_x, ans_y);
    }
    Vec getdir2(double t)
    {
        double ans_x = 0, ans_y = 0, t_pow = 1;
        for(int i = 2; i <= n; ++i)
        {
            ans_x += dx[i] * i * (i - 1) * t_pow;
            ans_y += dy[i] * i * (i - 1) * t_pow;
            t_pow *= t;
        }
        return Vec(ans_x, ans_y);
    }
};


class BezierSurface: public Object {
// the curve will rotate line (x=pos.x and z=pos.z) as pivot
public:
    BezierCurve curve;
    Vec pos; // the buttom center point
    BezierSurface(Vec pos_, BezierCurve c_):
            pos(pos_), curve(c_){}

    double solve_t(double yc) { // solve y(t)=yc
        // assert(0 <= yc && yc <= curve.height);
        double t = .5, ft, dft;
        for (int i = 10; i--; )
        {
            if (t < 0) t = 0;
            else if (t > 1) t = 1;
            ft = curve.getpos(t).y - yc, dft = curve.getdir(t).y;
            if (std::abs(ft) < eps)
                return t;
            t -= ft / dft;
        }
        return -1;
    }
    virtual void sample(Random *rng, Ray &ray, double &pdf) {}
    virtual Vec change_for_bezier(Vec inter_p) {
        double t = solve_t(inter_p.y - pos.y);
        double u = atan2(inter_p.z - pos.z, inter_p.x - pos.x); // between -pi ~ pi
        if (u < 0)
            u += 2 * pi;
        return Vec(u, t);
    }
    double get_sphere_intersect(Ray ray, Vec o, double r) {
        Vec ro = o - ray.origin;
        double b = ray.direction.dot(ro);
        double d = square(b) - ro.dot(ro) + square(r);
        if (d < 0) return -1;
        else d = sqrt(d);
        double t = b - d > eps ? b - d : b + d > eps? b + d : -1;
        if (t < 0)
            return -1;
        return t;
    }
    virtual Hitpoint intersect(const Ray &ray) {
        double final_dis = inf;
        // check for |dy|<eps
        if (std::abs(ray.direction.y) < 5e-4)
        {
            double dis_to_axis = (Vec(pos.x, ray.origin.y, pos.z) - ray.origin).len();
            double hit = ray.get(dis_to_axis).y;
            if (hit < pos.y + eps || hit > pos.y + curve.height - eps)
                return Hitpoint::null;
            // solve function pos.y+y(t)=ray.origin.y to get x(t)
            double t = solve_t(hit - pos.y);
            if (t < 0 || t > 1)
                return Hitpoint::null;
            Vec loc = curve.getpos(t);
            double ft = pos.y + loc.y - hit;
            if (std::abs(ft) > eps)
                return Hitpoint::null;
            // assume sphere (pos.x, pos.y + loc.y, pos.z) - loc.x
            final_dis = get_sphere_intersect(ray, Vec(pos.x, pos.y + loc.y, pos.z), loc.x);
            if (final_dis < 0)
                return Hitpoint::null;
            Vec inter_p = ray.get(final_dis);
            // printf("y %f small!!!",std::abs((inter_p - Vec(pos.x, inter_p.y, pos.z)).len2() - square(loc.x)));
            if (std::abs((inter_p - Vec(pos.x, inter_p.y, pos.z)).len2() - square(loc.x)) > 1e-1)
                return Hitpoint::null;
            // second iteration, more accuracy
            hit = inter_p.y;
            if (hit < pos.y + eps || hit > pos.y + curve.height - eps)
                return Hitpoint::null;
            t = solve_t(hit - pos.y);
            loc = curve.getpos(t);
            ft = pos.y + loc.y - hit;
            if (std::abs(ft) > eps)
                return Hitpoint::null;
            final_dis = get_sphere_intersect(ray, Vec(pos.x, pos.y + loc.y, pos.z), loc.x);
            if (final_dis < 0)
                return Hitpoint::null;
            inter_p = ray.get(final_dis);
            if (std::abs((inter_p - Vec(pos.x, hit, pos.z)).len2() - square(loc.x)) > 1e-2)
                return Hitpoint::null;
            // printf("---y %f small!!!",std::abs((inter_p - Vec(pos.x, inter_p.y, pos.z)).len2() - square(loc.x)));
            Hitpoint res;
            res.object = static_cast<Object*>(this);
            res.distance = final_dis;
            res.position = inter_p;
            res.norm = norm(inter_p);
            return res;
        }

        double a = 0, b = 0, c = 0, t1, t2;
        // (xo-x'+xd/yd*(y-yo))^2 -> (t1+t2*y)^2
        t1 = ray.origin.x - pos.x - ray.direction.x / ray.direction.y * ray.origin.y;
        t2 = ray.direction.x / ray.direction.y;
        a += t2 * t2;
        b += 2 * t1 * t2;
        c += t1 * t1;
        // (zo-z'+zd/yd*(y-yo))^2 -> (t1+t2*y)^2
        t1 = ray.origin.z - pos.z - ray.direction.z / ray.direction.y * ray.origin.y;
        t2 = ray.direction.z / ray.direction.y;
        a += square(t2);
        b += 2 * t1 * t2;
        c += square(t1);
        // ay^2+by+c -> a'(y-b')^2+c'
        c = c - b * b / 4 / a;
        b = -b / 2 / a - pos.y;
        // printf("%lf %lf %lf\n",a,b,c);
        if (0 <= b && b <= curve.height && c > curve.max2
            || (b < 0 || b > curve.height) && std::min(square(b), square(curve.height - b)) * a + c > curve.max2) // no intersect
            return Hitpoint::null;

        for(int ind = 0; ind <= curve.num; ++ind)
        {
            double t0 = curve.data[ind].t0, t1 = curve.data[ind].t1;
            {
                check(t0, t1, (t0 + t1 + t0) / 3, ray, a, b, c, final_dis);
                check(t0, t1, (t1 + t0 + t1) / 3, ray, a, b, c, final_dis);
            }
        }
        if (final_dis < inf / 2){
            Hitpoint res;
            res.distance = final_dis;
            res.position = ray.get(final_dis);
            res.norm = norm(res.position);
            res.object = static_cast<Object*>(this);
            return res;
        }

        else
            return Hitpoint::null;
    }
    double newton(double t, double a, double b, double c, double low=eps, double upp=1-eps)
    {
        // solve sqrt(a(y(t)+pos.y-b)^2+c)=x(t)
        // f(t) = x(t) - sqrt(a(y(t)+pos.y-b)^2+c)
        // f'(t) = x'(t) - a(y(t)+pos.y-b)*y'(t) / sqrt(...)
        // if t is not in [0, 1] then assume f(t) is a linear function
        double ft, dft, x, y, dx, dy, sq;
        Vec loc, dir;
        for (int i = 10; i--; )
        {
            if (t < 0) t = low;
            if (t > 1) t = upp;
            loc = curve.getpos(t), dir = curve.getdir(t);
            x = loc.x, dx = dir.x;
            y = loc.y, dy = dir.y;
            // printf("%lf %lf %lf\n",t,x,y);
            sq = sqrt(a * square(y - b) + c);
            ft = x - sq;
            dft = dx - a * (y - b) * dy / sq;
            if (std::abs(ft) < eps)
                return t;
            t -= ft / dft;
        }
        return -1;
    }

    bool check(double low, double upp, double init, Ray ray, double a, double b, double c, double&final_dis)
    {
        double t = newton(init, a, b, c, low, upp);
        if (t <= 0 || t >= 1)
            return false;
        Vec loc = curve.getpos(t);
        double x = loc.x, y = loc.y;
        double ft = x - sqrt(a * square(y - b) + c);
        if (std::abs(ft) > eps)
            return false;
        // calc t for ray
        double dis = (pos.y + y - ray.origin.y) / ray.direction.y;
        if (dis < eps)
            return false;
        Vec inter_p = ray.get(dis);
        if (std::abs((Vec(pos.x, pos.y + y, pos.z) - inter_p).len2() - x * x) > eps)
            return false;
        if (dis < final_dis)
        {
            final_dis = dis;
            return true;
        }
        return false;
    }
    double getft(double t, double a, double b, double c)
    {
        if (t < 0) t = eps;
        if (t > 1) t = 1 - eps;
        Vec loc = curve.getpos(t);
        double x = loc.x, y = loc.y;
        return x - sqrt(a * square(y - b) + c);
    }
    virtual std::pair<Vec, Vec> aabb() {
        return std::make_pair(Vec(pos.x - curve.max, pos.y, pos.z - curve.max), Vec(pos.x + curve.max, pos.y + curve.height, pos.z + curve.max));
    }
    virtual Vec norm(Vec p) {
        Vec tmp = change_for_bezier(p);
        Vec dir = curve.getdir(tmp.y);
        Vec d_surface = Vec(cos(tmp.x), dir.y / dir.x, sin(tmp.x));
        Vec d_circ = Vec(-sin(tmp.x), 0, cos(tmp.x));
        return cross(d_circ,d_surface).normalize();
    }
};




#endif //TRACER_FINAL_BEZIER_H

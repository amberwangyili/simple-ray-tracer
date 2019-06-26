//
// Created by YILI WANG on 6/24/19.
//

#ifndef TRACER_FINAL_VEC_H
#define TRACER_FINAL_VEC_H


#include "Common.h"
class Matrix;
class Vec;

class Vec {        // Usage: time ./explicit 16 && xv image.ppm
public:
    typedef double value_t;
    union {
        struct {
            double x, y, z;
        };
        struct {
            double p[3];
        };
    };


    Vec(double x_=0, double y_=0, double z_=0){ x=x_; y=y_; z=z_; }

    inline double &operator[](int index) {
        return p[index];
    }

    inline const double &operator[](int index) const {
        return p[index];
    }

    inline Vec normalize() {
        double norm = len();
        x /= norm; y /= norm; z /= norm;
        return *this;
    }

    double len2() const { return x * x + y * y + z * z; }
    double len() const {return sqrt(x * x + y * y + z * z);}
    double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; } // dot


    double max() const {
        return std::max(x, std::max(y, z));
    }
    double min() const {
        return std::min(x, std::min(y, z));
    }
    inline double mean(void) const {
        return (x + y + z) / 3;
    }
    void print() const {
        printf("%.5lf %.5lf %.5lf\n", x, y, z);
    }
    static const Vec  Zero;
    static const Vec XAxis;
    static const Vec YAxis;
    static const Vec ZAxis;
    Vec rotate(Vec axis, double theta) const {
        // the following function implements a 3d rotation
        // referenced Raina's implementation
        Vec ret;
        double cost = cos( theta );
        double sint = sin( theta );
        ret.x += x * ( axis.x * axis.x + ( 1 - axis.x * axis.x ) * cost );
        ret.x += y * ( axis.x * axis.y * ( 1 - cost ) - axis.z * sint );
        ret.x += z * ( axis.x * axis.z * ( 1 - cost ) + axis.y * sint );
        ret.y += x * ( axis.y * axis.x * ( 1 - cost ) + axis.z * sint );
        ret.y += y * ( axis.y * axis.y + ( 1 - axis.y * axis.y ) * cost );
        ret.y += z * ( axis.y * axis.z * ( 1 - cost ) - axis.x * sint );
        ret.z += x * ( axis.z * axis.x * ( 1 - cost ) - axis.y * sint );
        ret.z += y * ( axis.z * axis.y * ( 1 - cost ) + axis.x * sint );
        ret.z += z * ( axis.z * axis.z + ( 1 - axis.z * axis.z ) * cost );
        return ret;
    }

};

Vec operator*(const Matrix &a, const Vec &b);
Vec operator+(const Vec &a, const Vec &b);
Vec operator-(const Vec &a, const Vec &b);
Vec operator*(const Vec &a, const Vec &b);
Vec operator/(const Vec &a, const Vec &b);
Vec operator*(const Vec &a, const double &b);
Vec operator/(const Vec &a, const double &b);
Vec operator*( const double &b,const Vec &a);
Vec operator/( const double &b,const Vec &a);


inline Vec operator-(const Vec &op0) {
    return Vec(-op0.x, -op0.y, -op0.z);
}

bool operator<=(const Vec &a, const Vec &b);
bool operator>=(const Vec &a, const Vec &b);
double dot(const Vec &a, const Vec &b);
Vec cross(const Vec &a, const Vec &b) {
    return Vec(a.y * b.z - b.y * a.z, b.x * a.z - a.x * b.z, a.x * b.y - b.x * a.y);
}
double det(const Vec &a, const Vec &b, const Vec &c);
inline Vec reflect(const Vec &in, const Vec &norm) {
    return in - norm * 2 * dot(norm, in);
}
inline std::ostream &operator<<(std::ostream &os, const Vec &vec) {
    os << "Vector(" << vec.x << " " << vec.y << " " << vec.z << ")";
    return os;
}
class Matrix {
public:
    int n;
    double **v;
    Matrix(int n) {
        this->n = n;
        v = new double*[n];
        for (int i = 0; i < n; ++i)
            v[i] = new double[n];
    }
    ~Matrix() {
        for (int i = 0; i < n; ++i)
            delete[] v[i];
        delete[] v;
    }
    Matrix inv();
    double& at(int i, int j) {
        return v[i - 1][j - 1];
    }
};

Matrix Matrix::inv() {
    Matrix T(n);
    T.v[0][0] = at(2, 2) * at(3, 3) - at(2, 3) * at(3, 2);
    T.v[0][1] = at(1, 3) * at(3, 2) - at(1, 2) * at(3, 3);
    T.v[0][2] = at(1, 2) * at(2, 3) - at(1, 3) * at(2, 2);
    T.v[1][0] = at(2, 3) * at(3, 1) - at(2, 1) * at(3, 3);
    T.v[1][1] = at(1, 1) * at(3, 3) - at(1, 3) * at(3, 1);
    T.v[1][2] = at(2, 1) * at(1, 3) - at(1, 1) * at(2, 3);
    T.v[2][0] = at(2, 1) * at(3, 2) - at(2, 2) * at(3, 1);
    T.v[2][1] = at(1, 2) * at(3, 1) - at(1, 1) * at(3, 2);
    T.v[2][2] = at(1, 1) * at(2, 2) - at(2, 1) * at(1, 2);
    double d = at(1, 1) * (at(2, 2) * at(3, 3) - at(2, 3) * at(3, 2))
               - at(2, 1) * (at(1, 2) * at(3, 3) - at(1, 3) * at(3, 2))
               + at(3, 1) * (at(1, 2) * at(2, 3) - at(1, 3) * at(2, 2));
    for (int i = 0; i <= 2; ++i)
        for (int j = 0; j <= 2; ++j)
            T.v[i][j] /= d;
    return T;
}

Vec operator*(const Matrix &a, const Vec &b) {
    return Vec(
            a.v[0][0] * b.x + a.v[0][1] * b.y + a.v[0][2] * b.z,
            a.v[1][0] * b.x + a.v[1][1] * b.y + a.v[1][2] * b.z,
            a.v[2][0] * b.x + a.v[2][1] * b.y + a.v[2][2] * b.z
    );
}

Vec operator+(const Vec &a, const Vec &b) {
    return Vec(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vec operator-(const Vec &a, const Vec &b) {
    return Vec(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vec operator*(const Vec &a, const Vec &b) {
    return Vec(a.x * b.x, a.y * b.y, a.z * b.z);
}

Vec operator/(const Vec &a, const Vec &b) {
    return Vec(a.x / b.x, a.y / b.y, a.z / b.z);
}

Vec operator*(const Vec &a, const double &b) {
    return Vec(a.x * b, a.y * b, a.z * b);
}

Vec operator/(const Vec &a, const double &b) {
    return Vec(a.x / b, a.y / b, a.z / b);
}

Vec operator*( const double &b,const Vec &a){
    return Vec(a.x * b, a.y * b, a.z * b);
}
Vec operator/( const double &b,const Vec &a){
    return Vec(a.x / b, a.y / b, a.z / b);
}



bool operator<=(const Vec &a, const Vec &b) {
    return a.x <= b.x + epsilon && a.y <= b.y + epsilon && a.z <= b.z + epsilon;
}

bool operator>=(const Vec &a, const Vec &b) {
    return a.x + epsilon >= b.x && a.y + epsilon >= b.y && a.z + epsilon >= b.z;
}

double dot(const Vec &a, const Vec &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double det(const Vec &a, const Vec &b, const Vec &c) {
    return a.x * (b.y * c.z - b.z * c.y)
           - b.x * (a.y * c.z - a.z * c.y)
           + c.x * (a.y * b.z - a.z * b.y);
}
inline Vec &operator+=(Vec &op0, const Vec &op1) {
    op0.x += op1.x, op0.y += op1.y, op0.z += op1.z;
    return op0;
}

inline Vec &operator-=(Vec &op0, const Vec &op1) {
    op0.x -= op1.x, op0.y -= op1.y, op0.z -= op1.z;
    return op0;
}

const Vec Vec::XAxis(1, 0, 0);
const Vec Vec::YAxis(0, 1, 0);
const Vec Vec::ZAxis(0, 0, 1);
const Vec Vec::Zero(0, 0, 0);

#endif //TRACER_FINAL_VEC_H

#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <stack>
class point;
class mat4;
point transformPoint(const mat4&, const point&);
mat4 getRotationMatrix(point, double);
point applyRodrigues(const point&, const point&, const double&);
mat4 getTranslationMatrix(const point&);
mat4 getViewTransformationMatrix(const point&, const point&, const point&);
mat4 getProjectionMatrix(double, double, double, double);

class point {
   public:
    double x, y, z;
    point(double x, double y, double z)
        : x{x}, y{y}, z{z} {}
    point()
        : x{0}, y{0}, z{0} {}
    point operator+(const point& p) const {
        return point(x + p.x, y + p.y, z + p.z);
    }
    point operator-(const point& p) const {
        return point(x - p.x, y - p.y, z - p.z);
    }
    point operator*(const double& d) const {
        return point(x * d, y * d, z * d);
    }
    point operator/(const double& d) const {
        return point(x / d, y / d, z / d);
    }
    point operator-() const {
        return point(-x, -y, -z);
    }
    double dot(const point& p) const {
        return x * p.x + y * p.y + z * p.z;
    }
    point cross(const point& p) const {
        return point(y * p.z - z * p.y, z * p.x - x * p.z,
                     x * p.y - y * p.x);
    }
    double norm() const {
        return sqrt(x * x + y * y + z * z);
    }
    friend std::ostream& operator<<(std::ostream& os,
                                    const point& p) {
        os << p.x << " " << p.y << " " << p.z;
        return os;
    }
};

class mat4 {
   public:
    double m[4][4];
    mat4() {
        for (int i = 0; i < 4; i++) {
            for (int k = 0; k < 4; k++)
                m[i][k] = 0;
        }
    }
    mat4 operator*(const mat4& mat) const {
        mat4 res{};
        for (int i = 0; i < 4; i++) {
            for (int k = 0; k < 4; k++) {
                for (int j = 0; j < 4; j++) {
                    res.m[i][j] += m[i][k] * mat.m[k][j];
                }
            }
        }
        return res;
    }

    void setIdentity() {
        for (int i = 0; i < 4; i++) {
            for (int k = 0; k < 4; k++) {
                m[i][k] = (i == k) ? 1 : 0;
            }
        }
    }
    void setTranslate(double x, double y, double z) {
        setIdentity();
        m[0][3] = x;
        m[1][3] = y;
        m[2][3] = z;
    }
    void setScale(double x, double y, double z) {
        setIdentity();
        m[0][0] = x;
        m[1][1] = y;
        m[2][2] = z;
    }
    void setRotate(const point& c1, const point& c2, const point& c3) {
        setIdentity();
        m[0][0] = c1.x;
        m[1][0] = c1.y;
        m[2][0] = c1.z;
        m[0][1] = c2.x;
        m[1][1] = c2.y;
        m[2][1] = c2.z;
        m[0][2] = c3.x;
        m[1][2] = c3.y;
        m[2][2] = c3.z;
    }
    friend std::ostream& operator<<(std::ostream& os, const mat4& mat) {
        for (int i = 0; i < 4; i++) {
            os << "[";
            for (int k = 0; k < 4; k++) {
                os << mat.m[i][k];
                if (k != 3) os << ", ";
            }
            os << "]\n";
        }
        return os;
    }
};

point transformPoint(const mat4& mat, const point& p) {
    double x = p.x, y = p.y, z = p.z;
    double w = mat.m[3][0] * x + mat.m[3][1] * y + mat.m[3][2] * z + mat.m[3][3];
    return point((mat.m[0][0] * x + mat.m[0][1] * y + mat.m[0][2] * z + mat.m[0][3]) / w,
                 (mat.m[1][0] * x + mat.m[1][1] * y + mat.m[1][2] * z + mat.m[1][3]) / w,
                 (mat.m[2][0] * x + mat.m[2][1] * y + mat.m[2][2] * z + mat.m[2][3]) / w);
}

point applyRodrigues(const point& x, const point& a, const double& theta) {
    return (x * cos(theta) + (a * (a.dot(x))) * (1 - cos(theta)) + a.cross(x) * sin(theta));
}

mat4 getRotationMatrix(point a, double theta) {
    theta = theta * M_PI / 180;
    a = a / a.norm();
    point c1 = applyRodrigues({1, 0, 0}, a, theta);
    point c2 = applyRodrigues({0, 1, 0}, a, theta);
    point c3 = applyRodrigues({0, 0, 1}, a, theta);
    mat4 rotationMatrix;
    rotationMatrix.setRotate(c1, c2, c3);
    return rotationMatrix;
}

mat4 getTranslationMatrix(const point& p) {
    mat4 translationMatrix;
    translationMatrix.setTranslate(p.x, p.y, p.z);
    return translationMatrix;
}

mat4 getScaleMatrix(point p) {
    mat4 scaleMatrix;
    scaleMatrix.setScale(p.x, p.y, p.z);
    return scaleMatrix;
}

mat4 getViewTransformationMatrix(const point& eye, const point& look, const point& up) {
    point l = look - eye;
    l = l / l.norm();
    point r = l.cross(up);
    r = r / r.norm();
    point u = r.cross(l);
    mat4 T = getTranslationMatrix(-eye);
    mat4 R;
    R.setIdentity();
    R.m[0][0] = r.x;
    R.m[0][1] = r.y;
    R.m[0][2] = r.z;
    R.m[1][0] = u.x;
    R.m[1][1] = u.y;
    R.m[1][2] = u.z;
    R.m[2][0] = -l.x;
    R.m[2][1] = -l.y;
    R.m[2][2] = -l.z;
    return R * T;
}

mat4 getProjectionMatrix(double fovY, double aspectRatio, double near, double far) {
    mat4 P;
    double fovX = fovY * aspectRatio;
    double t = near * tan(fovY * M_PI / 360);
    double r = near * tan(fovX * M_PI / 360);
    P.m[0][0] = near / r;
    P.m[1][1] = near / t;
    P.m[2][2] = -(far + near) / (far - near);
    P.m[2][3] = -(2 * far * near) / (far - near);
    P.m[3][2] = -1;
    return P;
}
#pragma once
#include "utility.h"
#define EPSILON 0.000001

class point {
   public:
    GLdouble x, y, z;
    point(double x, double y, double z)
        : x{x}, y{y}, z{z} {}
    point()
        : x{0}, y{0}, z{0} {}
    point operator-(const point& p) const {
        point temp;
        temp.x = x - p.x;
        temp.y = y - p.y;
        temp.z = z - p.z;
        return temp;
    }
    point operator+(const point& p) const {
        point temp;
        temp.x = x + p.x;
        temp.y = y + p.y;
        temp.z = z + p.z;
        return temp;
    }
    point operator=(const point& p) {
        x = p.x;
        y = p.y;
        z = p.z;
        return *this;
    }
    point operator*(const double& d) const {
        point temp;
        temp.x = x * d;
        temp.y = y * d;
        temp.z = z * d;
        return temp;
    }
    point normalize() {
        GLdouble scale = sqrt(x * x + y * y + z * z);
        if (scale >= 0.01) {
            x /= scale;
            y /= scale;
            z /= scale;
        }
        return *this;
    }
    point cross(const point& p) const {
        return point(y * p.z - z * p.y, z * p.x - x * p.z,
                     x * p.y - y * p.x);
    }
    double norm() const {
        return sqrt(x * x + y * y + z * z);
    }
    double dot(point b) const {
        return (x * b.x + y * b.y + z * b.z);
    }
};

class Intersection {
   public:
    bool doesIntersect;
    double t;
    point ip;
};

class Ray {
   public:
    point origin;
    point direction;
    Ray(point origin, point direction)
        : origin{origin}, direction{direction} {}
    Ray()
        : origin{}, direction{} {}
};

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
    point operator-() const {
        point temp;
        temp.x = -x;
        temp.y = -y;
        temp.z = -z;
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
    double angleBetween(point b) const {
        return acos(dot(b) / (norm() * b.norm()));
    }
};
class Color {
   public:
    double r, g, b;
    Color(double r, double g, double b)
        : r{r}, g{g}, b{b} {}
    Color()
        : r{0}, g{0}, b{0} {}
    Color operator*(const double& d) const {
        Color temp;
        temp.r = r * d;
        temp.g = g * d;
        temp.b = b * d;
        return temp;
    }
    Color operator+(const Color& c) const {
        Color temp;
        temp.r = r + c.r;
        temp.g = g + c.g;
        temp.b = b + c.b;
        return temp;
    }
    Color operator+=(const Color& c) {
        r += c.r;
        g += c.g;
        b += c.b;
        return *this;
    }
};

class Coeffs {
   public:
    double ka, kd, ks, kr, shine;
    Coeffs(double ka, double kd, double ks, double kr, double shine)
        : ka{ka}, kd{kd}, ks{ks}, kr{kr}, shine{shine} {}
    Coeffs()
        : ka{0}, kd{0}, ks{0}, kr{0}, shine{0} {}
};
class Intersection {
   public:
    bool doesIntersect;
    double t;
    point ip;
    Color color;
    Coeffs coeffs;
    point normal;
    Intersection(bool doesIntersect, double t, point ip, Color color, Coeffs coeffs, point normal)
        : doesIntersect{doesIntersect}, t{t}, ip{ip}, color{color}, coeffs{coeffs}, normal{normal} {}
    Intersection()
        : doesIntersect{false}, t{0}, ip{}, color{}, coeffs{} {}
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

bool rayTriangleIntersect(const Ray& ray, const point& v0, const point& v1, const point& v2, double& t) {
    point e1 = v1 - v0;
    point e2 = v2 - v0;
    point h = ray.direction.cross(e2);
    double a = e1.dot(h);

    if (a > -1e-6 && a < 1e-6)
        return false;

    double f = 1.0 / a;
    point s = ray.origin - v0;
    double u = f * s.dot(h);

    if (u < 0.0 || u > 1.0)
        return false;

    point q = s.cross(e1);
    double v = f * ray.direction.dot(q);

    if (v < 0.0 || u + v > 1.0)
        return false;

    t = f * e2.dot(q);
    if (t > 1e-6)
        return true;

    return false;
}

bool rayQuadIntersect(const Ray& ray, const point& v0, const point& v1, const point& v2, const point& v3, double& t) {
    if (rayTriangleIntersect(ray, v0, v1, v2, t))
        return true;
    if (rayTriangleIntersect(ray, v0, v2, v3, t))
        return true;
    return false;
}
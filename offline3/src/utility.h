#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <GL/glut.h>
#include "bitmap_image.hpp"
#include "helper.h"

using namespace std;
// Class signatures
class CheckerBoard;
class point;
class Ray;
class Intersection;
class Sphere;
class Pyramid;
class Cube;
class NormalLightSource;
class SpotLightSource;
class Color;
class Coeffs;
// Function signatures
void readInputFile(string fileName);
void drawAxes();
void drawCheckerBoard();
void generateRays();
void capture();
Color traceRay(const Ray& ray, int levelOfRecursion);
Intersection closestIntersectionWithObjects(const Ray& ray);
// Helper variables
float PI = acos(-1);
const int LEN = 2000;  // axis length
// Global variables
extern point pos;  // position of the eye
extern point l;    // look/forward direction
extern point r;    // right direction
extern point u;    // up direction
extern bool isAxes;

double nearPlane, farPlane, fovY, aspectRatio;
int levelOfRecursion, width, height;
int numObjects, numNormLights, numSpotLights;
vector<Sphere> spheres;
vector<Pyramid> pyramids;
vector<Cube> cubes;
vector<NormalLightSource> normLights;
vector<SpotLightSource> spotLights;
vector<vector<Ray>> rays;

// class definitions
class CheckerBoard {
   public:
    int width;
    double ka, kd, kr;
    // override cout
    friend ostream& operator<<(ostream& os, const CheckerBoard& cboard) {
        os << "CheckerBoard: " << cboard.width << endl;
        os << "ka: " << cboard.ka << " kd: " << cboard.kd << " kr: " << cboard.kr << endl;
        return os;
    }

    Intersection intersect(Ray ray) const {
        Intersection intersection;
        point origin = ray.origin;
        point direction = ray.direction;
        point normal = point(0, 0, 1);
        double t = -origin.dot(normal) / direction.dot(normal);
        point ip = origin + direction * t;
        double dist = (ip - origin).norm();
        double theta = l.angleBetween(ip - origin);
        if (t < 0 || dist > (farPlane / cos(theta))) {
            intersection.doesIntersect = false;
            return intersection;
        }
        if (abs(ip.x) > LEN || abs(ip.y) > LEN) {
            intersection.doesIntersect = false;
            return intersection;
        }
        intersection.doesIntersect = true;
        intersection.t = t;
        intersection.ip = ip;
        int xSquare = (int)floor(intersection.ip.x / width);
        int ySquare = (int)floor(intersection.ip.y / width);
        bool isWhite = (xSquare + ySquare) % 2 == 0;
        intersection.color = isWhite ? Color(1, 1, 1) : Color(0, 0, 0);
        intersection.coeffs = Coeffs(this->ka, this->kd, 0, this->kr, 0);
        intersection.normal = normal;
        return intersection;
    }
};
CheckerBoard cboard;
class Sphere {
   public:
    double x, y, z, radius;
    double r, g, b;
    double ka, kd, ks, kr;
    double shine;

    Sphere()
        : x{0}, y{0}, z{0}, radius{10}, r{255}, g{255}, b{255} {}
    // override cout
    friend ostream& operator<<(ostream& os, const Sphere& sphere) {
        os << "Sphere: " << sphere.x << " " << sphere.y << " " << sphere.z << " " << sphere.radius << endl;
        os << "Color: " << sphere.r << " " << sphere.g << " " << sphere.b << endl;
        os << "ka: " << sphere.ka << " kd: " << sphere.kd << " ks: " << sphere.ks << " kr: " << sphere.kr << endl;
        os << "shine: " << sphere.shine << endl;
        return os;
    }

    void draw() const {
        glColor3f(r, g, b);
        glPushMatrix();
        glTranslated(x, y, z);
        glutSolidSphere(radius, 64, 64);
        glPopMatrix();
    }
    Intersection intersect(Ray ray) const {
        Intersection intersection;
        point origin = ray.origin;
        point direction = ray.direction;
        point center = point(x, y, z);
        point oc = origin - center;
        double a = direction.dot(direction);
        double b = 2 * oc.dot(direction);
        double c = oc.dot(oc) - radius * radius;
        double discriminant = b * b - 4 * a * c;
        if (discriminant < 0) {
            intersection.doesIntersect = false;
            return intersection;
        }
        double t1 = (-b + sqrt(discriminant)) / (2 * a);
        double t2 = (-b - sqrt(discriminant)) / (2 * a);
        if (t1 < 0 && t2 < 0) {
            intersection.doesIntersect = false;
            return intersection;
        }
        if (t1 < 0) {
            intersection.t = t2;
        } else if (t2 < 0) {
            intersection.t = t1;
        } else {
            intersection.t = min(t1, t2);
        }
        intersection.ip = origin + direction * intersection.t;
        intersection.doesIntersect = true;
        intersection.color = Color(this->r, this->g, this->b);
        intersection.coeffs = Coeffs(this->ka, this->kd, this->ks, this->kr, this->shine);
        intersection.normal = (intersection.ip - center).normalize();
        return intersection;
    }
};

class Pyramid {
   public:
    double x, y, z, width, height;
    double r, g, b;
    double ka, kd, ks, kr;
    double shine;
    // override cout
    friend ostream& operator<<(ostream& os, const Pyramid& pyramid) {
        os << "Pyramid: " << pyramid.x << " " << pyramid.y << " " << pyramid.z << " " << pyramid.width << " " << pyramid.height << endl;
        os << "Color: " << pyramid.r << " " << pyramid.g << " " << pyramid.b << endl;
        os << "ka: " << pyramid.ka << " kd: " << pyramid.kd << " ks: " << pyramid.ks << " kr: " << pyramid.kr << endl;
        os << "shine: " << pyramid.shine << endl;
        return os;
    }

    void draw() const {
        glColor3f(r, g, b);
        glBegin(GL_QUADS);
        // Bottom face
        glVertex3f(x, y, z);
        glVertex3f(x + width, y, z);
        glVertex3f(x + width, y + width, z);
        glVertex3f(x, y + width, z);
        glEnd();

        glBegin(GL_TRIANGLES);
        // Front
        glVertex3f(x, y, z);
        glVertex3f(x + width, y, z);
        glVertex3f(x + width / 2, y + width / 2, z + height);
        // Right
        glVertex3f(x + width, y, z);
        glVertex3f(x + width / 2, y + width / 2, z + height);
        glVertex3f(x + width, y + width, z);
        // Left
        glVertex3f(x, y, z);
        glVertex3f(x + width / 2, y + width / 2, z + height);
        glVertex3f(x, y + width, z);
        // Back
        glVertex3f(x, y + width, z);
        glVertex3f(x + width, y + width, z);
        glVertex3f(x + width / 2, y + width / 2, z + height);
        glEnd();
    }

    Intersection intersect(const Ray& ray) const {
        Intersection intersection;
        intersection.doesIntersect = false;
        double tMin = std::numeric_limits<double>::max();

        point apex = point(x + width / 2, y + width / 2, z + height);
        point v0 = apex;
        point v1 = point(x, y, z);
        point v2 = point(x + width, y, z);
        point v3 = point(x + width, y + width, z);
        point v4 = point(x, y + width, z);
        point normal;

        // Check intersection with each triangle
        double t;
        if (rayTriangleIntersect(ray, v0, v1, v2, t) && t < tMin) {
            tMin = t;
            intersection.doesIntersect = true;
            normal = (v1 - v0).cross(v2 - v0).normalize();
        }
        if (rayTriangleIntersect(ray, v0, v2, v3, t) && t < tMin) {
            tMin = t;
            intersection.doesIntersect = true;
            normal = (v2 - v0).cross(v3 - v0).normalize();
        }
        if (rayTriangleIntersect(ray, v0, v3, v4, t) && t < tMin) {
            tMin = t;
            intersection.doesIntersect = true;
            normal = (v3 - v0).cross(v4 - v0).normalize();
        }
        if (rayTriangleIntersect(ray, v0, v4, v1, t) && t < tMin) {
            tMin = t;
            intersection.doesIntersect = true;
            normal = (v4 - v0).cross(v1 - v0).normalize();
        }

        // Check intersection with base quad (2 triangles)
        if (rayTriangleIntersect(ray, v1, v2, v3, t) && t < tMin) {
            tMin = t;
            intersection.doesIntersect = true;
            normal = (v2 - v1).cross(v3 - v1).normalize();
        }
        if (rayTriangleIntersect(ray, v1, v3, v4, t) && t < tMin) {
            tMin = t;
            intersection.doesIntersect = true;
            normal = (v3 - v1).cross(v4 - v1).normalize();
        }

        if (intersection.doesIntersect) {
            intersection.t = tMin;
            intersection.ip = ray.origin + ray.direction * tMin;
            intersection.color = Color(this->r, this->g, this->b);
            intersection.coeffs = Coeffs(this->ka, this->kd, this->ks, this->kr, this->shine);
            intersection.normal = normal;
        }

        return intersection;
    }
};

class Cube {
   public:
    double x, y, z;
    double side;
    double r, g, b;
    double ka, kd, ks, kr;
    double shine;
    // override cout
    friend ostream& operator<<(ostream& os, const Cube& cube) {
        os << "Cube: " << cube.x << " " << cube.y << " " << cube.z << " " << cube.side << endl;
        os << "Color: " << cube.r << " " << cube.g << " " << cube.b << endl;
        os << "ka: " << cube.ka << " kd: " << cube.kd << " ks: " << cube.ks << " kr: " << cube.kr << endl;
        os << "shine: " << cube.shine << endl;
        return os;
    }

    void draw() const {
        glColor3f(r, g, b);
        glBegin(GL_QUADS);  // Begin drawing the color cube with 6 quads
        // Draw the cube where x,y,z are the bottom lower left coordinate, side is the length of each arm
        // Bottom face
        glVertex3f(x, y, z);
        glVertex3f(x + side, y, z);
        glVertex3f(x + side, y + side, z);
        glVertex3f(x, y + side, z);

        // Top face
        glVertex3f(x, y, z + side);
        glVertex3f(x + side, y, z + side);
        glVertex3f(x + side, y + side, z + side);
        glVertex3f(x, y + side, z + side);

        // Right face
        glVertex3f(x + side, y, z);
        glVertex3f(x + side, y + side, z);
        glVertex3f(x + side, y + side, z + side);
        glVertex3f(x + side, y, z + side);

        // Left face
        glVertex3f(x, y, z);
        glVertex3f(x, y + side, z);
        glVertex3f(x, y + side, z + side);
        glVertex3f(x, y, z + side);

        // Front face
        glVertex3f(x, y, z);
        glVertex3f(x + side, y, z);
        glVertex3f(x + side, y, z + side);
        glVertex3f(x, y, z + side);

        // Rear face
        glVertex3f(x, y + side, z);
        glVertex3f(x + side, y + side, z);
        glVertex3f(x + side, y + side, z + side);
        glVertex3f(x, y + side, z + side);

        glEnd();  // End of drawing color-cube
    }
    Intersection intersect(Ray ray) const {
        point v0 = point(x, y, z);
        point v1 = point(x + side, y, z);
        point v2 = point(x + side, y + side, z);
        point v3 = point(x, y + side, z);
        point v4 = point(x, y, z + side);
        point v5 = point(x + side, y, z + side);
        point v6 = point(x + side, y + side, z + side);
        point v7 = point(x, y + side, z + side);
        Intersection intersection;
        // check if ray intersects with bottom face
        double t;
        if (rayQuadIntersect(ray, v0, v1, v2, v3, t)) {
            intersection.doesIntersect = true;
            intersection.t = t;
            intersection.ip = ray.origin + ray.direction * t;
            intersection.color = Color(this->r, this->g, this->b);
            intersection.coeffs = Coeffs(this->ka, this->kd, this->ks, this->kr, this->shine);
            intersection.normal = point(0, 0, -1);
            return intersection;
        }
        // check if ray intersects with top face
        if (rayQuadIntersect(ray, v4, v5, v6, v7, t)) {
            intersection.doesIntersect = true;
            intersection.t = t;
            intersection.ip = ray.origin + ray.direction * t;
            intersection.color = Color(this->r, this->g, this->b);
            intersection.coeffs = Coeffs(this->ka, this->kd, this->ks, this->kr, this->shine);
            intersection.normal = point(0, 0, 1);
            return intersection;
        }
        // check if ray intersects with right face
        if (rayQuadIntersect(ray, v1, v5, v6, v2, t)) {
            intersection.doesIntersect = true;
            intersection.t = t;
            intersection.ip = ray.origin + ray.direction * t;
            intersection.color = Color(this->r, this->g, this->b);
            intersection.coeffs = Coeffs(this->ka, this->kd, this->ks, this->kr, this->shine);
            intersection.normal = point(1, 0, 0);
            return intersection;
        }
        // check if ray intersects with left face
        if (rayQuadIntersect(ray, v0, v4, v7, v3, t)) {
            intersection.doesIntersect = true;
            intersection.t = t;
            intersection.ip = ray.origin + ray.direction * t;
            intersection.color = Color(this->r, this->g, this->b);
            intersection.coeffs = Coeffs(this->ka, this->kd, this->ks, this->kr, this->shine);
            intersection.normal = point(-1, 0, 0);
            return intersection;
        }
        // check if ray intersects with front face
        if (rayQuadIntersect(ray, v0, v1, v5, v4, t)) {
            intersection.doesIntersect = true;
            intersection.t = t;
            intersection.ip = ray.origin + ray.direction * t;
            intersection.color = Color(this->r, this->g, this->b);
            intersection.coeffs = Coeffs(this->ka, this->kd, this->ks, this->kr, this->shine);
            intersection.normal = point(0, -1, 0);
            return intersection;
        }
        // check if ray intersects with rear face
        if (rayQuadIntersect(ray, v3, v2, v6, v7, t)) {
            intersection.doesIntersect = true;
            intersection.t = t;
            intersection.ip = ray.origin + ray.direction * t;
            intersection.color = Color(this->r, this->g, this->b);
            intersection.coeffs = Coeffs(this->ka, this->kd, this->ks, this->kr, this->shine);
            intersection.normal = point(0, 1, 0);
            return intersection;
        }
        intersection.doesIntersect = false;
        return intersection;
    }
};

class NormalLightSource {
   public:
    double x, y, z;
    double decay;
    // override cout
    friend ostream& operator<<(ostream& os, const NormalLightSource& light) {
        os << "NormalLightSource: " << light.x << " " << light.y << " " << light.z << " " << light.decay << endl;
        return os;
    }
    void draw() const {
        glPushMatrix();

        glTranslated(x, y, z);

        int numSegments = 20;                              // Number of cylindrical segments
        double radius = 10.0;                              // Radius of the sphere
        double segmentWidth = 2.0 * radius / numSegments;  // Width of each segment

        for (int i = 0; i < numSegments; ++i) {
            // Calculate the z-coordinate and radii for the front and back of this segment
            double zFront = radius - i * segmentWidth;
            double zBack = zFront - segmentWidth;
            double rFront = sqrt(radius * radius - zFront * zFront);
            double rBack = sqrt(radius * radius - zBack * zBack);

            // Calculate the gray value for this segment
            double grayValue = 0.5 * (std::abs(2.0 * i / static_cast<double>(numSegments) - 1.0));

            glColor3f(grayValue, grayValue, grayValue);

            glBegin(GL_QUAD_STRIP);
            for (int j = 0; j <= 360; j += (360 / numSegments)) {
                double angle = j * M_PI / 180.0;
                glVertex3f(rFront * cos(angle), rFront * sin(angle), zFront);
                glVertex3f(rBack * cos(angle), rBack * sin(angle), zBack);
            }
            glEnd();
        }

        glPopMatrix();
    }

    bool illuminates(const Ray& objToLight) const {
        Intersection intersection = closestIntersectionWithObjects(objToLight);
        return !intersection.doesIntersect;
    }
};

class SpotLightSource {
   public:
    double x, y, z;
    double decay;
    double dx, dy, dz;
    double cutOffAngle;
    // override cout
    friend ostream& operator<<(ostream& os, const SpotLightSource& light) {
        os << "SpotLightSource: " << light.x << " " << light.y << " " << light.z << " " << light.decay << endl;
        os << "Direction: " << light.dx << " " << light.dy << " " << light.dz << endl;
        os << "cutOffAngle: " << light.cutOffAngle << endl;
        return os;
    }
    void draw() const {
        // Save the current transformation matrix
        glPushMatrix();

        // Translate to the position of the light source
        glTranslated(x, y, z);

        double angle = acos(dz / sqrt(dx * dx + dy * dy + dz * dz)) * 180.0 / M_PI;
        double axisLength = sqrt(dx * dx + dy * dy);
        glRotated(angle, -dy / axisLength, dx / axisLength, 0);

        // Draw the cone
        double height = 25.0;
        // double baseRadius = tan(cutOffAngle * M_PI / 180.0) * height;
        double baseRadius = tan(cutOffAngle * M_PI / 180.0) * 10;
        int slices = 32;

        // Draw curved surface using individual triangles
        glBegin(GL_TRIANGLES);
        for (int i = 0; i < slices; ++i) {
            double theta1 = i * 2.0 * M_PI / slices;
            double theta2 = (i + 1) * 2.0 * M_PI / slices;

            double x1 = baseRadius * cos(theta1);
            double y1 = baseRadius * sin(theta1);
            double x2 = baseRadius * cos(theta2);
            double y2 = baseRadius * sin(theta2);

            // Apex
            glColor3f(1.0f, 1.0f, 1.0f);
            glVertex3f(0.0f, 0.0f, 0.0f);
            glColor3f(0.5f, 0.5f, 0.5f);
            glVertex3f(x1, y1, height);
            glVertex3f(x2, y2, height);
        }
        glEnd();

        // Restore the transformation matrix to its previous state
        glPopMatrix();
    }

    bool illuminates(const Ray& objToLight) const {
        Intersection intersection = closestIntersectionWithObjects(objToLight);
        if (intersection.doesIntersect) {
            return false;
        }
        point lightDir = point(x, y, z) - objToLight.origin;
        lightDir.normalize();
        double angle = (lightDir.angleBetween(point(dx, dy, dz))) * (180 / M_PI);
        if (angle > cutOffAngle) {
            return false;
        }
        return true;
    }
};

void readInputFile(string fileName) {
    ifstream fin = ifstream(fileName);
    string line;
    fin >> nearPlane >> farPlane >> fovY >> aspectRatio;
    fin >> levelOfRecursion;
    fin >> width;
    height = width / aspectRatio;
    fin >> cboard.width;
    fin >> cboard.ka >> cboard.kd >> cboard.kr;

    fin >> numObjects;

    for (int i = 0; i < numObjects; i++) {
        string type;
        fin >> type;
        if (type == "sphere") {
            Sphere sphere;
            fin >> sphere.x >> sphere.y >> sphere.z >> sphere.radius;
            fin >> sphere.r >> sphere.g >> sphere.b;
            fin >> sphere.ka >> sphere.kd >> sphere.ks >> sphere.kr;
            fin >> sphere.shine;
            spheres.push_back(sphere);
        } else if (type == "pyramid") {
            Pyramid pyramid;
            fin >> pyramid.x >> pyramid.y >> pyramid.z >> pyramid.width >> pyramid.height;
            fin >> pyramid.r >> pyramid.g >> pyramid.b;
            fin >> pyramid.ka >> pyramid.kd >> pyramid.ks >> pyramid.kr;
            fin >> pyramid.shine;
            pyramids.push_back(pyramid);
        } else if (type == "cube") {
            Cube cube;
            fin >> cube.x >> cube.y >> cube.z >> cube.side;
            fin >> cube.r >> cube.g >> cube.b;
            fin >> cube.ka >> cube.kd >> cube.ks >> cube.kr;
            fin >> cube.shine;
            cubes.push_back(cube);
        }
    }
    fin >> numNormLights;
    for (int i = 0; i < numNormLights; i++) {
        NormalLightSource light;
        fin >> light.x >> light.y >> light.z >> light.decay;
        normLights.push_back(light);
    }
    fin >> numSpotLights;
    for (int i = 0; i < numSpotLights; i++) {
        SpotLightSource light;
        fin >> light.x >> light.y >> light.z >> light.decay;
        fin >> light.dx >> light.dy >> light.dz;
        fin >> light.cutOffAngle;
        spotLights.push_back(light);
    }
    fin.close();

    cout << "nearPlane: " << nearPlane << " farPlane: " << farPlane << " fovY: " << fovY << " aspectRatio: " << aspectRatio << endl;
    cout << "levelOfRecursion: " << levelOfRecursion << endl;
    cout << "width: " << width << " height: " << height << endl;
    cout << cboard;
    cout << "numObjects: " << numObjects << endl;
    cout << "====================================" << endl;
    for (const auto& sphere : spheres) {
        cout << sphere;
    }
    cout << "====================================" << endl;
    for (const auto& pyramid : pyramids) {
        cout << pyramid;
    }
    cout << "====================================" << endl;
    for (const auto& cube : cubes) {
        cout << cube;
    }
    cout << "====================================" << endl;
    cout << "numNormLights: " << numNormLights << endl;
    for (int i = 0; i < numNormLights; i++) {
        cout << normLights[i];
    }
    cout << "====================================" << endl;
    cout << "numSpotLights: " << numSpotLights << endl;
    for (int i = 0; i < numSpotLights; i++) {
        cout << spotLights[i];
    }
}

void drawAxes() {
    glLineWidth(3);
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);  // Red
    // X axis
    glVertex3f(-LEN, 0, 0);
    glVertex3f(LEN, 0, 0);

    glColor3f(0, 1, 0);  // Green
    // Y axis
    glVertex3f(0, -LEN, 0);
    glVertex3f(0, LEN, 0);

    glColor3f(0, 0, 1);  // Blue
    // Z axis
    glVertex3f(0, 0, -LEN);
    glVertex3f(0, 0, LEN);
    glEnd();
}

void drawCheckerBoard() {
    bool black = 0;
    for (int i = -LEN; i <= LEN; i += cboard.width) {
        glBegin(GL_QUADS);
        for (int j = -LEN; j <= LEN; j += cboard.width) {
            if (black)
                glColor3f(0, 0, 0);
            else
                glColor3f(1, 1, 1);
            black = !black;
            glVertex3f(i, j, 0);
            glVertex3f(i + cboard.width, j, 0);
            glVertex3f(i + cboard.width, j + cboard.width, 0);
            glVertex3f(i, j + cboard.width, 0);
        }
        glEnd();
    }
}

void generateRays() {
    double screenWidth = 2 * nearPlane * tan(fovY * M_PI / 360);
    double fovX = fovY * aspectRatio;
    double screenHeight = 2 * nearPlane * tan(fovX * M_PI / 360);

    point midpoint = pos + l * nearPlane;
    // Calculate the top-left corner of the screen in the 3D world
    point topLeft = midpoint + (u * (screenHeight / 2)) - (r * (screenWidth / 2));
    double du = screenWidth / width;
    double dv = screenHeight / height;
    topLeft = topLeft + r * (0.5 * du) - u * (0.5 * dv);

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            // Calculate the 3D point for this pixel in the OpenGL world
            point pt = topLeft + (r * i * du) - (u * j * dv);
            // Generate the ray for this pixel
            rays[i][j].origin = pos;
            rays[i][j].direction = (pt - pos).normalize();
        }
    }
}

void capture() {
    generateRays();
    bitmap_image img(width, height);
    img.set_all_channels(0, 0, 0);
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            Ray ray = rays[i][j];
            Color c = traceRay(ray, levelOfRecursion);
            img.set_pixel(i, j, c.r * 255, c.g * 255, c.b * 255);
        }
    }
    img.save_image("out.bmp");
    img.clear();
}

Color traceRay(const Ray& ray, int depth) {
    if (depth == 0) {
        return Color(0, 0, 0);  // Return black if depth is zero
    }
    Intersection closestIntersection = closestIntersectionWithObjects(ray);
    // closestIntersection.color = closestIntersection.color * closestIntersection.coeffs.ka;
    // double phong = 0, lambert = 0;
    // for (const auto& normLight : normLights) {
    //     point lightSource = point(normLight.x, normLight.y, normLight.z);
    //     point ip = closestIntersection.ip;
    //     Ray objToLight = Ray(ip, (lightSource - ip).normalize());
    //     if (!normLight.illuminates(objToLight)) continue;

    //     point toSource = objToLight.direction;
    //     point normal = closestIntersection.ip - point(0, 0, 0);
    //     double distance = (lightSource - ip).norm();
    //     double scaleFactor = exp(-normLight.decay * distance * distance);
    // }
    // for (const auto& spotLight : spotLights) {
    // }
    return closestIntersection.color;
}

Intersection closestIntersectionWithObjects(const Ray& ray) {
    double tMin = std::numeric_limits<double>::max();
    Intersection closestIntersection;
    // Find the closest intersection
    for (const auto& sphere : spheres) {
        Intersection is = sphere.intersect(ray);
        if (!is.doesIntersect) continue;
        if (is.t < tMin) {
            tMin = is.t;
            closestIntersection = is;
        }
    }
    for (const auto& pyramid : pyramids) {
        Intersection is = pyramid.intersect(ray);
        if (!is.doesIntersect) continue;
        if (is.t < tMin) {
            tMin = is.t;
            closestIntersection = is;
        }
    }
    for (const auto& cube : cubes) {
        Intersection is = cube.intersect(ray);
        if (!is.doesIntersect) continue;
        if (is.t < tMin) {
            tMin = is.t;
            closestIntersection = is;
        }
    }
    Intersection is = cboard.intersect(ray);
    if (is.doesIntersect && is.t < tMin) {
        tMin = is.t;
        closestIntersection = is;
    }
    return closestIntersection;
}
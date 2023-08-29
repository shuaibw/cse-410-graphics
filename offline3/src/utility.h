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
// Function signatures
void readInputFile(string fileName);
void drawAxes();
void drawCheckerBoard();
// Helper variables
float PI = acos(-1);
const int LEN = 2000;  // axis length
float angle = 0.0f;    // rotation control
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
vector<vector<point>> pointBuffer;
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
        double t = (point(0, 0, 0) - origin).dot(normal) / direction.dot(normal);
        point ip = origin + direction * t;
        if (t < 0) {
            intersection.doesIntersect = false;
            return intersection;
        }
        // if (abs(ip.x) > LEN || abs(ip.y) > LEN) {
        //     intersection.doesIntersect = false;
        //     return intersection;
        // }
        intersection.doesIntersect = true;
        intersection.t = t;
        intersection.ip = ip;
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

    // Needed for drawing sphere
    vector<float> verticesXPos;
    int subdivision;
    Sphere()
        : subdivision{4}, x{0}, y{0}, z{0}, radius{10}, r{255}, g{255}, b{255} {}
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
        drawPositiveX();  // +X

        glPushMatrix();
        glRotated(90, 0, 1, 0);
        drawPositiveX();  // -Z
        glPopMatrix();

        glPushMatrix();
        glRotated(180, 0, 1, 0);
        drawPositiveX();  // -X
        glPopMatrix();

        glPushMatrix();
        glRotated(-90, 0, 1, 0);
        drawPositiveX();  // +Z
        glPopMatrix();

        glPushMatrix();
        glRotated(90, 0, 0, 1);
        drawPositiveX();  // +Y
        glPopMatrix();

        glPushMatrix();
        glRotated(-90, 0, 0, 1);
        drawPositiveX();  // -Y
        glPopMatrix();
        glPopMatrix();
    }
    void drawPositiveX() const {
        // draw the filled sphere from vertices
        if (verticesXPos.size() == 0) {
            throw("verticesXPos is empty");
        }
        vector<float> vertices = verticesXPos;
        int pointsPerRow = (int)pow(2, subdivision) + 1;
        int totalPoints = vertices.size() / 3;
        int numRows = totalPoints / pointsPerRow;
        glBegin(GL_QUADS);
        for (int i = 0; i < totalPoints - 1 - pointsPerRow; i++) {
            float x1 = vertices[i * 3];
            float y1 = vertices[i * 3 + 1];
            float z1 = vertices[i * 3 + 2];
            float x2 = vertices[(i + 1) * 3];
            float y2 = vertices[(i + 1) * 3 + 1];
            float z2 = vertices[(i + 1) * 3 + 2];
            float x3 = vertices[(i + pointsPerRow) * 3];
            float y3 = vertices[(i + pointsPerRow) * 3 + 1];
            float z3 = vertices[(i + pointsPerRow) * 3 + 2];
            float x4 = vertices[(i + pointsPerRow + 1) * 3];
            float y4 = vertices[(i + pointsPerRow + 1) * 3 + 1];
            float z4 = vertices[(i + pointsPerRow + 1) * 3 + 2];
            glVertex3f(x1, y1, z1);
            glVertex3f(x2, y2, z2);
            glVertex3f(x4, y4, z4);
            glVertex3f(x3, y3, z3);
        }
        glEnd();
    }
    void buildUnitPositiveX() {
        const float DEG2RAD = acos(-1) / 180.0f;

        vector<float> vertices;
        float n1[3];  // normal of longitudinal plane rotating along Y-axis
        float n2[3];  // normal of latitudinal plane rotating along Z-axis
        float v[3];   // direction vector intersecting 2 planes, n1 x n2
        float a1;     // longitudinal angle along Y-axis
        float a2;     // latitudinal angle along Z-axis

        // compute the number of vertices per row, 2^n + 1
        int pointsPerRow = (int)pow(2, subdivision) + 1;

        // rotate latitudinal plane from 45 to -45 degrees along Z-axis (top-to-bottom)
        for (unsigned int i = 0; i < pointsPerRow; ++i) {
            // normal for latitudinal plane
            // if latitude angle is 0, then normal vector of latitude plane is n2=(0,1,0)
            // therefore, it is rotating (0,1,0) vector by latitude angle a2
            a2 = DEG2RAD * (45.0f - 90.0f * i / (pointsPerRow - 1));
            n2[0] = -sin(a2);
            n2[1] = cos(a2);
            n2[2] = 0;

            // rotate longitudinal plane from -45 to 45 along Y-axis (left-to-right)
            for (unsigned int j = 0; j < pointsPerRow; ++j) {
                // normal for longitudinal plane
                // if longitude angle is 0, then normal vector of longitude is n1=(0,0,-1)
                // therefore, it is rotating (0,0,-1) vector by longitude angle a1
                a1 = DEG2RAD * (-45.0f + 90.0f * j / (pointsPerRow - 1));
                n1[0] = -sin(a1);
                n1[1] = 0;
                n1[2] = -cos(a1);

                // find direction vector of intersected line, n1 x n2
                v[0] = n1[1] * n2[2] - n1[2] * n2[1];
                v[1] = n1[2] * n2[0] - n1[0] * n2[2];
                v[2] = n1[0] * n2[1] - n1[1] * n2[0];

                // normalize direction vector
                float scale = 1 / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
                v[0] *= (scale * radius);
                v[1] *= (scale * radius);
                v[2] *= (scale * radius);

                // add a vertex into array
                vertices.push_back(v[0]);
                vertices.push_back(v[1]);
                vertices.push_back(v[2]);
            }
        }

        verticesXPos = vertices;
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

    Intersection intersect(Ray ray) const {
        Intersection intersection;
        point origin = ray.origin;
        point direction = ray.direction;
        point p1 = point(x, y, z);
        point p2 = point(x + width, y, z);
        point p3 = point(x + width, y + width, z);
        point p4 = point(x, y + width, z);
        point p5 = point(x + width / 2, y + width / 2, z + height);
        point normal = (p2 - p1).cross(p3 - p1);
        double t = (p1 - origin).dot(normal) / direction.dot(normal);
        point ip = origin + direction * t;
        point v1 = p2 - p1;
        point v2 = p3 - p2;
        point v3 = p4 - p3;
        point v4 = p1 - p4;
        point v5 = p5 - p1;
        point v6 = p5 - p2;
        point v7 = p5 - p3;
        point v8 = p5 - p4;
        point n1 = v1.cross(v5);
        point n2 = v2.cross(v6);
        point n3 = v3.cross(v7);
        point n4 = v4.cross(v8);
        if (t < 0) {
            intersection.doesIntersect = false;
            return intersection;
        }
        if (n1.dot(ip - p1) < 0 || n2.dot(ip - p2) < 0 || n3.dot(ip - p3) < 0 || n4.dot(ip - p4) < 0) {
            intersection.doesIntersect = false;
            return intersection;
        }
        intersection.doesIntersect = true;
        intersection.t = t;
        intersection.ip = ip;
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
        Intersection intersection;
        point origin = ray.origin;
        point direction = ray.direction;
        point p1 = point(x, y, z);
        point p2 = point(x + side, y, z);
        point p3 = point(x + side, y + side, z);
        point p4 = point(x, y + side, z);
        point p5 = point(x + side, y, z + side);
        point p6 = point(x + side, y + side, z + side);
        point p7 = point(x, y + side, z + side);
        point p8 = point(x, y, z + side);
        point normal = (p2 - p1).cross(p3 - p1);
        double t = (p1 - origin).dot(normal) / direction.dot(normal);
        point ip = origin + direction * t;
        point v1 = p2 - p1;
        point v2 = p3 - p2;
        point v3 = p4 - p3;
        point v4 = p1 - p4;
        point v5 = p5 - p1;
        point v6 = p6 - p2;
        point v7 = p7 - p3;
        point v8 = p8 - p4;
        point v9 = p5 - p8;
        point v10 = p6 - p7;
        point v11 = p8 - p7;
        point v12 = p6 - p5;
        point n1 = v1.cross(v5);
        point n2 = v2.cross(v6);
        point n3 = v3.cross(v7);
        point n4 = v4.cross(v8);
        point n5 = v9.cross(v10);
        point n6 = v11.cross(v12);
        if (t < 0) {
            intersection.doesIntersect = false;
            return intersection;
        }
        if (n1.dot(ip - p1) < 0 || n2.dot(ip - p2) < 0 || n3.dot(ip - p3) < 0 || n4.dot(ip - p4) < 0 || n5.dot(ip - p8) < 0 || n6.dot(ip - p5) < 0) {
            intersection.doesIntersect = false;
            return intersection;
        }
        intersection.doesIntersect = true;
        intersection.t = t;
        intersection.ip = ip;
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
};

void readInputFile(string fileName) {
    ifstream fin = ifstream(fileName);
    string line;
    fin >> nearPlane >> farPlane >> fovY >> aspectRatio;
    fin >> levelOfRecursion;
    fin >> width;
    height = width;
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
            sphere.buildUnitPositiveX();
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
    double aspectRatio = (double)width / (double)height;
    double screenHeight = screenWidth / aspectRatio;  // Using aspect ratio

    point midpoint = pos + l * nearPlane;
    // Calculate the top-left corner of the screen in the 3D world
    point topLeft = midpoint + (u * (screenHeight / 2)) - (r * (screenWidth / 2));
    double du = screenWidth / width;
    double dv = screenHeight / height;
    topLeft = topLeft + r * (0.5 * du) - u * (0.5 * dv);

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            // Calculate the 3D point for this pixel in the OpenGL world
            pointBuffer[i][j] = topLeft + (r * i * du) - (u * j * dv);

            // Generate the ray for this pixel
            rays[i][j].origin = pos;
            point dir = pointBuffer[i][j] - pos;
            dir.normalize();
            rays[i][j].direction = dir;
        }
    }
}

void capture() {
    bitmap_image img(width, height);
    img.set_all_channels(0, 0, 0);
    // now populate pointBuffer with the color of the object that intersects witht the rays
    generateRays();
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            Ray ray = rays[i][j];
            double tMin = std::numeric_limits<double>::max();
            unsigned r = 0, g = 0, b = 0;
            for (const auto& sphere : spheres) {
                Intersection is = sphere.intersect(ray);
                if (!is.doesIntersect) continue;
                if (is.t < tMin) {
                    tMin = is.t;
                    r = sphere.r;
                    g = sphere.g;
                    b = sphere.b;
                }
            }
            // for (const auto& pyramid : pyramids) {
            //     Intersection is = pyramid.intersect(ray);
            //     if (!is.doesIntersect) continue;
            //     if (is.t < tMin) {
            //         tMin = is.t;
            //         r = pyramid.r;
            //         g = pyramid.g;
            //         b = pyramid.b;
            //     }
            // }
            // for (const auto& cube : cubes) {
            //     Intersection is = cube.intersect(ray);
            //     if (!is.doesIntersect) continue;
            //     if (is.t < tMin) {
            //         tMin = is.t;
            //         r = cube.r;
            //         g = cube.g;
            //         b = cube.b;
            //     }
            // }
            img.set_pixel(i, j, r * 255, g * 255, b * 255);
        }
    }
    img.save_image("out.bmp");
    img.clear();
}
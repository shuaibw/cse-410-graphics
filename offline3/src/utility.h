#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using namespace std;
// Function signatures
void readInputFile(string fileName);

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
};

class Sphere {
   public:
    double x, y, z, radius;
    double r, g, b;
    double ka, kd, ks, kr;
    double shine;
    // override cout
    friend ostream& operator<<(ostream& os, const Sphere& sphere) {
        os << "Sphere: " << sphere.x << " " << sphere.y << " " << sphere.z << " " << sphere.radius << endl;
        os << "Color: " << sphere.r << " " << sphere.g << " " << sphere.b << endl;
        os << "ka: " << sphere.ka << " kd: " << sphere.kd << " ks: " << sphere.ks << " kr: " << sphere.kr << endl;
        os << "shine: " << sphere.shine << endl;
        return os;
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

struct point {
    GLdouble x, y, z;
    point(double x, double y, double z)
        : x{x}, y{y}, z{z} {}
    point()
        : x{0}, y{0}, z{0} {}
    // overload - operator
    point operator-(const point& p) const {
        point temp;
        temp.x = x - p.x;
        temp.y = y - p.y;
        temp.z = z - p.z;
        return temp;
    }
    // overload + operator
    point operator+(const point& p) const {
        point temp;
        temp.x = x + p.x;
        temp.y = y + p.y;
        temp.z = z + p.z;
        return temp;
    }
    // overload = operator
    point operator=(const point& p) {
        x = p.x;
        y = p.y;
        z = p.z;
        return *this;
    }
    // overload * operator
    point operator*(const double& d) const {
        point temp;
        temp.x = x * d;
        temp.y = y * d;
        temp.z = z * d;
        return temp;
    }
    // normalize
    void normalize() {
        GLdouble scale = sqrt(x * x + y * y + z * z);
        if (scale >= 0.01) {
            x /= scale;
            y /= scale;
            z /= scale;
        }
    }
    point cross(const point& p) const {
        return point(y * p.z - z * p.y, z * p.x - x * p.z,
                     x * p.y - y * p.x);
    }
};

// Global variables
double nearPlane, farPlane, fovY, aspectRatio;
int levelOfRecursion, width, height;
int numObjects, numNormLights, numSpotLights;
vector<Sphere> spheres;
vector<Pyramid> pyramids;
vector<Cube> cubes;
vector<NormalLightSource> normLights;
vector<SpotLightSource> spotLights;
CheckerBoard cboard;

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
    for(const auto& sphere : spheres) {
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
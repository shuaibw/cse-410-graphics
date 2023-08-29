#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <GL/glut.h>
#include "bitmap_image.hpp"

using namespace std;
// Class signatures
class CheckerBoard;
class point;
class Ray;
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

class Intersection {
   public:
    bool doesIntersect;
    double t;
    point intersectionPoint;
    Intersection()
        : doesIntersect{false}, t{0}, intersectionPoint{} {}
    Intersection(bool doesIntersect, double t, point intersectionPoint)
        : doesIntersect{doesIntersect}, t{t}, intersectionPoint{intersectionPoint} {}
};
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
        point N = {0, 0, 1};  // Normal of the checkerboard (XY plane)
        point P = {0, 0, 0};  // A point on the checkerboard (XY plane)

        double denominator = N.dot(ray.direction);

        // If the ray is parallel to the checkerboard, no intersection
        if (fabs(denominator) < 1e-6) {
            return {false, 0, {}};
        }

        double t = (P - ray.origin).dot(N) / denominator;

        // If t < 0, the intersection is behind the ray's origin
        if (t < 0) {
            return {false, 0, {}};
        }

        point intersectionPoint = ray.origin + ray.direction * t;

        // Optional: Check if the intersection point is within a certain boundary of the checkerboard
        // if (fabs(intersectionPoint.x) > cboard.width || fabs(intersectionPoint.y) > cboard.width) {
        //     return {false, 0, {}};
        // }

        return {true, t, intersectionPoint};
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
    // generate vertices for +X face only by intersecting 2 circular planes
    // (longitudinal and latitudinal) at the given longitude/latitude angles
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
        point center = {x, y, y};
        point OC = ray.origin - center;
        double a = ray.direction.dot(ray.direction);
        double b = 2.0 * OC.dot(ray.direction);
        double c = OC.dot(OC) - r * r;

        double discriminant = b * b - 4 * a * c;
        if (discriminant < 0) {
            return {false, 0, {}};
        } else {
            double t1 = (-b - sqrt(discriminant)) / (2.0 * a);
            double t2 = (-b + sqrt(discriminant)) / (2.0 * a);
            if (t1 > 0) {
                return {true, t1, ray.origin + ray.direction * t1};
            } else if (t2 > 0) {
                return {true, t2, ray.origin + ray.direction * t2};
            } else {
                return {false, 0, {}};
            }
        }
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
        double tmin = (x - ray.origin.x) / ray.direction.x;
        double tmax = (x + side - ray.origin.x) / ray.direction.x;

        if (tmin > tmax) std::swap(tmin, tmax);

        double tymin = (y - ray.origin.y) / ray.direction.y;
        double tymax = (y + side - ray.origin.y) / ray.direction.y;

        if (tymin > tymax) std::swap(tymin, tymax);

        if ((tmin > tymax) || (tymin > tmax))
            return {false, 0, {}};

        if (tymin > tmin)
            tmin = tymin;

        if (tymax < tmax)
            tmax = tymax;

        double tzmin = (z - ray.origin.z) / ray.direction.z;
        double tzmax = (z + side - ray.origin.z) / ray.direction.z;

        if (tzmin > tzmax) std::swap(tzmin, tzmax);

        if ((tmin > tzmax) || (tzmin > tmax))
            return {false, 0, {}};

        if (tzmin > tmin)
            tmin = tzmin;

        if (tzmax < tmax)
            tmax = tzmax;

        if (tmin < 0) {
            if (tmax < 0) {
                return {false, 0, {}};
            } else {
                // We're inside the cube
                return {true, tmax, ray.origin + ray.direction * tmax};
            }
        }

        return {true, tmin, ray.origin + ray.direction * tmin};
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

class point {
   public:
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

class Ray {
   public:
    point origin;
    point direction;
    Ray(point origin, point direction)
        : origin{origin}, direction{direction} {}
    Ray()
        : origin{}, direction{} {}
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
    double screenWidth = 2 * nearPlane * tan(fovY / 2);
    double aspectRatio = (double)width / (double)height;
    double screenHeight = screenWidth / aspectRatio;  // Using aspect ratio

    point midpoint = pos + l * nearPlane;

    // Calculate the top-left corner of the screen in the 3D world
    point topLeft = midpoint + (u * (screenHeight / 2)) - (r * (screenWidth / 2));

    double pixelWidth = screenWidth / width;
    double pixelHeight = screenHeight / height;

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            // Calculate the 3D point for this pixel in the OpenGL world
            pointBuffer[i][j] = topLeft + (r * i * pixelWidth) - (u * j * pixelHeight);

            // Generate the ray for this pixel
            rays[i][j].origin = pos;
            point dir = pointBuffer[i][j] - pos;
            dir.normalize();
            rays[i][j].direction = dir;
        }
    }
}
void capture() {
    bitmap_image image(width, height);
    image.set_all_channels(0, 0, 0);  // set background to black

    double planeDistance = (height / 2.0) / tan(fovY / 2.0);
    point topLeft = pos + l * planeDistance - r * (width / 2) + u * (height / 2);
    topLeft = topLeft + r * 0.5 - u * 0.5;

    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            point curPixel = topLeft + r * i - u * j;
            Ray ray;
            ray.origin = pos;
            ray.direction = (curPixel - pos).normalize();

            Intersection nearestIntersection;
            nearestIntersection.t = std::numeric_limits<double>::max();

            // Check intersection with each Sphere
            for (const auto& sphere : spheres) {
                Intersection intersection = sphere.intersect(ray);
                if (intersection.doesIntersect && intersection.t < nearestIntersection.t) {
                    nearestIntersection = intersection;
                }
            }

            // Check intersection with each Cube
            for (const auto& cube : cubes) {
                Intersection intersection = cube.intersect(ray);
                if (intersection.doesIntersect && intersection.t < nearestIntersection.t) {
                    nearestIntersection = intersection;
                }
            }

            // Check intersection with each Pyramid
            // for (const auto& pyramid : pyramids) {
            //     Intersection intersection = pyramid.intersect(ray);
            //     if (intersection.doesIntersect && intersection.t < nearestIntersection.t) {
            //         nearestIntersection = intersection;
            //     }
            // }

            // Check intersection with Checkerboard
            Intersection intersection = cboard.intersect(ray);
            if (intersection.doesIntersect && intersection.t < nearestIntersection.t) {
                nearestIntersection = intersection;
            }

            if (nearestIntersection.doesIntersect) {
                // image.set_pixel(i, j, nearestIntersection.color.r, nearestIntersection.color.g, nearestIntersection.color.b);
            }
        }
    }

    image.save_image("out.bmp");
}

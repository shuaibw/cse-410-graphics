#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stack>
#include <vector>
#include "util.cpp"
#include "bitmap_image.hpp"
using namespace std;
bool debug = 0;
void readPoints(vector<triangle>& triangles, string input) {
    ifstream fin = ifstream(input);
    if (!fin.is_open()) {
        cout << "Error opening file: " << input << endl;
        return;
    }
    string line{};
    while (getline(fin, line)) {
        if (line == "") continue;
        point p1, p2, p3;
        stringstream ss(line);
        ss >> p1.x >> p1.y >> p1.z;
        getline(fin, line);
        stringstream ss2(line);
        ss2 >> p2.x >> p2.y >> p2.z;
        getline(fin, line);
        stringstream ss3(line);
        ss3 >> p3.x >> p3.y >> p3.z;
        triangles.push_back(triangle(p1, p2, p3));
    }
    fin.close();
}
int main(int argc, char** argv) {
    if (argc != 2) {
        cout << "Usage: main.exe <test_case>" << endl;
        return 1;
    }
    // Instructions for autograder: g++ -std=c++17 main.cpp -o main.exe
    // Input file should be in ../tests/<test_case>/scene.txt
    // Provide <test_case> as command line argument, i.e. main.exe 1
    // Output files will be in the same directory as main.exe
    string input = "../tests_2/" + string(argv[1]) + "/scene.txt";
    string config = "../tests_2/" + string(argv[1]) + "/config.txt";
    ifstream fin(input);
    ofstream stage1("stage1.txt");
    ofstream stage2("stage2.txt");
    ofstream stage3("stage3.txt");
    stage1 << fixed << setprecision(7);
    stage2 << fixed << setprecision(7);
    stage3 << fixed << setprecision(7);
    cout << fixed << setprecision(7);
    if (!fin.is_open()) {
        cout << "Error opening file: " << input << endl;
        return 1;
    }
    // Initialize stack of transformation matrices
    stack<mat4> s;
    mat4 identityMat{};
    identityMat.setIdentity();
    s.push(identityMat);  // initial identity matrix pushed to stack
    // this will allow us to reuse the stack top as current transformation matrix

    // Initialize vector for storing triangles
    vector<triangle> triangles;

    // Read in initial values
    point eye, look, up;
    double fovY, aspectRatio, near, far;
    fin >> eye.x >> eye.y >> eye.z;
    fin >> look.x >> look.y >> look.z;
    fin >> up.x >> up.y >> up.z;
    fin >> fovY >> aspectRatio >> near >> far;

    // Generate view and projection matrices for stage 2 and 3
    string line{};
    mat4 V = getViewTransformationMatrix(eye, look, up);
    mat4 P = getProjectionMatrix(fovY, aspectRatio, near, far);

    while (!fin.eof()) {
        fin >> line;
        if (line == "triangle") {
            point p1, p2, p3;
            fin >> p1.x >> p1.y >> p1.z;
            fin >> p2.x >> p2.y >> p2.z;
            fin >> p3.x >> p3.y >> p3.z;
            // stage 1: Modeling transformation
            p1 = transformPoint(s.top(), p1);
            p2 = transformPoint(s.top(), p2);
            p3 = transformPoint(s.top(), p3);
            stage1 << p1.x << " " << p1.y << " " << p1.z << endl;
            stage1 << p2.x << " " << p2.y << " " << p2.z << endl;
            stage1 << p3.x << " " << p3.y << " " << p3.z << endl
                   << endl;
            // stage 2: View transformation
            p1 = transformPoint(V, p1);
            p2 = transformPoint(V, p2);
            p3 = transformPoint(V, p3);
            stage2 << p1.x << " " << p1.y << " " << p1.z << endl;
            stage2 << p2.x << " " << p2.y << " " << p2.z << endl;
            stage2 << p3.x << " " << p3.y << " " << p3.z << endl
                   << endl;
            // stage 3: Projection transformation
            p1 = transformPoint(P, p1);
            p2 = transformPoint(P, p2);
            p3 = transformPoint(P, p3);
            stage3 << p1.x << " " << p1.y << " " << p1.z << endl;
            stage3 << p2.x << " " << p2.y << " " << p2.z << endl;
            stage3 << p3.x << " " << p3.y << " " << p3.z << endl
                   << endl;
            // Store triangle for stage 4
            if (!debug) {
                triangle t{p1, p2, p3};
                triangles.push_back(t);
            }
        } else if (line == "translate") {
            point translate;
            fin >> translate.x >> translate.y >> translate.z;
            mat4 translateMat = getTranslationMatrix(translate);
            s.top() = s.top() * translateMat;

        } else if (line == "scale") {
            point scale;
            fin >> scale.x >> scale.y >> scale.z;
            mat4 scaleMat = getScaleMatrix(scale);
            s.top() = s.top() * scaleMat;
        } else if (line == "rotate") {
            double theta;
            point a;
            fin >> theta >> a.x >> a.y >> a.z;
            mat4 rotateMat = getRotationMatrix(a, theta);
            s.top() = s.top() * rotateMat;
        } else if (line == "push") {
            s.push(s.top());
        } else if (line == "pop") {
            s.pop();
        } else if (line == "end") {
            break;
        }
    }
    fin.close();
    stage1.close();
    stage2.close();
    stage3.close();
    if (debug) {
        readPoints(triangles, "debug3.txt");
    }
    fin = ifstream(config);
    if (!fin.is_open()) {
        cout << "Error opening file: " << config << endl;
        return 1;
    }
    // Read in screen width and height
    int width, height;
    fin >> width >> height;
    fin.close();
    // Stage 4: Clipping and scan conversion using Z-buffer algorithm
    const double dx = 2.0 / width;
    const double dy = 2.0 / height;
    // Scan starts from top left corner
    const double topY = 1.0 - dy / 2;
    const double leftX = -1.0 + dx / 2;
    // Z-buffer initialization
    const double zNear = -1;
    const double zFar = 1;
    const double zMax = 2;
    std::vector<std::vector<double>> zBuffer(width, std::vector<double>(height));
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            zBuffer[i][j] = zMax;
        }
    }
    bitmap_image image(width, height);
    image.set_all_channels(0, 0, 0);  // set background to black
    for (const auto& t : triangles) {
        // Calulate topmost and bottommost rows of the triangle
        double top = t.p1.y;  // p1 is the topmost point
        double bottom = min(t.p2.y, t.p3.y);
        int topRow = (topY - top) / dy;
        int bottomRow = (topY - bottom) / dy;
        if (topRow < 0) topRow = 0;
        if (bottomRow >= height) bottomRow = height - 1;
        for (int i = topRow; i <= bottomRow; i++) {
            double y = topY - i * dy;
            if (y > top + dy / 2) continue;
            pair<pair<double, double>, pair<double, double>> xs = intersectionX(t, y);
            double left = xs.first.first;
            double zLeft = xs.first.second;
            double right = xs.second.first;
            double zRight = xs.second.second;
            int leftCol = (left - leftX + dx/2) / dx;
            int rightCol = (right - leftX) / dx;
            if (leftCol < 0) leftCol = 0;
            if (rightCol >= width) rightCol = width - 1;
            for (int j = leftCol; j <= rightCol; j++) {
                double x = leftX + j * dx;
                double z = zLeft + (zRight - zLeft) * (x - left) / (right - left);
                if (z < zNear || z > zFar) continue;
                if (z < zBuffer[i][j]) {
                    zBuffer[i][j] = z;
                    rgb_t color = t.color;
                    image.set_pixel(j, i, color.red, color.green, color.blue);
                }
            }
        }
    }
    image.save_image("output.bmp");
    ofstream zBufferOut("z_buffer.txt");
    zBufferOut << fixed << setprecision(6);
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            if (zBuffer[i][j] < zMax) zBufferOut << zBuffer[i][j] << "\t";
        }
        zBufferOut << endl;
    }
    zBufferOut.close();
    return 0;
}

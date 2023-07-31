#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stack>
#include "util.cpp"
using namespace std;

int main(int argc, char** argv) {
    if (argc != 2) {
        cout << "Usage: main.exe <test_case>" << endl;
        return 1;
    }
    // Instructions for autograder: g++ -std=c++17 main.cpp -o main.exe
    // Input file should be in ../tests/<test_case>/scene.txt
    // Provide <test_case> as command line argument, i.e. main.exe 1
    // Output files will be in the same directory as main.exe
    string input = "../tests_1/" + string(argv[1]) + "/scene.txt";
    ifstream fin(input);
    ofstream stage1("stage1.txt");
    ofstream stage2("stage2.txt");
    ofstream stage3("stage3.txt");
    stage1 << fixed << setprecision(7);
    stage2 << fixed << setprecision(7);
    stage3 << fixed << setprecision(7);
    cout << fixed << setprecision(2);
    if (!fin.is_open()) {
        cout << "Error opening file: " << input << endl;
        return 1;
    }
    stack<mat4> s;
    mat4 identityMat{};
    identityMat.setIdentity();
    s.push(identityMat);  // initial identity matrix pushed to stack
    // this will allow us to reuse the stack top as current transformation matrix

    point eye, look, up;
    double fovY, aspectRatio, near, far;
    fin >> eye.x >> eye.y >> eye.z;
    fin >> look.x >> look.y >> look.z;
    fin >> up.x >> up.y >> up.z;
    fin >> fovY >> aspectRatio >> near >> far;

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
    return 0;
}
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stack>
#include "util.cpp"
using namespace std;

void transformPoints(const mat4& M, string input, string output) {
    ifstream fin = ifstream(input);
    ofstream fout = ofstream(output);
    fout << fixed << setprecision(7);
    if (!fin.is_open()) {
        cout << "Error opening file" << endl;
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
        p1 = transformPoint(M, p1);
        p2 = transformPoint(M, p2);
        p3 = transformPoint(M, p3);
        fout << p1.x << " " << p1.y << " " << p1.z << endl;
        fout << p2.x << " " << p2.y << " " << p2.z << endl;
        fout << p3.x << " " << p3.y << " " << p3.z << endl
             << endl;
    }
}

int main(int argc, char** argv) {
    if (argc != 2) {
        cout << "Usage: main.exe <test_case>" << endl;
        return 1;
    }
    string input = "../tests/" + string(argv[1]) + "/scene.txt";
    ifstream fin(input);
    ofstream fout("stage1.txt");
    fout << fixed << setprecision(7);
    cout << fixed << setprecision(2);
    if (!fin.is_open()) {
        cout << "Error opening file" << endl;
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

    // stage 1: Modeling transformation
    string line{};
    while (!fin.eof()) {
        fin >> line;
        if (line == "triangle") {
            point p1, p2, p3;
            fin >> p1.x >> p1.y >> p1.z;
            fin >> p2.x >> p2.y >> p2.z;
            fin >> p3.x >> p3.y >> p3.z;
            p1 = transformPoint(s.top(), p1);
            p2 = transformPoint(s.top(), p2);
            p3 = transformPoint(s.top(), p3);
            fout << p1.x << " " << p1.y << " " << p1.z << endl;
            fout << p2.x << " " << p2.y << " " << p2.z << endl;
            fout << p3.x << " " << p3.y << " " << p3.z << endl
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
    fout.close();

    // stage 2: View transformation
    mat4 V = getViewTransformationMatrix(eye, look, up);
    transformPoints(V, "stage1.txt", "stage2.txt");
    // stage 3: Projection transformation
    mat4 P = getProjectionMatrix(fovY, aspectRatio, near, far);
    transformPoints(P, "stage2.txt", "stage3.txt");
    return 0;
}
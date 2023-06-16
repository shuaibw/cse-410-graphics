#pragma once
#include <GL/glut.h>
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;
void drawAxes();
void drawCube();
void drawPyramid();
void drawCylinder();
void drawTriangle();
void drawSphere();
void drawPositiveX();
void drawHalfOctahedron();
void drawOctahedron();
void drawCurvedSurface();
void transformCylinder();
void fillUpperCyl();
void fillPosXCyl();
vector<float> buildUnitPositiveX();
// Global variables for triangle
float centroidX;
float centroidY;
float centroidZ;
float shrinkFactor;

// Global variables for sphere
vector<float> verticesXPos;
int subdivision;
float radius;

// Global variables for cylinder
float cr;  // cylinder radius
float ch;  // cylinder height

/* Draw axes: X in Red, Y in Green and Z in Blue */
void drawAxes() {
    glLineWidth(3);
    glBegin(GL_LINES);
    glColor3f(1, 0, 0);  // Red
    // X axis
    glVertex3f(0, 0, 0);
    glVertex3f(2, 0, 0);

    glColor3f(0, 1, 0);  // Green
    // Y axis
    glVertex3f(0, 0, 0);
    glVertex3f(0, 2, 0);

    glColor3f(0, 0, 1);  // Blue
    // Z axis
    glVertex3f(0, 0, 0);
    glVertex3f(0, 0, 2);
    glEnd();
}

/* Draw a cube centered at the origin */
void drawCube() {
    glBegin(GL_QUADS);  // Begin drawing the color cube with 6 quads
    // Top face (y = 1.0f)
    // Define vertices in counter-clockwise (CCW) order with normal pointing out
    glColor3f(0.0f, 1.0f, 0.0f);  // Green
    glVertex3f(1.0f, 1.0f, -1.0f);
    glVertex3f(-1.0f, 1.0f, -1.0f);
    glVertex3f(-1.0f, 1.0f, 1.0f);
    glVertex3f(1.0f, 1.0f, 1.0f);

    // Bottom face (y = -1.0f)
    glColor3f(1.0f, 0.5f, 0.0f);  // Orange
    glVertex3f(1.0f, -1.0f, 1.0f);
    glVertex3f(-1.0f, -1.0f, 1.0f);
    glVertex3f(-1.0f, -1.0f, -1.0f);
    glVertex3f(1.0f, -1.0f, -1.0f);

    // Front face  (z = 1.0f)
    glColor3f(1.0f, 0.0f, 0.0f);  // Red
    glVertex3f(1.0f, 1.0f, 1.0f);
    glVertex3f(-1.0f, 1.0f, 1.0f);
    glVertex3f(-1.0f, -1.0f, 1.0f);
    glVertex3f(1.0f, -1.0f, 1.0f);

    // Back face (z = -1.0f)
    glColor3f(1.0f, 1.0f, 0.0f);  // Yellow
    glVertex3f(1.0f, -1.0f, -1.0f);
    glVertex3f(-1.0f, -1.0f, -1.0f);
    glVertex3f(-1.0f, 1.0f, -1.0f);
    glVertex3f(1.0f, 1.0f, -1.0f);

    // Left face (x = -1.0f)
    glColor3f(0.0f, 0.0f, 1.0f);  // Blue
    glVertex3f(-1.0f, 1.0f, 1.0f);
    glVertex3f(-1.0f, 1.0f, -1.0f);
    glVertex3f(-1.0f, -1.0f, -1.0f);
    glVertex3f(-1.0f, -1.0f, 1.0f);

    // Right face (x = 1.0f)
    glColor3f(1.0f, 0.0f, 1.0f);  // Magenta
    glVertex3f(1.0f, 1.0f, -1.0f);
    glVertex3f(1.0f, 1.0f, 1.0f);
    glVertex3f(1.0f, -1.0f, 1.0f);
    glVertex3f(1.0f, -1.0f, -1.0f);
    glEnd();  // End of drawing color-cube
}

/* Draw a pyramid centered at the origin */
void drawPyramid() {
    glBegin(GL_TRIANGLES);  // Begin drawing the pyramid with 4 triangles
    // Front
    glColor3f(1.0f, 0.0f, 0.0f);  // Red
    glVertex3f(0.0f, 1.0f, 0.0f);
    glColor3f(0.0f, 1.0f, 0.0f);  // Green
    glVertex3f(-1.0f, -1.0f, 1.0f);
    glColor3f(0.0f, 0.0f, 1.0f);  // Blue
    glVertex3f(1.0f, -1.0f, 1.0f);

    // Right
    glColor3f(1.0f, 0.0f, 0.0f);  // Red
    glVertex3f(0.0f, 1.0f, 0.0f);
    glColor3f(0.0f, 0.0f, 1.0f);  // Blue
    glVertex3f(1.0f, -1.0f, 1.0f);
    glColor3f(0.0f, 1.0f, 0.0f);  // Green
    glVertex3f(1.0f, -1.0f, -1.0f);

    // Back
    glColor3f(1.0f, 0.0f, 0.0f);  // Red
    glVertex3f(0.0f, 1.0f, 0.0f);
    glColor3f(0.0f, 1.0f, 0.0f);  // Green
    glVertex3f(1.0f, -1.0f, -1.0f);
    glColor3f(0.0f, 0.0f, 1.0f);  // Blue
    glVertex3f(-1.0f, -1.0f, -1.0f);

    // Left
    glColor3f(1.0f, 0.0f, 0.0f);  // Red
    glVertex3f(0.0f, 1.0f, 0.0f);
    glColor3f(0.0f, 0.0f, 1.0f);  // Blue
    glVertex3f(-1.0f, -1.0f, -1.0f);
    glColor3f(0.0f, 1.0f, 0.0f);  // Green
    glVertex3f(-1.0f, -1.0f, 1.0f);
    glEnd();  // Done drawing the pyramid
}

void drawTriangle() {
    glBegin(GL_TRIANGLES);
    glVertex3f((1.0 - centroidX) * shrinkFactor + centroidX, (0.0 - centroidY) * shrinkFactor + centroidY, (0.0 - centroidZ) * shrinkFactor + centroidZ);
    glVertex3f((0.0 - centroidX) * shrinkFactor + centroidX, (1.0 - centroidY) * shrinkFactor + centroidY, (0.0 - centroidZ) * shrinkFactor + centroidZ);
    glVertex3f((0.0 - centroidX) * shrinkFactor + centroidX, (0.0 - centroidY) * shrinkFactor + centroidY, (1.0 - centroidZ) * shrinkFactor + centroidZ);
    glEnd();
}
void drawHalfOctahedron() {
    glColor3f(0, 1, 1);  // cyan
    drawTriangle();
    glRotated(90, 0, 1, 0);
    glColor3f(1, 0, 1);  // magenta
    drawTriangle();
    glRotated(90, 0, 1, 0);
    glColor3f(0, 1, 1);  // cyan
    drawTriangle();
    glRotated(90, 0, 1, 0);
    glColor3f(1, 0, 1);  // magenta
    drawTriangle();
    glRotated(90, 0, 1, 0);
}
void drawOctahedron() {
    drawHalfOctahedron();
    glPushMatrix();
    glRotated(90, 0, 1, 0);
    glScalef(1, -1, 1);  // mirror along XZ plane
    drawHalfOctahedron();
    glPopMatrix();
}
void drawSphere() {
    glColor3f(0, 1, 0);
    glPushMatrix();
    glTranslated(shrinkFactor, 0, 0);
    glScaled(1 - shrinkFactor, 1 - shrinkFactor, 1 - shrinkFactor);
    drawPositiveX();  // +X
    glPopMatrix();

    glColor3f(0, 0, 1);
    glPushMatrix();
    glRotated(90, 0, 1, 0);
    glTranslated(shrinkFactor, 0, 0);
    glScaled(1 - shrinkFactor, 1 - shrinkFactor, 1 - shrinkFactor);
    drawPositiveX();  // -Z
    glPopMatrix();

    glColor3f(0, 1, 0);
    glPushMatrix();
    glRotated(180, 0, 1, 0);
    glTranslated(shrinkFactor, 0, 0);
    glScaled(1 - shrinkFactor, 1 - shrinkFactor, 1 - shrinkFactor);
    drawPositiveX();  // -X
    glPopMatrix();

    glColor3f(0, 0, 1);
    glPushMatrix();
    glRotated(-90, 0, 1, 0);
    glTranslated(shrinkFactor, 0, 0);
    glScaled(1 - shrinkFactor, 1 - shrinkFactor, 1 - shrinkFactor);
    drawPositiveX();  // +Z
    glPopMatrix();

    glPushMatrix();
    glColor3f(1, 0, 0);
    glRotated(90, 0, 0, 1);
    glTranslated(shrinkFactor, 0, 0);
    glScaled(1 - shrinkFactor, 1 - shrinkFactor, 1 - shrinkFactor);
    drawPositiveX();  // +Y
    glPopMatrix();

    glPushMatrix();
    glRotated(-90, 0, 0, 1);
    glTranslated(shrinkFactor, 0, 0);
    glScaled(1 - shrinkFactor, 1 - shrinkFactor, 1 - shrinkFactor);
    drawPositiveX();  // -Y
    glPopMatrix();
}
void drawPositiveX() {
    // draw the filled sphere from vertices
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
    // glPointSize(3);
    // float r = 1;
    // float g = 0;
    // float b = 0;
    // glBegin(GL_POINTS);
    // for (int i = 0; i < vertices.size(); i += 3) {
    //     // change color every 2 vertices
    //     if (i % 6 == 0) {
    //         if (r == 1) {
    //             r = 0;
    //             g = 1;
    //             b = 0;
    //         } else if (g == 1) {
    //             r = 0;
    //             g = 0;
    //             b = 1;
    //         } else {
    //             r = 1;
    //             g = 0;
    //             b = 0;
    //         }
    //         glColor3f(r, g, b);
    //     }
    //     if (i == 0) glColor3f(1, 1, 0);
    //     glVertex3f(vertices[i], vertices[i + 1], vertices[i + 2]);
    // }
    // glEnd();
}
void drawCylinder() {
    fillUpperCyl();
    glPushMatrix();
    glScaled(1, -1, 1);  // mirror along XZ plane
    fillUpperCyl();
    glPopMatrix();

    fillPosXCyl();
    glPushMatrix();
    glScaled(-1, 1, 1); // mirror along YZ plane
    fillPosXCyl();
    glPopMatrix();
}
void fillUpperCyl(){
    transformCylinder();
    glRotated(90, 0, 1, 0);
    transformCylinder();
    glRotated(90, 0, 1, 0);
    transformCylinder();
    glRotated(90, 0, 1, 0);
    transformCylinder();
    glRotated(90, 0, 1, 0);
}
void fillPosXCyl(){
    glPushMatrix();
    glRotated(90, 1, 0 ,0 );
    transformCylinder();
    glPopMatrix();
    glPushMatrix();
    glRotated(-90, 1, 0 ,0 );
    transformCylinder();
    glPopMatrix();
}
void transformCylinder() {
    cr = sqrt(1.0 / 3.0) * (1 - shrinkFactor);
    ch = sqrt(2) * shrinkFactor;
    glPushMatrix();
    glRotated(45, 0, 0, 1);
    glTranslated(shrinkFactor / sqrt(2), 0, 0);
    drawCurvedSurface();
    glPopMatrix();
}
void drawCurvedSurface() {
    float dTheta = 5.0;
    float PI = acos(-1);
    for (int theta = -35; theta < 35.2; theta += dTheta) {
        if (theta % 10 == 0)
            glColor3f(1, 1, 1);
        else
            glColor3f(1, 0.59, 0.67);
        float start = theta * PI / 180.0;
        float end = (theta + dTheta) * PI / 180.0;
        glBegin(GL_QUADS);
        glVertex3f(cr * cos(start), ch / 2, cr * sin(start));
        glVertex3f(cr * cos(end), ch / 2, cr * sin(end));
        glVertex3f(cr * cos(end), -ch / 2, cr * sin(end));
        glVertex3f(cr * cos(start), -ch / 2, cr * sin(start));
        glEnd();
    }
    float theta = -35;
    dTheta = -0.2644;
    float start = theta * PI / 180.0;
    float end = (theta + dTheta) * PI / 180.0;
    glBegin(GL_QUADS);
    glVertex3f(cr * cos(start), ch / 2, cr * sin(start));
    glVertex3f(cr * cos(end), ch / 2, cr * sin(end));
    glVertex3f(cr * cos(end), -ch / 2, cr * sin(end));
    glVertex3f(cr * cos(start), -ch / 2, cr * sin(start));
    glEnd();
    theta = 35;
    dTheta = 0.2644;
    start = theta * PI / 180.0;
    end = (theta + dTheta) * PI / 180.0;
    glBegin(GL_QUADS);
    glVertex3f(cr * cos(start), ch / 2, cr * sin(start));
    glVertex3f(cr * cos(end), ch / 2, cr * sin(end));
    glVertex3f(cr * cos(end), -ch / 2, cr * sin(end));
    glVertex3f(cr * cos(start), -ch / 2, cr * sin(start));
    glEnd();
}
// generate vertices for +X face only by intersecting 2 circular planes
// (longitudinal and latitudinal) at the given longitude/latitude angles
vector<float> buildUnitPositiveX() {
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

    return vertices;
}
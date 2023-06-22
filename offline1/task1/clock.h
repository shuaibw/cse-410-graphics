#pragma once
#include <GL/glut.h>
#include <math.h>

void line(double x1, double y1, double x2, double y2);
void circle(double x, double y, double r, int slices);
void filled_circle(double x, double y, double r, int slices);
void text(GLdouble x, GLdouble y, char* str, void* font);
void setcolor(double r, double g, double b);
void filled_triangle(double x1, double y1, double x2, double y2, double x3, double y3);

void line(double x1, double y1, double x2, double y2) {
    glBegin(GL_LINE_STRIP);
    glVertex2f(x1, y1);
    glVertex2f(x2, y2);
    glEnd();
}

void circle(double x, double y, double r, int slices = 100) {
    double t, PI = acos(-1.0), dt, x1, y1, xp, yp;
    dt = 2 * PI / slices;
    xp = x + r;
    yp = y;
    for (t = 0; t <= 2 * PI; t += dt) {
        x1 = x + r * cos(t);
        y1 = y + r * sin(t);
        line(xp, yp, x1, y1);
        xp = x1;
        yp = y1;
    }
}

void filled_circle(double x, double y, double r, int slices = 100) {
    double t, PI = acos(-1.0), dt, x1, y1, xp, yp;
    dt = 2 * PI / slices;
    xp = x + r;
    yp = y;
    glBegin(GL_POLYGON);
    for (t = 0; t <= 2 * PI; t += dt) {
        x1 = x + r * cos(t);
        y1 = y + r * sin(t);

        glVertex2f(xp, yp);
        xp = x1;
        yp = y1;
    }
    glEnd();
}

void text(GLdouble x, GLdouble y, char* str, void* font = GLUT_BITMAP_8_BY_13) {
    glRasterPos3d(x, y, 0);
    int i;
    for (i = 0; str[i]; i++) {
        glutBitmapCharacter(font, str[i]);  //,GLUT_BITMAP_8_BY_13, GLUT_BITMAP_TIMES_ROMAN_24
    }
}
void setcolor(double r, double g, double b) {
    double mmx;
    mmx = r;
    if (g > mmx) mmx = g;
    if (b > mmx) mmx = b;
    mmx = 255;
    if (mmx > 0) {
        r /= mmx;
        g /= mmx;
        b /= mmx;
    }
    glColor3f(r, g, b);
}

void filled_triangle(double x1, double y1, double x2, double y2, double x3, double y3) {
    glBegin(GL_POLYGON);
    glVertex2f(x1, y1);
    glVertex2f(x2, y2);
    glVertex2f(x3, y3);
    glEnd();
}
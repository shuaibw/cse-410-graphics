// #include <windows.h>  // for MS Windows
#include <GL/glut.h>  // GLUT, include glu.h and gl.h
#include <cmath>
#include "magic_cube.h"
#include <iostream>
using namespace std;
void handleFixedUpwardRotation();
void handleFixedDownwardRotation();
/* Initialize OpenGL Graphics */
void initGL() {
    // Set "clearing" or background color
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);  // Black and opaque
    glEnable(GL_DEPTH_TEST);               // Enable depth testing for z-culling
}
// Global variables
struct point {
    GLfloat x, y, z;
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
};
struct point pos;  // position of the eye
struct point l;    // look/forward direction
struct point r;    // right direction
struct point u;    // up direction
bool isAxes = true;

/*  Handler for window-repaint event. Call back when the window first appears and
    whenever the window needs to be re-painted. */
void display() {
    // glClear(GL_COLOR_BUFFER_BIT);            // Clear the color buffer (background)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);  // To operate on Model-View matrix
    glLoadIdentity();            // Reset the model-view matrix

    // default arguments of gluLookAt
    // gluLookAt(0,0,0, 0,0,-100, 0,1,0);

    // control viewing (or camera)
    gluLookAt(pos.x, pos.y, pos.z,
              pos.x + l.x, pos.y + l.y, pos.z + l.z,
              u.x, u.y, u.z);
    // draw
    if (isAxes)
        drawAxes();
    glRotatef(angle, 0.0f, 1.0f, 0.0f);
    drawOctahedron();
    drawSphere();
    drawCylinder();

    glutSwapBuffers();  // Render now
}

void reshapeListener(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
    if (height == 0)
        height = 1;  // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
    glLoadIdentity();             // Reset the projection matrix
    gluPerspective(45.0f, aspect, 0.1f, 100.0f);
}

float angleBetween(point p, float dy) {
    point tp = p - l;
    point tq = tp;
    tq.y += dy;
    float dot = tp.x * tq.x + tp.y * tq.y;
    float mag = sqrt(tp.x * tp.x + tp.y * tp.y) *
                sqrt(tq.x * tq.x + tq.y * tq.y);
    float angle = acos(dot / mag);
    return angle;
}
void rotateCameraLU(float angle) {
    l.x = l.x * cos(angle) + u.x * sin(angle);
    l.y = l.y * cos(angle) + u.y * sin(angle);
    l.z = l.z * cos(angle) + u.z * sin(angle);
    u.x = u.x * cos(angle) - l.x * sin(angle);
    u.y = u.y * cos(angle) - l.y * sin(angle);
    u.z = u.z * cos(angle) - l.z * sin(angle);
}
void handleFixedUpwardRotation() {
    float dy = 0.1;
    pos.y += dy;
    float angle = angleBetween(pos, dy);
    // look down by angle
    rotateCameraLU(-angle);
}
void handleFixedDownwardRotation() {
    float dy = 0.1;
    pos.y -= dy;
    float angle = angleBetween(pos, dy);
    // look up by angle
    rotateCameraLU(angle);
}
/* Callback handler for normal-key event */
void keyboardListener(unsigned char key, int xx, int yy) {
    double rate = 0.02;
    double v = 0.05;
    switch (key) {
            // Control shrinkFactor for triangle
        case ',':
            if (shrinkFactor > 0.0) shrinkFactor -= v;
            break;
        case '.':
            if (shrinkFactor < 1.0) shrinkFactor += v;
            break;
        case 'a':
            // rotate the object in clockwise direction
            angle -= 5;
            if (angle <= -360) angle = 0;
            break;
        case 'd':
            // rotate the object in anti-clockwise direction
            angle += 5;
            if (angle >= 360) angle = 0;
            break;
        case 'w':
            handleFixedUpwardRotation();
            break;
        case 's':
            handleFixedDownwardRotation();
            break;
        case '1':
            r.x = r.x * cos(rate) + l.x * sin(rate);
            r.y = r.y * cos(rate) + l.y * sin(rate);
            r.z = r.z * cos(rate) + l.z * sin(rate);

            l.x = l.x * cos(rate) - r.x * sin(rate);
            l.y = l.y * cos(rate) - r.y * sin(rate);
            l.z = l.z * cos(rate) - r.z * sin(rate);
            break;

        case '2':
            r.x = r.x * cos(-rate) + l.x * sin(-rate);
            r.y = r.y * cos(-rate) + l.y * sin(-rate);
            r.z = r.z * cos(-rate) + l.z * sin(-rate);

            l.x = l.x * cos(-rate) - r.x * sin(-rate);
            l.y = l.y * cos(-rate) - r.y * sin(-rate);
            l.z = l.z * cos(-rate) - r.z * sin(-rate);
            break;

        case '3':
            l.x = l.x * cos(rate) + u.x * sin(rate);
            l.y = l.y * cos(rate) + u.y * sin(rate);
            l.z = l.z * cos(rate) + u.z * sin(rate);

            u.x = u.x * cos(rate) - l.x * sin(rate);
            u.y = u.y * cos(rate) - l.y * sin(rate);
            u.z = u.z * cos(rate) - l.z * sin(rate);
            break;

        case '4':
            l.x = l.x * cos(-rate) + u.x * sin(-rate);
            l.y = l.y * cos(-rate) + u.y * sin(-rate);
            l.z = l.z * cos(-rate) + u.z * sin(-rate);

            u.x = u.x * cos(-rate) - l.x * sin(-rate);
            u.y = u.y * cos(-rate) - l.y * sin(-rate);
            u.z = u.z * cos(-rate) - l.z * sin(-rate);
            break;

        case '5':
            u.x = u.x * cos(rate) + r.x * sin(rate);
            u.y = u.y * cos(rate) + r.y * sin(rate);
            u.z = u.z * cos(rate) + r.z * sin(rate);

            r.x = r.x * cos(rate) - u.x * sin(rate);
            r.y = r.y * cos(rate) - u.y * sin(rate);
            r.z = r.z * cos(rate) - u.z * sin(rate);
            break;

        case '6':
            u.x = u.x * cos(-rate) + r.x * sin(-rate);
            u.y = u.y * cos(-rate) + r.y * sin(-rate);
            u.z = u.z * cos(-rate) + r.z * sin(-rate);

            r.x = r.x * cos(-rate) - u.x * sin(-rate);
            r.y = r.y * cos(-rate) - u.y * sin(-rate);
            r.z = r.z * cos(-rate) - u.z * sin(-rate);
            break;
        case 'p':
            // print all vectors
            cout << "--------------------------------------" << endl;
            cout << "pos: " << pos.x << " " << pos.y << " " << pos.z << endl;
            cout << "l: " << l.x << " " << l.y << " " << l.z << endl;
            cout << "r: " << r.x << " " << r.y << " " << r.z << endl;
            cout << "u: " << u.x << " " << u.y << " " << u.z << endl;
            cout << "angle: " << angle << endl;
            cout << "lambda: " << sqrt(u.x * u.x + u.y * u.y + u.z * u.z) << endl;
            break;

        default:
            break;
    }
    glutPostRedisplay();
}

/* Callback handler for special-key event */
void specialKeyListener(int key, int x, int y) {
    float rate = 0.1;
    switch (key) {
        case GLUT_KEY_UP:
            pos.x += l.x * rate;
            pos.y += l.y * rate;
            pos.z += l.z * rate;
            break;
        case GLUT_KEY_DOWN:
            pos.x -= l.x * rate;
            pos.y -= l.y * rate;
            pos.z -= l.z * rate;
            break;

        case GLUT_KEY_RIGHT:
            pos.x += r.x * rate;
            pos.y += r.y * rate;
            pos.z += r.z * rate;
            break;
        case GLUT_KEY_LEFT:
            pos.x -= r.x * rate;
            pos.y -= r.y * rate;
            pos.z -= r.z * rate;
            break;

        case GLUT_KEY_PAGE_UP:
            pos.x += u.x * rate;
            pos.y += u.y * rate;
            pos.z += u.z * rate;
            break;
        case GLUT_KEY_PAGE_DOWN:
            pos.x -= u.x * rate;
            pos.y -= u.y * rate;
            pos.z -= u.z * rate;
            break;

        case GLUT_KEY_INSERT:
            break;

        case GLUT_KEY_HOME:
            break;
        case GLUT_KEY_END:
            break;

        default:
            break;
    }
    glutPostRedisplay();
}
void initGlobalVars() {
    // Triangle
    centroidX = (1.0 + 0.0 + 0.0) / 3.0;
    centroidY = (0.0 + 1.0 + 0.0) / 3.0;
    centroidZ = (0.0 + 0.0 + 1.0) / 3.0;
    shrinkFactor = 1.0;
    // Sphere
    subdivision = 4;
    radius = sqrt(centroidX * centroidX + centroidY * centroidY + centroidZ * centroidZ);
    verticesXPos = buildUnitPositiveX();
    // cylinder
    cr = radius;
    ch = sqrt(2);
    // camera vectors
    pos.x = 4;
    pos.y = 4;
    pos.z = 4;
    l.x = -1 / sqrt(3);
    l.y = -1 / sqrt(3);
    l.z = -1 / sqrt(3);
    u.x = 0;
    u.y = 1;
    u.z = 0;
    r.x = l.y * u.z - l.z * u.y;
    r.y = l.z * u.x - l.x * u.z;
    r.z = l.x * u.y - l.y * u.x;
    float scale = sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
    r.x /= scale;
    r.y /= scale;
    r.z /= scale;
    u.x = r.y * l.z - r.z * l.y;
    u.y = r.z * l.x - r.x * l.z;
    u.z = r.x * l.y - r.y * l.x;
    scale = sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
    u.x /= scale;
    u.y /= scale;
    u.z /= scale;
}
/* Main function: GLUT runs as a console application starting at main()  */
int main(int argc, char** argv) {
    initGlobalVars();
    glutInit(&argc, argv);                                     // Initialize GLUT
    glutInitWindowSize(640, 640);                              // Set the window's initial width & height
    glutInitWindowPosition(50, 50);                            // Position the window's initial top-left corner
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);  // Depth, Double buffer, RGB color
    glutCreateWindow("OpenGL 3D Drawing");                     // Create a window with the given title
    glutDisplayFunc(display);                                  // Register display callback handler for window re-paint
    glutReshapeFunc(reshapeListener);                          // Register callback handler for window re-shape
    glutKeyboardFunc(keyboardListener);                        // Register callback handler for normal-key event
    glutSpecialFunc(specialKeyListener);                       // Register callback handler for special-key event
    initGL();                                                  // Our own OpenGL initialization
    glutMainLoop();                                            // Enter the event-processing loop
    return 0;
}

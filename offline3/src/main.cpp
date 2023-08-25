// #include <windows.h>  // for MS Windows
#include <GL/glut.h>  // GLUT, include glu.h and gl.h
#include <cmath>
#include "generate_object.h"
#include "utility.h"
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
    for (const auto& sphere : spheres) {
        drawSphere(sphere);
    }
    // drawCube();
    // drawSphere();
    // drawPyramid();

    glutSwapBuffers();  // Render now
}

void reshapeListener(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
    if (height == 0)
        height = 1;  // To prevent divide by 0
    GLdouble aspect = (GLdouble)width / (GLdouble)height;
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
    glLoadIdentity();             // Reset the projection matrix
    gluPerspective(fovY, aspect, nearPlane, farPlane);
}

double angleBetween(point p, double dy) {
    point tp = p - l;
    point tq = tp;
    tq.y += dy;
    double dot = tp.x * tq.x + tp.y * tq.y;
    double mag = sqrt(tp.x * tp.x + tp.y * tp.y) *
                sqrt(tq.x * tq.x + tq.y * tq.y);
    double angle = acos(dot / mag);
    return angle;
}
void rotateCameraLU(double angle) {
    l = l * cos(angle) + u * sin(angle);
    u = u * cos(angle) - l * sin(angle);
}
void handleFixedUpwardRotation() {
    double dy = 0.1;
    pos.y += dy;
    double angle = angleBetween(pos, dy);
    // look down by angle
    rotateCameraLU(-angle);
}
void handleFixedDownwardRotation() {
    double dy = 0.1;
    pos.y -= dy;
    double angle = angleBetween(pos, dy);
    // look up by angle
    rotateCameraLU(angle);
}
/* Callback handler for normal-key event */
void keyboardListener(unsigned char key, int xx, int yy) {
    double rate = 0.02;
    switch (key) {
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
            r = r * cos(rate) + l * sin(rate);
            l = l * cos(rate) - r * sin(rate);
            break;

        case '2':
            r = r * cos(-rate) + l * sin(-rate);
            l = l * cos(-rate) - r * sin(-rate);
            break;

        case '3':
            l = l * cos(rate) + u * sin(rate);
            u = u * cos(rate) - l * sin(rate);
            break;

        case '4':
            l = l * cos(-rate) + u * sin(-rate);
            u = u * cos(-rate) - l * sin(-rate);
            break;

        case '5':
            u = u * cos(rate) + r * sin(rate);
            r = r * cos(rate) - u * sin(rate);
            break;

        case '6':
            u = u * cos(-rate) + r * sin(-rate);
            r = r * cos(-rate) - u * sin(-rate);
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
    double rate = 0.5;
    switch (key) {
        case GLUT_KEY_UP:
            pos = pos + l * rate;
            break;
        case GLUT_KEY_DOWN:
            pos = pos - l * rate;
            break;

        case GLUT_KEY_RIGHT:
            pos = pos + r * rate;
            break;
        case GLUT_KEY_LEFT:
            pos = pos - r * rate;
            break;

        case GLUT_KEY_PAGE_UP:
            pos = pos + u * rate;
            break;
        case GLUT_KEY_PAGE_DOWN:
            pos = pos - u * rate;
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
    // Sphere
    subdivision = 4;
    radius = 5;
    verticesXPos = buildUnitPositiveX();
    // camera vectors
    pos = {40, 40, 40};
    l = {-1, -1, -1};
    u = {0, 1, 0};
    l.normalize();
    r = l.cross(u);
    r.normalize();
    u = r.cross(l);
    u.normalize();
}
/* Main function: GLUT runs as a console application starting at main()  */
int main(int argc, char** argv) {
    readInputFile("description.txt");
    initGlobalVars();
    glutInit(&argc, argv);                                     // Initialize GLUT
    glutInitWindowSize(width, height);                              // Set the window's initial width & height
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

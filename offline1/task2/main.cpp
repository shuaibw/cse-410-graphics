// #include <windows.h>  // for MS Windows
#include <GL/glut.h>  // GLUT, include glu.h and gl.h
#include <cmath>
#include "util.cpp"
#include <iostream>
using namespace std;
/* Initialize OpenGL Graphics */
void initGL() {
    // Set "clearing" or background color
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);  // Black and opaque
    glEnable(GL_DEPTH_TEST);               // Enable depth testing for z-culling
}
// Global variables
GLfloat eyex = 4, eyey = 4, eyez = 4;
GLfloat centerx = 0, centery = 0, centerz = 0;
GLfloat upx = 0, upy = 1, upz = 0;
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
    gluLookAt(eyex, eyey, eyez,
              centerx, centery, centerz,
              upx, upy, upz);
    // draw
    if (isAxes)
        drawAxes();
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

/* Callback handler for normal-key event */
void keyboardListener(unsigned char key, int x, int y) {
    float v = 0.05;
    switch (key) {
        // Control shrinkFactor for triangle
        case ',':
            if (shrinkFactor > 0.0) shrinkFactor -= v;
            break;
        case '.':
            if (shrinkFactor < 1.0) shrinkFactor += v;
            break;
        // Control eye (location of the eye)
        // control eyex
        case '1':
            eyex += v;
            break;
        case '2':
            eyex -= v;
            break;
        // control eyey
        case '3':
            eyey += v;
            break;
        case '4':
            eyey -= v;
            break;
        // control eyez
        case '5':
            eyez += v;
            break;
        case '6':
            eyez -= v;
            break;

        // Control center (location where the eye is looking at)
        // control centerx
        case 'q':
            centerx += v;
            break;
        case 'w':
            centerx -= v;
            break;
        // control centery
        case 'e':
            centery += v;
            break;
        case 'r':
            centery -= v;
            break;
        // control centerz
        case 't':
            centerz += v;
            break;
        case 'y':
            centerz -= v;
            break;

        // Control what is shown
        case 'a':
            isAxes = !isAxes;  // show/hide Axes if 'a' is pressed
            break;

        // Control exit
        case 27:      // ESC key
            exit(0);  // Exit window
            break;
    }
    glutPostRedisplay();  // Post a paint request to activate display()
}

/* Callback handler for special-key event */
void specialKeyListener(int key, int x, int y) {
    double v = 0.25;
    double lx = centerx - eyex;
    double lz = centerz - eyez;
    double s;
    switch (key) {
        case GLUT_KEY_LEFT:
            eyex += v * (upy * lz);
            eyez += v * (-lx * upy);
            s = sqrt(eyex * eyex + eyez * eyez) / (4 * sqrt(2));
            eyex /= s;
            eyez /= s;
            break;
        case GLUT_KEY_RIGHT:
            eyex += v * (-upy * lz);
            eyez += v * (lx * upy);
            s = sqrt(eyex * eyex + eyez * eyez) / (4 * sqrt(2));
            eyex /= s;
            eyez /= s;
            break;
        case GLUT_KEY_UP:
            eyey += v;
            break;
        case GLUT_KEY_DOWN:
            eyey -= v;
            break;

        default:
            return;
    }
    glutPostRedisplay();  // Post a paint request to activate display()
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

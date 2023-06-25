#include <iostream>
#include <GL/glut.h>
#include "clock.h"
#include <math.h>
#include <string.h>
#include <ctime>
#define pi 3.14159265359

int win_width = 800;
int win_height = 800;

unsigned int timer_cnt = 0;

int r = 120;  // Circle radius
int x = 400;  // Circle cx
int y = 500;  // Circle cy

float ang1 = pi / 4, ang2 = pi / 2;
float theta_max = pi / 4;

int i;

int s = 5;   // set seconds
int m = 30;  // set minutes
int h = 10;  // set hour

double stheta;
double mtheta;
double htheta;

char hour[3];
char min[3];
char sec[3];
char* AM = "AM";
char* buff;

double z = 0;  // angle
double d = 0;  // angle

double sec_x;
double sec_y;
double min_x;
double min_y;
double hour_x;
double hour_y;
double sec_width = 2;
double min_width = 4;
double hour_width = 6;
struct triangle {
    double x1, y1, x2, y2, x3, y3;
};
triangle sec_t, min_t, hour_t;

void showTime();
void hand_anim();
void pendulum_anim();
void timer(int value);

void reshape(GLsizei width, GLsizei height) {
    if (height == 0) height = 1;  // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
    glLoadIdentity();             // Reset the projection matrix
    if (width >= height) {
        gluOrtho2D(0, win_width * aspect, 0, win_height);
    } else {
        gluOrtho2D(0, win_width, 0, win_height / aspect);
    }
}
void set_clock_angles() {
    // initialize angles
    stheta = (pi * s / 30.0);
    mtheta = (pi * m / 30.0) + (pi * s / 1800.0);
    htheta = (pi * h / 6.0) + (pi * m / 360.0) + (pi * s / 21600.0);
}
void calculate_triangle_coords(struct triangle* t, char type, int offset, float thickness) {
    double midx, midy, theta;
    switch (type) {
        case 's':
            midx = sec_x;
            midy = sec_y;
            theta = d + stheta;
            break;
        case 'm':
            midx = min_x;
            midy = min_y;
            theta = d / 60 + mtheta;
            break;
        case 'h':
            midx = hour_x;
            midy = hour_y;
            theta = d / 720 + htheta;
            break;
        default:
            break;
    }
    thickness = thickness / 2;
    t->x1 = midx - thickness * cos(theta);
    t->y1 = midy + thickness * sin(theta);
    t->x2 = midx + thickness * cos(theta);
    t->y2 = midy - thickness * sin(theta);
    t->x3 = x + (r - offset) * sin(theta);
    t->y3 = y + (r - offset) * cos(theta);
};
void update_clock_coords() {
    int moff = 30;
    int hoff = 40;
    int soff = 20;
    int overlap = 4;
    sec_x = x + (r - soff + overlap) * sin(d + stheta);
    sec_y = y + (r - soff + overlap) * cos(d + stheta);
    min_x = x + (r - moff + overlap) * sin(d / 60 + mtheta);
    min_y = y + (r - moff + overlap) * cos(d / 60 + mtheta);
    hour_x = x + (r - hoff + overlap) * sin(d / 720 + htheta);
    hour_y = y + (r - hoff + overlap) * cos(d / 720 + htheta);
    calculate_triangle_coords(&sec_t, 's', 0, sec_width);
    calculate_triangle_coords(&min_t, 'm', 10, min_width);
    calculate_triangle_coords(&hour_t, 'h', 20, hour_width);
}
void hand_anim() {
    s++;
    m += s / 60;
    h += m / 60;
    if (s == 60) s = 0;
    if (m == 60) m = 0;
    if (h == 13) h = 1;
    d += pi / 30;  // 6 degrees for each second
    update_clock_coords();
}
void pendulum_anim() {
    ang2 += pi / 100.0;
    ang1 = theta_max * sin(ang2 + theta_max);
    // printf("%f\n", 180*ang2/pi);
}
void timer(int value) {
    pendulum_anim();
    glutTimerFunc(10, timer, 0);  // subsequent timer call at milliseconds
    timer_cnt++;
    if (timer_cnt % 100 == 0) {
        // std::cout << "Timer: " << timer_cnt << " counts" << std::endl;
        // show h,m,s,stheata,mtheta,htheta,d
        // std::cout << "h: " << h << " m: " << m << " s: " << s << std::endl;
        // std::cout << "stheta: " << stheta << " mtheta: " << mtheta << " htheta: " << htheta << std::endl;
        // std::cout << "d: " << d << std::endl;
        timer_cnt = 0;
        hand_anim();
    }
    glutPostRedisplay();  // Post a paint request to activate display()
}

void showTime() {
    // convert h1, m1, s1 to str
    snprintf(hour, 3, "%02d", h);
    snprintf(min, 3, "%02d", m);
    snprintf(sec, 3, "%02d", s);
    for (int i = 0; i < 2; i++) {
        buff[i] = hour[i];
        buff[i + 3] = min[i];
        buff[i + 6] = sec[i];
    }
    buff[2] = ':';
    buff[5] = ':';
    buff[8] = ' ';
    buff[9] = AM[0];
    buff[10] = AM[1];
    buff[11] = '\0';
    text(335, 650, buff, GLUT_BITMAP_TIMES_ROMAN_24);
}

void drawSineWave(int x1, int x2, int height) {
    float amplitude = 20.0;
    float frequency = 0.05;

    glBegin(GL_LINE_STRIP);
    glColor4f(1, 0.5451, 0.5804, 0.5);

    for (int x = x1; x < x2; ++x) {
        float y = amplitude * sin(x * frequency);
        glVertex2f(x, y + height);
    }

    glEnd();
}
void display() {
    glClear(GL_COLOR_BUFFER_BIT);  // Clear the color buffer
    glMatrixMode(GL_MODELVIEW);    // To operate on the model-view matrix
    glLoadIdentity();              // Reset model-view matrix
    z = 0;
    drawSineWave(150, 650, 750);
    line(153, 770, 153, 400);
    line(647, 770, 647, 400);
    line(153, 400, 400, 182);
    line(647, 400, 400, 182);
    setcolor(168, 230, 207);
    glLineWidth(4);
    circle(x, y, r);
    glLineWidth(2);
    // iCircle(x, y, r - 10);
    setcolor(168, 230, 207);
    filled_circle(x, y, 5);
    for (i = 0; i < 60; i++) {
        setcolor(185, 255, 229);
        line(x + (r - 10) * sin(z), y + (r - 10) * cos(z), x + r * sin(z), y + r * cos(z));
        if (!(i % 5)) {
            int hh = i / 5;
            if (hh == 0) hh = 12;
            buff = itoa(hh, buff, 10);
            int offset = 30;
            if (i > 30) offset = 25;
            setcolor(117, 223, 185);
            text(x + (r - offset) * sin(z), y + (r - offset) * cos(z), buff, GLUT_BITMAP_HELVETICA_12);
            setcolor(255, 214, 92);
            line(x + (r - 20) * sin(z), y + (r - 20) * cos(z), x + r * sin(z), y + r * cos(z));
        }
        z += pi / 30;
    }
    setcolor(0, 100, 255);
    glLineWidth(10);
    filled_circle(x, y - r, 5);
    line(x, y - r, x + r * sin(ang1), y - r * cos(ang1) - r);
    filled_circle(x + r * sin(ang1), y - r * cos(ang1) - r, 25);
    setcolor(97, 204, 255);
    filled_circle(x + r * sin(ang1), y - r * cos(ang1) - r, 10);

    setcolor(255, 139, 148);
    glLineWidth(sec_width);
    line(x, y, sec_x, sec_y);
    filled_triangle(sec_t.x1, sec_t.y1, sec_t.x2, sec_t.y2, sec_t.x3, sec_t.y3);

    setcolor(97, 204, 255);
    glLineWidth(min_width);
    line(x, y, min_x, min_y);
    filled_triangle(min_t.x1, min_t.y1, min_t.x2, min_t.y2, min_t.x3, min_t.y3);

    setcolor(255, 214, 92);
    glLineWidth(hour_width);
    line(x, y, hour_x, hour_y);
    filled_triangle(hour_t.x1, hour_t.y1, hour_t.x2, hour_t.y2, hour_t.x3, hour_t.y3);
    showTime();
    glutSwapBuffers();  // Swap front and back buffers (of double buffered mode)
}

/* Main function: GLUT runs as a console application starting at main() */
int main(int argc, char** argv) {
    buff = (char*)malloc(12 * sizeof(char));
    std::time_t t = std::time(0);  // get time now
    std::tm* now = std::localtime(&t);
    h = now->tm_hour;
    m = now->tm_min;
    s = now->tm_sec;
    if (h == 0) {
        h = 12;
        AM = "AM";
    } else if (h > 12) {
        h -= 12;
        AM = "PM";
    }
    set_clock_angles();
    glutInit(&argc, argv);                      // Initialize GLUT
    glutInitDisplayMode(GLUT_DOUBLE);           // Enable double buffered mode
    glutInitWindowSize(win_width, win_height);  // Set the window's initial width & height - non-square
    glutInitWindowPosition(100, 200);           // Position the window's initial top-left corner
    glutCreateWindow("ClockTickTock");          // Create window with the given title
    glClearColor(0.0, 0.0, 0.0, 1.0);           // Set background color to black
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, win_width, 0.0, win_height, -1.0, 1.0);

    update_clock_coords();
    glutDisplayFunc(display);  // Register callback handler for window re-paint event
    glutReshapeFunc(reshape);  // Register callback handler for window re-size event
    glutTimerFunc(0, timer, 0);
    glutMainLoop();  // Enter the infinite event-processing loop
    return 0;
}
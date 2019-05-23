//header for your drawing functions

#ifndef GLTEMPLATE_DISPLAY_H_
#define GLTEMPLATE_DISPLAY_H_

#include <GL/freeglut.h>

void initialize();

void reshape(int _w, int _h);
void keyEvent(unsigned char _key, int _x, int _y);
void specKeyEvent(int _key, int _x, int _y);
#ifdef SIMULATION_2D
void mouseClick(int _button, int _state, int _x, int _y);
void mouseDrag(int _x, int _y);
#endif
void mouseMove(int _x, int _y);
void refresh();
void timer(int value);
#endif

#include "stdafx.h"
#include "display.h"
#include <math.h>
#include <ctime>
#include <cstdio>

#ifdef SIMULATION_2D
#include "fluidCube2D.h"
	FluidCube2D *cube;
#else
#include "fluidCube3D.h"
	FluidCube3D *cube;

float ll = 2*LENGTH, seita = 0, fai = 0;
float px = ll * cosf(seita) * cosf(fai);
float py = ll * sinf(seita);
float pz = ll * cosf(seita) * sinf(fai);
#endif

int wide, height;

void initialize(){

#ifdef SIMULATION_2D

	wide = (_W+2) * GRIDSIZE + 20;
	height = (_H+2) * GRIDSIZE + 20;
	cube = new FluidCube2D(VISCOSITY, FRAMERATE, MYSCENE, MYRENDER);

	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE |GLUT_RGBA  | GLUT_STENCIL
                      | GLUT_ACCUM);
	glutInitWindowPosition(10, 10); //initial position of the window on the screen
	glutInitWindowSize(wide, height);
	glutCreateWindow("2D Fluid Simulation"); //create a window and set its focus

	//for the current window, set the callback functions:
	glutDisplayFunc(refresh); //infinite loop to draw on the window
	glutTimerFunc(0, timer, 0);
	//glutReshapeFunc(reshape); //called when the window is resized
	//glutMouseFunc(mouseClick); //called when the mouse is clicked in the window
	//glutMotionFunc(mouseDrag); //called when the mouse is dragged after being clicked

#else 

	cube = new FluidCube3D(VISCOSITY, FRAMERATE, MYSCENE, MYRENDER);

	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE |GLUT_RGBA  | GLUT_STENCIL
                      | GLUT_ACCUM);
	glutInitWindowPosition(10, 10); //initial position of the window on the screen
	glutInitWindowSize(700, 700);
	glutCreateWindow("3D Fluid Simulation"); //create a window and set its focus

	glutDisplayFunc(refresh); //infinite loop to draw on the window
	glutTimerFunc(0, timer, 0);
	glutReshapeFunc(reshape); //called when the window is resized
	glutKeyboardFunc(keyEvent); //called when a standard key is pressed
	glutSpecialFunc(specKeyEvent); //called when a special key is pressed (ie. enter);
	
	srand(time(0));

#endif
	
	//glutIdleFunc(refresh);
	//glutPassiveMotionFunc(mouseMove); //called when the mouse moves (with/without click)
}

void reshape(int _w, int _h){
	PRINT("start - reshape(" <<_w << "," << _h <<")");
	//glutPostRedisplay();
	wide = _w;
	height = _h;

	// Prevent a divide by zero, when window is too short
	// (you cant make a window of zero width).
	if (_h == 0)
		_h = 1;

	float ratio =  _w * 1.0 / _h;

	// Use the Projection Matrix
	glMatrixMode(GL_PROJECTION);

	// Reset Matrix
	glLoadIdentity();

	// Set the viewport to be the entire window
	glViewport(0, 0, _w, _h);

	// Set the correct perspective.
	gluPerspective(45,ratio,1,100);

	// Get Back to the Modelview
	glMatrixMode(GL_MODELVIEW);
	
	
	PRINT("done - reshape(" <<_w << "," << _h <<")");
}

#ifdef SIMULATION_2D

int lastx, lasty;
void mouseClick(int _button, int _state, int _x, int _y){
	PRINT("mouseClick(" <<_button << "," << _state << "," << _x << "," << _y <<")");
	if(_state)
	{
		PRINT("released");
	}
	else
	{
		PRINT("clicked");
		lastx = _x;
		lasty = height -_y;
	}
}
#else

void keyEvent(unsigned char _key, int _x, int _y){
	
	PRINT("keyEvent(" << _key << "," << _x << "," << _y <<")");

	switch (_key)
	{
		case 27:  // ESC
			exit(0);
			break;
		case 'x':
			seita = 0;
			fai = 0;
			ll = 2 * LENGTH;
			px = ll * cosf(seita) * cosf(fai);
			py = ll * sinf(seita);
			pz = ll * cosf(seita) * sinf(fai);
			break;
		case 'y':
			seita = PI / 2;
			fai = 0;
			ll = 2 * LENGTH;
			px = ll * cosf(seita) * cosf(fai);
			py = ll * sinf(seita);
			pz = ll * cosf(seita) * sinf(fai);
			break;
		case 'z':
			seita = 0; 
			fai = PI / 2;
			ll = 2 * LENGTH;
			px = ll * cosf(seita) * cosf(fai);
			py = ll * sinf(seita);
			pz = ll * cosf(seita) * sinf(fai);
			break;
		default:
			break;
	}
}

void specKeyEvent(int _key, int _x, int _y){

	float fraction = 0.01f * LENGTH;

	PRINT("speckeyEvent(" << _key << "," << _x << "," << _y <<")");
	switch (_key)
	{
		case GLUT_KEY_LEFT:
			fai -= 0.01f;
			px = ll * cosf(seita) * cosf(fai);
			py = ll * sinf(seita);
			pz = ll * cosf(seita) * sinf(fai);
			break;
		case GLUT_KEY_RIGHT:
			fai += 0.01f;
			px = ll * cosf(seita) * cosf(fai);
			py = ll * sinf(seita);
			pz = ll * cosf(seita) * sinf(fai);
			break;
		case GLUT_KEY_UP:
			seita += 0.01f;
			px = ll * cosf(seita) * cosf(fai);
			py = ll * sinf(seita);
			pz = ll * cosf(seita) * sinf(fai);
			break;
		case GLUT_KEY_DOWN:
			seita -= 0.01f;
			px = ll * cosf(seita) * cosf(fai);
			py = ll * sinf(seita);
			pz = ll * cosf(seita) * sinf(fai);
			break;
		case GLUT_KEY_PAGE_UP:
			ll -= fraction;
			px = ll * cosf(seita) * cosf(fai);
			py = ll * sinf(seita);
			pz = ll * cosf(seita) * sinf(fai);
			break;
		case GLUT_KEY_PAGE_DOWN:
			ll += fraction;
			px = ll * cosf(seita) * cosf(fai);
			py = ll * sinf(seita);
			pz = ll * cosf(seita) * sinf(fai);
			break;
		default:
			break;
	}
	REPORT(ll);
	REPORT(seita);
	REPORT(fai);
	//cube->render();
}

#endif

void mouseDrag(int _x, int _y){

}

void mouseMove(int _x, int _y){

}
void refresh(){

	PRINT("start - refresh()");
	cube->simulate();
	PRINT("done - refresh()");

}

void timer(int value) {	

	cube->simulate();
	glutTimerFunc(1, timer, 0); // next timer call milliseconds later
}

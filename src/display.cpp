#include "stdafx.h"
#include "display.h"
#include "Configer.h"
#include "Logger.h"
#include "fluidCube2D.h"
#include "fluidCube3D.h"

#include <math.h>

// global variables
bool SIMULATION2D = true;
int LENGTH;
int windowWidth, windowHeight;
float ll, seita = 0, fai = 0;
float px, py, pz;
bool sim_pause = false;
FluidCube *cube;


void initialize(){

	assert(Configer::getConfiger()->getBool("Base", "Simulation2D", SIMULATION2D));

	if(SIMULATION2D)
	{
		int NumGridW, NumGridH, GRIDSIZE;
		assert(Configer::getConfiger()->getInt("Simulation2D", "NumGridW", NumGridW)); 
		assert(Configer::getConfiger()->getInt("Simulation2D", "NumGridH", NumGridH)); 
		assert(Configer::getConfiger()->getInt("Simulation2D", "GridSize", GRIDSIZE)); 		
		windowWidth = (NumGridW + 2) * GRIDSIZE + 20;
		windowHeight = (NumGridH + 2) * GRIDSIZE + 20;
		cube = new FluidCube2D();

		glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE |GLUT_RGBA  | GLUT_STENCIL
	                      | GLUT_ACCUM);
		glutInitWindowPosition(10, 10); //initial position of the window on the screen
		glutInitWindowSize(windowWidth, windowHeight);
		glutCreateWindow("2D Fluid Simulation"); //create a window and set its focus

		//for the current window, set the callback functions:
		glutDisplayFunc(refresh); //infinite loop to draw on the window
		glutTimerFunc(0, timer, 0);
		glutKeyboardFunc(keyEvent); //called when a standard key is pressed
	}
	else 
	{
		int NUMGRIDX, GRIDSIZE;
		assert(Configer::getConfiger()->getInt("Simulation3D", "NumGridX", NUMGRIDX)); 
		assert(Configer::getConfiger()->getInt("Simulation3D", "GridSize", GRIDSIZE)); 		
		assert(Configer::getConfiger()->getInt("Simulation3D", "WindowWidth", windowWidth)); 
		assert(Configer::getConfiger()->getInt("Simulation3D", "WindowHeight", windowHeight)); 
		LENGTH = NUMGRIDX * GRIDSIZE;
		ll = 2 * LENGTH;
		px = ll * cosf(seita) * cosf(fai);
		py = ll * sinf(seita);
		pz = ll * cosf(seita) * sinf(fai);
		cube = new FluidCube3D();

		glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE |GLUT_RGBA  | GLUT_STENCIL
	                      | GLUT_ACCUM);
		glutInitWindowPosition(10, 10); //initial position of the window on the screen
		glutInitWindowSize(windowWidth, windowHeight);
		glutCreateWindow("3D Fluid Simulation"); //create a window and set its focus

		glutDisplayFunc(refresh); //infinite loop to draw on the window
		glutTimerFunc(0, timer, 0);
		glutReshapeFunc(reshape); //called when the window is resized
		glutKeyboardFunc(keyEvent); //called when a standard key is pressed
		glutSpecialFunc(specKeyEvent); //called when a special key is pressed (ie. enter);
	}
}

void reshape(int _w, int _h){
	PRINT("start - reshape(" <<_w << "," << _h <<")");
	//glutPostRedisplay();
	windowWidth = _w;
	windowHeight = _h;

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
		lasty = windowHeight -_y;
	}
}

void keyEvent(unsigned char _key, int _x, int _y){
	
	PRINT("keyEvent(" << _key << "," << _x << "," << _y <<")");

	switch (_key)
	{
		case 27:  // ESC
			exit(0);
			break;
		case 'e':
			sim_pause = !sim_pause;
			break;
		default:
			break;
	}

	if(!SIMULATION2D)
	{
		switch (_key)
		{
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
		}
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

	if(!sim_pause)
		cube->simulate();
	else
		cube->render();

	glutTimerFunc(0, timer, 0); // next timer call milliseconds later
}

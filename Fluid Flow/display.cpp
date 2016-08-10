//implement your drawing functions here
#include "stdafx.h"
#include "fluidCube.h"
#include "display.h"
#include <math.h>
#include <ctime>
#include <cstdio>

#ifdef SIMULATION_2D
FluidCube2D *cube;
#else
FluidCube3D *cube;
float ll = 2 * LENGTH;
float seita = 0, fai = 0;
float px = 2 * LENGTH, py = 0, pz = 0;
float dx = cosf(seita) * cosf(fai);
float dy = sinf(seita);
float dz = cosf(seita) * sinf(fai);
#endif

int count = 0;
int wide, height;

void initialize(){

#ifdef SIMULATION_2D

	wide = (_W+2) * GRIDSIZE + 20;
	height = (_H+2) * GRIDSIZE + 20;
	cube = new FluidCube2D(VISCOSITY, TIMESTEP, CONTAINER);

	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE |GLUT_RGBA  | GLUT_STENCIL
                      | GLUT_ACCUM);
	glutInitWindowPosition(100, 100); //initial position of the window on the screen
	glutInitWindowSize(wide, height);
	glutCreateWindow("2D Fluid Simulation"); //create a window and set its focus

	//for the current window, set the callback functions:
	glutDisplayFunc(refresh); //infinite loop to draw on the window
	//glutTimerFunc(0, timer, 0);
	//glutReshapeFunc(reshape); //called when the window is resized
	//glutMouseFunc(mouseClick); //called when the mouse is clicked in the window
	//glutMotionFunc(mouseDrag); //called when the mouse is dragged after being clicked

#else 

	cube = new FluidCube3D(DIFFUSION, VISCOSITY, TIMESTEP);

	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE |GLUT_RGBA  | GLUT_STENCIL
                      | GLUT_ACCUM);
	glutInitWindowPosition(100, 100); //initial position of the window on the screen
	glutInitWindowSize(800, 800);
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

int lastx,lasty;
clock_t start;
void mouseClick(int _button, int _state, int _x, int _y){
	PRINT("mouseClick(" <<_button << "," << _state << "," << _x << "," << _y <<")");
	if(_state)
	{
		PRINT("released");
		memset(cube->Vx0, 0, sizeof(float) * cube->size);
		memset(cube->Vy0, 0, sizeof(float) * cube->size);

		//there may be some errors with the position
		int xx = _x / GRIDSIZE + 1;
		int yy = (height - _y) / GRIDSIZE ;
		clock_t time = clock() - start;
		float dxx = _x - lastx;
		float dyy = (height -_y) - lasty;
		float dxy = sqrtf(dxx*dxx + dyy*dyy);
		if(dxy == 0)
			return;

		cube->Vx0[IX(xx, yy)] = dxx / dxy * time * DRAGSCALE;
		cube->Vy0[IX(xx, yy)] = dyy / dxy * time * DRAGSCALE;

		cube->simulate();
	}
	else
	{
		PRINT("clicked");
		lastx = _x;
		lasty = height -_y;
		start = clock();
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
			dx = cosf(seita) * cosf(fai);
			dy = sinf(seita);
			dz = cosf(seita) * sinf(fai);
			break;
		case 'y':
			seita = PI / 2;
			fai = 0;
			ll = 2 * LENGTH;
			px = ll * cosf(seita) * cosf(fai);
			py = ll * sinf(seita);
			pz = ll * cosf(seita) * sinf(fai);
			dx = cosf(seita) * cosf(fai);
			dy = sinf(seita);
			dz = cosf(seita) * sinf(fai);
			break;
		case 'z':
			seita = 0; 
			fai = PI / 2;
			ll = 2 * LENGTH;
			px = ll * cosf(seita) * cosf(fai);
			py = ll * sinf(seita);
			pz = ll * cosf(seita) * sinf(fai);
			dx = cosf(seita) * cosf(fai);
			dy = sinf(seita);
			dz = cosf(seita) * sinf(fai);
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
			dx = cosf(seita) * cosf(fai);
			dy = sinf(seita);
			dz = cosf(seita) * sinf(fai);
			break;
		case GLUT_KEY_RIGHT:
			fai += 0.01f;
			px = ll * cosf(seita) * cosf(fai);
			py = ll * sinf(seita);
			pz = ll * cosf(seita) * sinf(fai);
			dx = cosf(seita) * cosf(fai);
			dy = sinf(seita);
			dz = cosf(seita) * sinf(fai);
			break;
		case GLUT_KEY_UP:
			seita += 0.01f;
			px = ll * cosf(seita) * cosf(fai);
			py = ll * sinf(seita);
			pz = ll * cosf(seita) * sinf(fai);
			dx = cosf(seita) * cosf(fai);
			dy = sinf(seita);
			dz = cosf(seita) * sinf(fai);
			break;
		case GLUT_KEY_DOWN:
			seita -= 0.01f;
			px = ll * cosf(seita) * cosf(fai);
			py = ll * sinf(seita);
			pz = ll * cosf(seita) * sinf(fai);
			dx = cosf(seita) * cosf(fai);
			dy = sinf(seita);
			dz = cosf(seita) * sinf(fai);
			break;
		case GLUT_KEY_PAGE_UP:
			px -= dx * fraction;
			py -= dy * fraction;
			pz -= dz * fraction;
			ll -= fraction;
			break;
		case GLUT_KEY_PAGE_DOWN:
			px += dx * fraction;
			py += dy * fraction;
			pz += dz * fraction;
			ll += fraction;
			break;
		default:
			break;
	}
	REPORT(px);
	REPORT(py);
	REPORT(pz);
	cube->draw_dens();
}

#endif

void mouseDrag(int _x, int _y){
	PRINT("mouseDrag(" << _x << "," << _y <<"): displacement from click: (" << _x-lastx << "," << _y-lasty << ")");
}

void mouseMove(int _x, int _y){

}
void refresh(){

	PRINT("start - refresh()");
	cube->simulate();
	PRINT("done - refresh()");

}

void timer(int value) {
/*
#ifdef SIMULATION_2D

	if(count %FLOWTIME == 0)
	{
		memset(cube->Vx0, 0, sizeof(float) * cube->size);
		memset(cube->Vy0, 0, sizeof(float) * cube->size);
		//for(int x = 1; x <= 5; x++)
		for(int y = 1; y <= _H; y++)
		{
			cube->Vx0[IX(1, y)] = SPEED;  //10000~50000 for 2 vertexes
			//cube->Vy0[IX(1, y)] = 0;
		}
		cube->simulate(false);
	}
	else
	{
		glutPostRedisplay();
	}
	if(count < 2025)
		count ++;
#endif
*/	

	glutPostRedisplay();
	glutTimerFunc(FRAMERATE, timer, 0); // next timer call milliseconds later
}

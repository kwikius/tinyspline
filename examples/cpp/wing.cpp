#include "tinyspline.h"
#include <quan/three_d/vect.hpp>
#include <vector>

#ifdef TINYSPLINE_DOUBLE_PRECISION
#error OpenGL requires tinyspline library to be complied for floats
#endif

#if defined(__APPLE__)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h> 
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include <stdio.h>
#include <debugging.h>
#include <assert.h>
#include <math.h>

/*
   show a spar e.g a spline and then another spline at a given distance normal to curve.
   In fact with a spline of 5 points, the curve isnt a great fit except hopefully at interpolated points
*/

extern "C" tsError ts_bspline_clean(tsBSpline * spline);

/********************************************************
*                                                       *
* Modify these lines for experimenting.                 *
*                                                       *
********************************************************/

static tsBSpline spline_orig;
static tsBSpline spline_copy;
static tsBSpline deriv;
static tsBSpline spline_spar;


static GLUnurbsObj *theNurb;

static void perp(const tsReal* x_in, const tsReal * y_in, tsReal* resultx, tsReal* resulty)
{
     // copy as may be used in place
     const tsReal x = *x_in;
     const tsReal y = *y_in;
     *resultx = -y;
     *resulty = x;
}

static tsReal norm (const tsReal * x,const tsReal * y)
{
   return *x * *x + *y * *y;
}

static tsReal magnitude(const tsReal * x, const  tsReal * y)
{
   return sqrt(norm(x,y));
}

static void multiply_vector(const tsReal* x_in, const tsReal * y_in, tsReal val, tsReal* resultx, tsReal* resulty)
{
     *resultx = *x_in * val ;
     *resulty = *y_in * val;
}

static void unit_vector(const tsReal* x, const tsReal * y, tsReal* resultx, tsReal* resulty)
{
   const tsReal mag = magnitude(x,y);
   //check size
   assert ( mag > 1.e-6);
   *resultx = *x / mag;
   *resulty = *y / mag;
}

tsReal* create_spline(std::vector<quan::three_d::vect<double> > const & points_in)
{
   std:;size_t const size = points_in.size();
   std::size_t const length = 3 * size;
   tsReal* arr = new tsReal [length];
   
   for (size_t i = 0U; i < size; ++i){
      quan::three_d::vect<double> const &  vect = points_in[i];
      arr[3 * i] = vect.x;
      arr[3 * i + 1U] = vect.y;
      arr[3 * i +2] = vect.z;
   }
   return arr;
}


void setup()
{
   printf ("spar test\n\n");

   constexpr size_t dimension = 3U;
   
   std::vector<quan::three_d::vect<double> > vect =
   {
      {1,-1,0}
      ,{-1,2,0}
      ,{1,4,0}
      ,{4,3,0}
      ,{7,5,0}
      ,{8,6,0}
   };

   const size_t num_points = vect.size(); 

   tsReal * points = create_spline(vect);

//	tsReal points[num_points * dimension];
//	points[0] = 1;
//	points[1] = -1;
//	points[2] = 0;

//	points[3] = -1;
//	points[4] = 2;
//	points[5] = 0;

//	points[6] = 1;
//	points[7] = 4;
//	points[8] = 0;

//	points[9] = 4;
//	points[10] = 3;
//	points[11] = 0;

//	points[12] = 7;
//	points[13] = 5;
//	points[14] = 0;

	assert(ts_bspline_interpolate_cubic(points, num_points, dimension, &spline_orig)== TS_SUCCESS);
  // ts_bspline_print(&spline_orig);
   assert(ts_bspline_copy(&spline_orig,&spline_copy)== TS_SUCCESS);
   assert(ts_bspline_clean(&spline_copy) == TS_SUCCESS);
   assert(ts_bspline_derive(&spline_copy,&deriv) == TS_SUCCESS);

   // create the spar 
  // const size_t num_segments = sizeof(points) / (3 * sizeof(points[0])) - 1U;

   const size_t num_segments = num_points -1U;
   
   for( size_t i = 0U; i <= num_segments; ++i){
      tsReal u = (tsReal)i / num_segments; // 0 < 0 ; num_segments / num_segments -> 1
      if ( u < 0.f){
         u = 0.f;
      }else {
         if ( u > 1.f){
            u = 1.f;
         }
      }
      tsDeBoorNet spline_net;
      tsDeBoorNet deriv_net;
      assert(ts_bspline_evaluate(&spline_orig, u, &spline_net)== TS_SUCCESS);

      tsReal px = spline_net.result[0];
      tsReal py = spline_net.result[1];
      tsReal pz = spline_net.result[2];

      assert(ts_bspline_evaluate(&deriv, u, &deriv_net)== TS_SUCCESS);

      tsReal dx = deriv_net.result[0];
      tsReal dy = deriv_net.result[1];
      tsReal dz = spline_net.result[2];
      
      perp(&dx,&dy,&dx,&dy);
      unit_vector(&dx,&dy,&dx,&dy);
      const tsReal spar_width = 0.05;
      multiply_vector(&dx,&dy,spar_width,&dx,&dy);

      points[3 * i] = px + dx;
      points[3 * i + 1U] = py + dy;
      points[3 * i + 2U] = pz + dz;
     
      ts_deboornet_free(&spline_net);
      ts_deboornet_free(&deriv_net);
   }
   assert(ts_bspline_interpolate_cubic(points, 5, 3, &spline_spar) == TS_SUCCESS);
   
}

void tear_down()
{
	ts_bspline_free(&spline_copy);
   ts_bspline_free(&spline_orig);
   ts_bspline_free(&deriv);
}

// function to show the two splines in varying colours, where each colour is a array of dim 3
static void draw_spar(GLfloat* spline_colour,GLfloat* ctrlp_colour)
{

	glLineWidth(3);
   glColor3f(spline_colour[0],spline_colour[1],spline_colour[2]);
	gluBeginCurve(theNurb);
		gluNurbsCurve(
			theNurb, 
			(GLint)spline_copy.n_knots,
			spline_copy.knots, 
			(GLint)spline_copy.dim,
			spline_copy.ctrlp, 
			(GLint)spline_copy.order,
			GL_MAP1_VERTEX_3
		);
	gluEndCurve(theNurb);

   glColor3f(ctrlp_colour[0],ctrlp_colour[1],ctrlp_colour[2]);
   gluBeginCurve(theNurb);
		gluNurbsCurve(
			theNurb, 
			(GLint)spline_spar.n_knots,
			spline_spar.knots, 
			(GLint)spline_spar.dim,
			spline_spar.ctrlp, 
			(GLint)spline_spar.order,
			GL_MAP1_VERTEX_3
		);
	gluEndCurve(theNurb);
}

/*
   draw the tangent of the spline at u in colour
   derivative has been evaluated by calling ts_bspline_derive on the spline
*/

static GLfloat red[]= { 1.0,0.0,0.0};
static GLfloat green[]= { 0,1.0,0.0};
static GLfloat blue[]= { 0.0,0.0,1.0};
static GLfloat grey[]= { 0.5,0.5,0.5};
static GLfloat yellow[]= { 1.0,1.0,0.0};
static GLfloat orange[]= { 1.0,0.75,0.0};

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	

   draw_spar(red,green);

	glutSwapBuffers();
	glutPostRedisplay();
}

/********************************************************
*                                                       *
* Framework                                             *
*                                                       *
********************************************************/
void nurbsError(GLenum errorCode)
{
   const GLubyte *estring;

   estring = gluErrorString(errorCode);
   fprintf (stderr, "Nurbs Error: %s\n", estring);
   exit (0);
}
   
void init(void)
{
	glClearColor (0.0, 0.0, 0.0, 0.0);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	theNurb = gluNewNurbsRenderer();
	gluNurbsProperty (theNurb, GLU_SAMPLING_TOLERANCE, 10.0);
	gluNurbsCallback(theNurb, GLU_ERROR, (GLvoid (*)()) nurbsError);
	setup();
}

void reshape(int w, int h)
{
   glViewport(0, 0, (GLsizei) w, (GLsizei) h);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   gluPerspective (80.0, (GLdouble)w/(GLdouble)h, 3.0, 8.0);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   glTranslatef (-2.5, -1.0, -7.0);
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize (500, 500);
	glutInitWindowPosition (100, 100);
	glutCreateWindow(argv[0]);
	init();
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutMainLoop();
	tear_down();
	return 0; 
}

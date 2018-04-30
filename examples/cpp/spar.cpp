#include "tinyspline.h"

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
   In fact with a spline of 5 points, the curve isnt a great fit.
   Howvere increase the evaluated points greatly and a good fit results
*/

extern "C" tsError ts_bspline_clean(tsBSpline * spline);

/********************************************************
*                                                       *
* Modify these lines for experimenting.                 *
*                                                       *
********************************************************/

// wrap raw_points

#include "andy/points_wrapper.hpp"

namespace {

   tsBSpline spline_orig;
   tsBSpline spline_copy;
   tsBSpline deriv;
   tsBSpline spline_spar;

   GLUnurbsObj *theNurb;

   quan::three_d::vect<tsReal> get_result(tsDeBoorNet const & net)
   {
      assert(net.dim == 3U);
      return { net.result[0],net.result[1],net.result[2]};
   }

   void interpolate(ts_points_wrapper const & points,tsBSpline & spline)
   {
      assert(ts_bspline_interpolate_cubic(points.get_points(), points.get_num_points() / 3, 3, &spline) == TS_SUCCESS);
   }

   void interpolate(std::vector<quan::three_d::vect<tsReal> > const & vect,tsBSpline & spline)
   {
      ts_points_wrapper points(vect);
      interpolate(points,spline);
   }

   /*
   evaluate the spline at u where u is between 0 and 1
*/
   quan::three_d::vect<tsReal> evaluate(tsBSpline const & spline,tsReal const & u)
   {
      tsDeBoorNet spline_net;
      assert(ts_bspline_evaluate(&spline, u, &spline_net)== TS_SUCCESS);
      auto const p = get_result(spline_net);
      ts_deboornet_free(&spline_net);
      return p;
   }

} // namespace

void setup()
{
   printf ("spar test\n\n");

   std::vector<quan::three_d::vect<tsReal> > vect =
   {
      {1,-1,0}
      ,{-1,2,0}
      ,{1,4,0}
      ,{4,3,0}
      ,{7,5,0}
   };

   const size_t num_segments = vect.size() - 1U;
   constexpr tsReal spar_thickness = 0.8;

   interpolate(vect,spline_orig);

   assert(ts_bspline_copy(&spline_orig,&spline_copy) == TS_SUCCESS);
   assert(ts_bspline_clean(&spline_copy) == TS_SUCCESS);
   assert(ts_bspline_derive(&spline_copy,&deriv) == TS_SUCCESS);
 
   vect.clear();
   constexpr size_t split = 100;
   for( size_t i = 0U; i <= num_segments * split; ++i){

      tsReal const u = quan::constrain<tsReal>((tsReal)i / (num_segments * split),0.f,1.f); 

      auto const p = evaluate(spline_orig,u);
      auto const p1 = evaluate(deriv,u);

      auto const d = unit_vector(perp_vector_axis_z(p1)) * spar_thickness;
      
      vect.push_back(p + d);
   }
  // points = vect;
   interpolate(vect,spline_spar);
   
}

void tear_down()
{
	ts_bspline_free(&spline_copy);
   ts_bspline_free(&spline_orig);
   ts_bspline_free(&deriv);
   ts_bspline_free(&spline_spar);
}

static void do_glNurbsCurve(tsBSpline const &  spline, GLfloat * spline_colour, int line_width = 3)
{
   glLineWidth(line_width);
   glColor3f(spline_colour[0],spline_colour[1],spline_colour[2]);
	gluBeginCurve(theNurb);
		gluNurbsCurve(
			theNurb, 
			(GLint)spline.n_knots,
			spline.knots, 
			(GLint)spline.dim,
			spline.ctrlp, 
			(GLint)spline.order,
			GL_MAP1_VERTEX_3
		);
	gluEndCurve(theNurb);
}

// function to show the two splines in varying colours, where each colour is a array of dim 3
static void draw_spar(GLfloat* spline_colour,GLfloat* ctrlp_colour)
{
   do_glNurbsCurve(spline_copy,spline_colour);
   do_glNurbsCurve(spline_spar,ctrlp_colour);
}

static GLfloat red[]= { 1.0,0.0,0.0};
static GLfloat green[]= { 0,1.0,0.0};

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

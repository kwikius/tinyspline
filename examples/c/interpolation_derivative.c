#include "tinyspline.h"

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

tsBSpline spline;
tsBSpline deriv;
GLUnurbsObj *theNurb;

/*

A simple solution is to reduce to multiplicity of an internal knot by 1 
if it is equals to order (deg + 1) and if the two corresponding control points 
(https://github.com/msteinbeck/tinyspline/blob/master/library/tinyspline.h#L240) 
are equals
 (you can use: https://github.com/msteinbeck/tinyspline/blob/master/library/tinyspline.h#L889 
to check if two control points are "equals"). 
When reducing the multiplicity of such a knot we have to merge the two control points by 
shifting the control point array to left so that the second control point corresponding 
to our knot is replaced by the next one and so forth:

// control points 
[p0, p1, p2, *p3, p4, p5, p6] // p2 and p3 need to be merged
[p0, p1, p2, p4, p5, p6] // fixed

Similarly, we can shift the knot array:

[0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4, 5, 6, 7, 8, 9, 9, 9, 9] // 4 needs to be merged
[0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 5, 6, 7, 8, 9, 9, 9, 9] // fixed

This step can be done before performing the actual derivation.

*/

// find multiplicity of a knot
static tsError get_knot_multiplicity(tsBSpline* spline,size_t knot_idx, size_t * result)
{
     assert(spline != NULL);
     assert(result != NULL);
     assert(knot_idx < spline->n_knots);

     tsDeBoorNet net;
     const tsError e = ts_bspline_evaluate(spline,spline->knots[knot_idx],&net);
     if ( e == TS_SUCCESS){
        *result = net.s;
     }
     return e;
}

/*
   get the 2 ? control points corresponding to a knot at index knot_idx in the spline
*/
static tsError get_corresponding_control_points(
      tsBSpline* spline,size_t knot_idx,
      tsReal const * result_ctrlpa, tsReal const * result_ctrlpb)
{
    return TS_DIM_ZERO;// error for now
}

/********************************************************
*                                                       *
* Modify these lines for experimenting.                 *
*                                                       *
********************************************************/
void setup()
{
	tsReal points[15];
	points[0] = 1;
	points[1] = -1;
	points[2] = 0;
	points[3] = -1;
	points[4] = 2;
	points[5] = 0;
	points[6] = 1;
	points[7] = 4;
	points[8] = 0;
	points[9] = 4;
	points[10] = 3;
	points[11] = 0;
	points[12] = 7;
	points[13] = 5;
	points[14] = 0;

	ts_bspline_interpolate_cubic(points, 5, 3, &spline);

   const size_t n_knots = spline.n_knots;
   const size_t deg = spline.deg;
   const size_t dim = spline.dim;
   for ( size_t knot_idx = 0U; knot_idx < n_knots; ++knot_idx){
      size_t multiplicity = 0U;
      tsError e = get_knot_multiplicity(&spline,knot_idx,&multiplicity);
      assert ( e == TS_SUCCESS);
      if ( multiplicity == (deg + 1)){
         tsReal const * ctrl_pointsa;
         tsReal const * ctrl_pointsb;
         tsError e = get_corresponding_control_points(&spline,knot_idx,ctrl_pointsa,ctrl_pointsb);
         assert ( e == TS_SUCCESS);
         // compare the points
         const tsReal min_useful_distance = (tsReal)1.e-6;// ? for now
         if (ts_ctrlp_dist2(ctrl_pointsa,ctrl_pointsb, dim) < min_useful_distance){
            printf("problem knot at %u, value %f:\n",(unsigned int)knot_idx,spline.knots[knot_idx]) ;
         }
      }
   }

//########################################
   //N.B Will assert here until this is sorted
   assert (ts_bspline_derive(&spline,&deriv) == TS_SUCCESS);
//##########################################
	ts_bspline_print(&spline);
   ts_bspline_print(&deriv);
}

void tear_down()
{
	ts_bspline_free(&spline);
}

void display(void)
{
	tsDeBoorNet net;
	size_t i;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	/* draw spline */
	glColor3f(1.0, 1.0, 1.0);
	glLineWidth(3);
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

	/* draw control points */
	glColor3f(1.0, 0.0, 0.0);
	glPointSize(5.0);
	glBegin(GL_POINTS);
	  for (i = 0; i < spline.n_ctrlp; i++) 
		 glVertex3fv(&spline.ctrlp[i * spline.dim]);
	glEnd();
	
	/* draw evaluation */
	glColor3f(0.0, 0.0, 1.0);
	glPointSize(5.0);
	ts_bspline_evaluate(&spline, 0.5f, &net);
	glBegin(GL_POINTS);
		glVertex3fv(net.result);
	glEnd();
	ts_deboornet_free(&net);
	
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

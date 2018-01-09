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

[0, 0, 0, 0, 1, 2, 3, 4, 4, 4, *4, 5, 6, 7, 8, 9, 9, 9, 9] // 4 needs to be merged
[0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 5, 6, 7, 8, 9, 9, 9, 9] // fixed

This step can be done before performing the actual derivation.

static tsError remove_knot(tsBSpline* spline_in,size_t knot_index,tsBSpline* spline_out )
{
   make a knot vector of size one less than spline and copy all < knot_index and all > knot_index to it
   assert( spline_in != spline_out);

   

}

*/

/********************************************************
*                                                       *
* Modify these lines for experimenting.                 *
*                                                       *
********************************************************/
void setup()
{

   printf ("interpolated spline derivative test with control point cleaning\n\n");
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

   printf("original spline:\n");
   ts_bspline_print(&spline);

   printf("\nchecking spline for derivative issues\n");

  // const size_t n_knots = spline.n_knots;
   const size_t order = spline.order;
   const size_t dim = spline.dim;
   const tsReal min_useful_distance = (tsReal)1.e-6;// ? for now

   size_t knot_idx = 0U;
   
   // n.b that number of  knots and control points will change dynamically in the loop
   while(knot_idx < spline.n_knots){

      // count the multiplicity of the knot
      size_t multiplicity = 1U;
      const tsReal cur_knot_value = spline.knots[knot_idx];
      while ( (knot_idx < (spline.n_knots-1)) && ( spline.knots[knot_idx +1] == cur_knot_value) ){
         ++knot_idx;
         ++multiplicity;
      }
      if ( multiplicity == order){
         const size_t idxa = (knot_idx - multiplicity) * dim;
         const size_t idxb = idxa + dim;
         // get and compare the control points
         const tsReal* ctrl_pointsa = &spline.ctrlp[idxa];
         const tsReal* ctrl_pointsb = &spline.ctrlp[idxb];
     
         if (ts_ctrlp_dist2(ctrl_pointsa,ctrl_pointsb, dim) < min_useful_distance){

            printf("found problem knot at index %u, removing...\n\n",(unsigned int)knot_idx) ;

            // create and init a temporary spline for the mods
            tsBSpline mod_spline;
            {
               mod_spline.deg = spline.deg;
               mod_spline.order = spline.order;
               mod_spline.dim = spline.dim;
               mod_spline.n_ctrlp = spline.n_ctrlp-1U;
               mod_spline.n_knots = spline.n_knots-1U;
               const size_t n_raw_bytes_for_real_array = sizeof(tsReal) *  (mod_spline.n_knots + mod_spline.n_ctrlp * mod_spline.dim);
               mod_spline.ctrlp = (tsReal*) malloc(n_raw_bytes_for_real_array);
               // TODO check malloc result
               mod_spline.knots = mod_spline.ctrlp + mod_spline.n_ctrlp * mod_spline.dim;
            }
            // copy control points prior to and including first control point coresponding to current knot
            const tsReal* ctrlp_src = spline.ctrlp;
            tsReal* ctrlp_dest = mod_spline.ctrlp;

            size_t reals_copied = 0U;

            while ( ctrlp_src != ctrl_pointsb){
               *ctrlp_dest++ = *ctrlp_src++;
                ++reals_copied;
            }
            // skip copying 2nd control point coresponding to current knot
            ctrlp_src +=  mod_spline.dim;
            // copy the rest of the control points over
            while ( reals_copied < (mod_spline.n_ctrlp * mod_spline.dim) ){
                *ctrlp_dest++ = *ctrlp_src++;
                ++reals_copied;
            }
            // copy all knots before current
            const tsReal* knots_src = spline.knots;
            tsReal* knots_dest = mod_spline.knots;
            size_t idx = 0U;
            while ( idx < knot_idx){
               knots_dest[idx] = knots_src[idx];
               ++idx;
            }
            // copy all knots after current
            while (idx < mod_spline.n_knots){
               knots_dest[idx] = knots_src[idx + 1];
               ++idx;
            } 
            // clear the original spline and replace the contents with the new one  
            ts_bspline_free(&spline);
            ts_bspline_move(&mod_spline,&spline);
            
            printf("... spline after problem knot removed:\n");
            ts_bspline_print(&spline);
         }
      }else{
         // nothing to do
         ++knot_idx;
      }
   }

   printf("spline cleaned:\n\n");

   assert (ts_bspline_derive(&spline,&deriv) == TS_SUCCESS);

   printf("show spline derivative:\n");
	//ts_bspline_print(&spline);
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

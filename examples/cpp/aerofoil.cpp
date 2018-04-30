#include "tinyspline.h"

#ifndef TINYSPLINE_FLOAT_PRECISION
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
#include <fstream>
#include <sstream>


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
#include "aerofoil.hpp"

namespace {

   tsBSpline spline_orig;
   tsBSpline spline_copy;
   tsBSpline deriv;
   tsBSpline spline_spar;

   std::string const aerofoil_filename = "/home/andy/cpp/projects/aerofoil/Sections/selig/mh32.dat";

   GLUnurbsObj *theNurb;
}

void setup()
{
   printf ("aerofoil test\n\n");

   aerofoil foil;
   if ( !foil.load(aerofoil_filename,std::cout)){
      printf("aerofoil load failed\n");
      return;
   }

   // vector to get the coords from the foil
   std::vector<quan::three_d::vect<tsReal> > vect;

   double constexpr chord =  10;

   // create the section from the aerofoil and chord and put in vect
   for ( std::size_t i = 0U; i < foil.get_num_coords();++i){
        auto p = foil.get_coord(i) * chord;
        vect.push_back({p.x,p.y,0.0});
   }

   const size_t num_segments = vect.size() - 1U;
   constexpr tsReal spar_thickness = 0.1;

   // create the spline from the vect
   interpolate(vect,spline_orig);

   assert(ts_bspline_copy(&spline_orig,&spline_copy) == TS_SUCCESS);
   // cleabn the spline to do the derivative
   assert(ts_bspline_clean(&spline_copy) == TS_SUCCESS);
   // create the derivative spline
   assert(ts_bspline_derive(&spline_copy,&deriv) == TS_SUCCESS);
 
   // create a vector of points at the skin thickness from the original section
   vect.clear();
   // If there arent many points on the original curve, then evaluate 
   constexpr size_t split = 1;
   for( size_t i = 0U; i <= num_segments * split; ++i){

      tsReal const u = quan::constrain<tsReal>((tsReal)i / (num_segments * split),0.f,1.f); 

      auto const p = evaluate(spline_orig,u);
      auto const p1 = evaluate(deriv,u);

      auto const d = unit_vector(perp_vector_axis_z(p1)) * spar_thickness;
      
      vect.push_back(p + d);
   }
   // create the spline representing the position of the skin
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

namespace {

void strip_leading(std::string & str)
{
   while (isspace(str[0])){
     str = str.substr(1,std::string::npos);
   }
}

void strip_trailing(std::string & str)
{
   for(;;){
      auto len = str.length();
      if (isspace(str[len-1])){
        str = str.substr(0,len-1);
      }else{
        break;
      }
   }
}
}

bool aerofoil::load(
    std::string const & data_file,
    std::ostream & e)
{

   std::ifstream f(data_file);
   if ( ! f || f.fail()){
     e << "Invalid input file";
     return false;
   }

   std::string aerofoil_name;
   for(;;){
     if ( !f || f.eof() ){
      e << "expected aerofoil name and coords";
      return false;
     }
     getline(f,aerofoil_name);
     strip_leading(aerofoil_name);
     strip_trailing(aerofoil_name);
     if ( aerofoil_name != "") {
          break;
     }
     else continue;
   }
   std::vector<quan::two_d::vect<double> > coords;
   
   for(;;){
    if ( !f || f.eof() ) break;
     std::string line;
     getline(f,line);
     strip_leading(line);
     strip_trailing(line);
     if( line !=""){
     std::istringstream is(line);
     quan::two_d::vect<double> coord;
     
     is >> coord.x;
     // fix some selig files seem to have commas etc here
     is >> coord.y;
     coords.push_back(coord);
     }
   }
   if (coords.size() < 4){
    e << "insufficient coords in aerofoilfile";
    return false;
   }
   this->m_name = aerofoil_name;
   this->m_coords.clear();
   this->m_coords = coords;
   return true;

}


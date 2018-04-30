#ifndef TINY_SPLINE_LIBRARY_ANDY_POINTS_WRAPPER_HPP_INCLUDED
#define TINY_SPLINE_LIBRARY_ANDY_POINTS_WRAPPER_HPP_INCLUDED


#include <cassert>
#include <vector>
#include <quan/three_d/vect.hpp>
#include <quan/two_d/vect.hpp>
#include <quan/constrain.hpp>

#include "tinyspline.h"

struct ts_points_wrapper;

quan::three_d::vect<tsReal> get_result(tsDeBoorNet const & net);

// create the spline as a result of interpolating the points in points
void interpolate(ts_points_wrapper const & points,tsBSpline & spline);

// create the spline as a result of interpolating the points in vect
void interpolate(std::vector<quan::three_d::vect<tsReal> > const & vect,tsBSpline & spline);

// evaluate the spline at "knot" u
quan::three_d::vect<tsReal> evaluate(tsBSpline const & spline,tsReal const & u);

struct ts_points_wrapper{

   ts_points_wrapper(size_t num_points);
   ts_points_wrapper(ts_points_wrapper const & in);
   ts_points_wrapper(std::vector<quan::three_d::vect<tsReal> > const & vect);

   ts_points_wrapper(ts_points_wrapper && temp);

   ts_points_wrapper & operator = (ts_points_wrapper const & temp);
   ts_points_wrapper & operator = (std::vector<quan::three_d::vect<tsReal> > const & vect);
   ts_points_wrapper & operator = (ts_points_wrapper && temp);

   ~ts_points_wrapper();

   operator tsReal * () { assert(m_points); return m_points;}

   tsReal* get_points()const { return m_points;}
   size_t get_num_points() const {return m_num_points;}
   tsReal  & operator[] (size_t i)
   {
     return m_points[i];
   }

   tsReal const & operator[] (size_t i)const
   {
      return m_points[i];
   }

   private:
   void null_points() { m_num_points = 0; m_points = nullptr; }
   tsReal * move_points()
   {
      tsReal * temp = m_points; m_points = nullptr; m_num_points = 0U; return temp;
   }

   tsReal * m_points;
   size_t m_num_points;
};

#endif // TINY_SPLINE_LIBRARY_ANDY_SPLINE_WRAPPER_HPP_INCLUDED

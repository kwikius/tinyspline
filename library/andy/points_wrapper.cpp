
#include "points_wrapper.hpp"

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

namespace {
void set_point(tsReal * pt, quan::three_d::vect<tsReal> const & vect)
{
    pt[0] = vect.x;
    pt[1] = vect.y;
    pt[2] = vect.z;
}

tsReal* convert_to_ts_bspline_points(std::vector<quan::three_d::vect<tsReal> > const & points_in)
{
   std::size_t const size = points_in.size();
   tsReal* arr = new tsReal [3 * size];
   
   for (size_t i = 0U; i < size; ++i){
      set_point(&arr[3 *i],points_in[i]);
   }
   return arr;
}
}

ts_points_wrapper::ts_points_wrapper(size_t num_points) 
:m_points{ new tsReal [num_points] } ,m_num_points{num_points}{}

ts_points_wrapper::ts_points_wrapper(ts_points_wrapper const & in)
: m_points{ new tsReal [in.m_num_points] } ,m_num_points{in.m_num_points}
{
   for ( size_t i = 0U; i < m_num_points; ++i){
      m_points[i] = in.m_points[i];
   }
}

ts_points_wrapper::ts_points_wrapper(std::vector<quan::three_d::vect<tsReal> > const & vect)
: m_points{convert_to_ts_bspline_points(vect)},m_num_points{vect.size() * 3}{}

ts_points_wrapper::ts_points_wrapper(ts_points_wrapper && temp): m_points{temp.get_points()}, m_num_points{temp.get_num_points()}
{
   temp.null_points();
}

ts_points_wrapper::~ts_points_wrapper()
{
   if (m_points != nullptr){ delete [] m_points;}
}

ts_points_wrapper & ts_points_wrapper::operator = (ts_points_wrapper const & temp)
{
    if ( m_num_points != temp.m_num_points){
      delete [] m_points;
      m_points = new tsReal [temp.m_num_points];
      m_num_points = temp.m_num_points;
    }
    for (size_t i = 0U; i < m_num_points; ++i){
         m_points[i] = temp.m_points[i];
    }
    return *this;
}

ts_points_wrapper & ts_points_wrapper::operator = (std::vector<quan::three_d::vect<tsReal> > const & vect)
{
    
    auto const vect_num_points = vect.size() * 3;
    if ( m_num_points != vect_num_points){
      delete [] m_points;
      m_points = new tsReal [vect_num_points];
      m_num_points = vect_num_points;
    }
    for (size_t i = 0U; i < m_num_points/3; ++i){
       //m_points[i] = temp.m_points[i];
      set_point(&m_points[3 *i],vect[i]);
    }
    return *this;
}

ts_points_wrapper & ts_points_wrapper::operator = (ts_points_wrapper && temp)
{
    if (m_points != nullptr){ delete [] m_points;}
    m_num_points = temp.m_num_points;
    m_points = temp.move_points();
    return *this;
}



#include "tinyspline.h"


/*
   sams interface as the quan::two_d::cubic_spline?

   This differs since for a given x there is no unique y = f(x)
   firstly the point on the curve is a position between 0 and 1
   second a curve can cross x at multiple positions

   Therefore no, dont make same interface as quan::two_d::cubic_spline

   possibly we can verify that the interpolation points on the spline increase in x

   we will need to work to find the y value is a function of  a particular x point.
*/

namespace quan{ namespace two_d{

   // use a double tinyspline
   // we have to copy values to the reuired type
   template <typename ValueType>
   


}}
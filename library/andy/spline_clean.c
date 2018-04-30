#include "tinyspline.h"
#include <stdlib.h>

/*
  clean a spline created by ts_bspline_interpolate_cubic
  so that ts_bspline_derive works correctly
  The spline is cleaned in place
*/

tsError ts_bspline_clean(tsBSpline * spline)
{
   const size_t order = spline->order;
   const size_t dim = spline->dim;
   const tsReal min_useful_distance = (tsReal)1.e-6;// ? for now

   size_t knot_idx = 0U;
   // n.b that number of knots and control points in the spline_copy may change dynamically in the loop
   while(knot_idx < spline->n_knots){
      // count the multiplicity of the knot
      size_t multiplicity = 1U;
      const tsReal cur_knot_value = spline->knots[knot_idx];
      while ( (knot_idx < (spline->n_knots-1)) && ( spline->knots[knot_idx +1] == cur_knot_value) ){
         ++knot_idx;
         ++multiplicity;
      }
      if ( multiplicity == order){
         const size_t idxa = (knot_idx - multiplicity) * dim;
         const size_t idxb = idxa + dim;
         // get and compare the control points
         const tsReal* const ctrl_pointsa = &spline->ctrlp[idxa];
         const tsReal* const ctrl_pointsb = &spline->ctrlp[idxb];
     
         if (ts_ctrlp_dist2(ctrl_pointsa,ctrl_pointsb, dim) < min_useful_distance){

            // create and init a temporary for the modified spline
            tsBSpline mod_spline;
            
            mod_spline.deg = spline->deg;
            mod_spline.order = spline->order;
            mod_spline.dim = spline->dim;
            mod_spline.n_ctrlp = spline->n_ctrlp-1U;
            mod_spline.n_knots = spline->n_knots-1U;

            const size_t n_mod_ctrlp_reals = mod_spline.n_ctrlp * mod_spline.dim;
            const size_t n_mod_raw_bytes_for_real_array = sizeof(tsReal) *  (mod_spline.n_knots + n_mod_ctrlp_reals);
            mod_spline.ctrlp = (tsReal*) malloc(n_mod_raw_bytes_for_real_array);
            if ( mod_spline.ctrlp == NULL){
               return TS_MALLOC;
            }
            // knots array is at an offset in the ctrlpoints array
            mod_spline.knots = mod_spline.ctrlp + n_mod_ctrlp_reals;
            
            // copy control points prior to and including first control point corresponding to current knot
            const tsReal* ctrlp_src = spline->ctrlp;
            tsReal* ctrlp_dest = mod_spline.ctrlp;
            size_t reals_copied = 0U;
            while ( ctrlp_src != ctrl_pointsb){
               *ctrlp_dest++ = *ctrlp_src++;
                ++reals_copied;
            }
            // skip copying 2nd control point corresponding to current knot
            ctrlp_src +=  mod_spline.dim;
            // copy the rest of the control points over
            while ( reals_copied < n_mod_ctrlp_reals ){
                *ctrlp_dest++ = *ctrlp_src++;
                ++reals_copied;
            }
            // copy all knots before current
            const tsReal* const knots_src = spline->knots;
            tsReal* const knots_dest = mod_spline.knots;
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
            // clear the original spline_copy and replace the contents with the new one  
            ts_bspline_free(spline);
            ts_bspline_move(&mod_spline,spline);
         }
      }else{
         ++knot_idx;
      }
   }
   return TS_SUCCESS;
}

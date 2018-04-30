#ifndef SELIG_AEROFOIL_HPP_INCLUDED
#define SELIG_AEROFOIL_HPP_INCLUDED

#include <iostream>
#include <vector>
#include <typeinfo>

   struct  aerofoil{

      bool load(std::string const & filename,std::ostream & e);
      std::string get_name() const {return m_name;}
      std::size_t get_num_coords() const
      {
         return this->m_coords.size();
      }
      quan::two_d::vect<double> get_coord(int i) const
      {
         return this->m_coords.at(i);
      }

      void  add_coord(quan::two_d::vect<double> const & cd)
      {
         this->m_coords.push_back(cd);
      }
      //   void set_name( const char* const name);

      //        quan::two_d::vect<double> get_upper_coord(double const & percent) const;
      //        quan::two_d::vect<double> get_lower_coord(double const & percent) const;
      private:
      std::string m_name;
      std::vector<quan::two_d::vect<double> > m_coords;
   };

#endif // SELIG_AEROFOIL_HPP_INCLUDED

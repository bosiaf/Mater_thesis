#ifndef EPIDEMICS_HPP
#define EPIDEMICS_HPP
#include<iostream>


using namespace std;

namespace epi
{
  class epidemics
  {
    public:
      static unsigned epi_time;
      vector<*host> h_vec;
      const array<const array<const double, 4>, 4 > tmat = 
      {{
      {{0.545, 0.2475, 0.1137, 0.0937}},
      {{0.3836, 0.4357, 0.0891, 0.0915}},
      {{0.2192, 0.1108, 0.4274, 0.2425}},
      {{0.1823, 0.1149, 0.2448, 0.4581}}
      }};


      epidemics(vector<*host> h)
      : epi_time(0), h_vec(h), is_over(false), 
        host_nr(h.size())
      {}
      //Create a vector containing all the strains in a host
      
      ~epidemics()
      {
        h_vec.clear();
#ifdef NDEBUG
        cout << "host vector is clear? " << h_vec.empty() << endl;
#endif
      }
    private:
     bool is_over; 
     unsigned host_nr;

  };
  
}


#endif

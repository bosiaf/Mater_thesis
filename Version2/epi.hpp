#ifndef EPIDEMICS_HPP
#define EPIDEMICS_HPP
#include<iostream>
#include<array>
#include<vector>
#include<random>
#include<string>
#include<memory>
#include<list>


using namespace std;

namespace epi
{
  class epidemics
  {
    public:

      hash<string> = hs_sim;
      bool is_over;//bool to see if the simulation is over
      unsigned epi_time;//the global simulation time
      vector<unique_ptr<host> > h_vec; //vector collecting all the infected hosts

      //Add and delete hosts
      void add_host(const unique_ptr<host> & h)
      {
        h_vec.push_back(move(h));
        ++host_nr;
      }
      
      void delete_host(const unsigned ind)
      {
        swap(h_vec[ind], h_vec[h_vec.size()-1]);
        h_vec.pop_back();
        --host_nr;
      }

      //TRANSITION MATRIX
      const array<const array<const double, 4>, 4 > tmat = 
      {{
      {{0.545, 0.2475, 0.1137, 0.0937}},
      {{0.3836, 0.4357, 0.0891, 0.0915}},
      {{0.2192, 0.1108, 0.4274, 0.2425}},
      {{0.1823, 0.1149, 0.2448, 0.4581}}
      }}; //hardcoded transition matrix
      
      //Hardcoded cumulative sums
      // precalculated sums of the row of the
      // transition matrix.
      const array<const double, 4> cs_a = 
      {
      {tmat[0][0], tmat[0][0] + tmat[0][1], 
       tmat[0][0] + tmat[0][1] + tmat[0][2], 1.}
      }; // A nucleotide cumsum

      const array<const double, 4> cs_g = 
      {
      {tmat[1][0], tmat[1][0] + tmat[1][1], 
       tmat[1][0] + tmat[1][1] + tmat[1][2], 1.}
      }; // G nucleotide cumsum

      const array<const double, 4> cs_c = 
      {
      {tmat[2][0], tmat[2][0] + tmat[2][1], 
       tmat[2][0] + tmat[2][1] + tmat[2][2], 1.}
      }; // C nucleotide cumsum
      
      const array<const double, 4> cs_t = 
      {
      {tmat[3][0], tmat[3][0] + tmat[3][1], 
       tmat[3][0] + tmat[3][1] + tmat[3][2], 1.}
      }; // T nucleotide cumsum
     
      //Matrix containing the 4 cumulative sums
      const array<const array<const double, 4>, 4> cs_tmat =
      {
      {cs_a, cs_g, cs_c, cs_t}
      };
      //END TRANSITION MATRIX

      epidemics(vector<unique_ptr<host> > & h)
      : is_over(false), epi_time(0), h_vec(h), 
        host_nr(h.size())
      {}//constructor
      
      ~epidemics()
      {
        h_vec.clear();//make sure memory is released
#ifdef NDEBUG
        cout << "host vector is clear? " << h_vec.empty() << endl;
#endif
      }
    private:
     unsigned host_nr;//save number of hosts

  };
  
}


#endif

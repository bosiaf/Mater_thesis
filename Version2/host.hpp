#ifndef HOST_HPP
#define HOST_HPP

#include"strain.hpp"
#include<iostream>
#include<string>
#include<vector>

using namespace std;

namespace epi
{
  class host
  {
    public: 
      typedef long count_t;     
      
      //static variable for counting the hosts
      static count_t total_hosts;

      //Develop methods for working with h_vec
      unsigned get_time() const
      { return time; }//get the host creation time
      count_t get_ID() const
      { return ID; }//get the host ID
      count_t get_nr_strains() const
      { return nr_strains; }//returns the number of strains in the host

      //TODO
      //methods to propagate time
      //dynamics according to basic viral model
      bool wi_host_dyn_std(mt19937 & rng, const epi::par & par);
      //evolve sequences
      void evolve();
      
      
      //container for the strains
      vector<*strain> V;
      
      //add and delete a strain
      void add_strain(strain * s) //add a strain and update nr_strains
      {
        V.push_back(s);
        nr_strains = V.size();
      }
      // strain pointer gets destroyed, and the last
      // element is swapped with the destroyed one.
      //this means the order in the vector is not garanteeed
      bool delete_strain(count_t index) //selete a strain without moving too much
      {
        iter_swap(V.begin() + index, V.back());
        V.pop_back();
        --nr_strains;
      }

      count_t healthy_cells;
      
      //constructor
      host(count_t hc, vector<*strain> A)
      : ID(total_hosts++), healthy_cells(hc), V(A),
      nr_strains(A.size()), time(epi_time) 
      {}

      //destructor
      ~host()
      {
        V.clear();
#ifdef NDEBUG
        cout << "V has been cleared? " << V.empty() << endl;
#endif
      }
 
    private:
      count_t nr_strains; //number of strains in the host
      const count_t ID; //ID of the host
      const unsigned time; //time of infection
  }
}

#endif

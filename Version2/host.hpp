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
      //Variable to track number of strains over time
      long unsigned tot_strains;

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
      void evolve(mt19937 & rng, const array<const array <const double> > & tmat,
                  const vector<const unsigned> & SNPs_list, const double p_mut,
                  const unsigned t, const unsigned hash, const par & par);
      
      
      //container for the strains
      vector<strain> V;
      
      //add and delete a strain
      void add_strain(const strain::count_t v, const strain::count_t i_c,
                      const strain::count_t t_c, const strain::count_t l_c,
                      const double fit, const string & s,
                      const unsigned t, const unsigned hash) //add a strain and update nr_strains
      {
        V.emplace_back(v, i_c, t_c, l_c, fit, s, tot_strains, t, hash);
        nr_strains = V.size();
        ++tot_strains;
      }
      // strain pointer gets destroyed, and the last
      // element is swapped with the destroyed one.
      //this means the order in the vector is not garanteeed
      bool delete_strain(const count_t index) //selete a strain without moving too much
      {
        iter_swap(V.begin() + index, V.back());
        V.pop_back();
        --nr_strains;
      }

      //number of healthy cells in the host
      count_t healthy_cells;

      //spread rate
      const double k_spread;
      
      //constructor
      host(const count_t hc, const vector<strain> & A, const double k_spread,
           const unsigned t)
      : tot_strains(0), V(A), healthy_cells(hc), k_spread(k_spread),
      nr_strains(A.size()), ID(total_hosts++), time(t) 
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

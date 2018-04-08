#ifndef HOST_HPP
#define HOST_HPP

#include"strain.h"
#include "par_class.h"
#include "TransitionMatrix.h"
#include<iostream>
#include<string>
#include<vector>

using namespace std;

namespace epi {
class host {
 public:
  typedef unsigned count_t;

  host(count_t hc, const vector<strain>& A, double k_spread, unsigned int t, hash<string>& hs_sim);


  //static variable for counting the hosts
  static count_t total_hosts;
  //Variable to track number of strains over time
  long unsigned tot_strains;

  //Develop methods for working with h_vec
  unsigned get_time() const { return time; }//get the host creation time
  count_t get_ID() const { return ID; }//get the host ID
  count_t get_nr_strains() const { return nr_strains; }//returns the number of strains in the host

  //TODO
  //methods to propagate time
  //dynamics according to basic viral model
  bool wi_host_dyn_std(mt19937& rng, const epi::par& par);
  //evolve sequences
  void evolve(mt19937& rng,
              const vector<unsigned>& SNPs_list, double p_mut,
              unsigned t, const par& par);

  void add_strain(const count_t v, const count_t i_c, const count_t t_c, const count_t l_c, const double fit,
             const string& s,
             const unsigned int t, const unsigned int hash);
  bool delete_strain(const count_t index);
  //container for the strains
  vector<strain> V;

  ~host();

  //number of healthy cells in the host
  count_t healthy_cells;

  //spread rate
  const double k_spread;

 private:
  count_t nr_strains; //number of strains in the host
  const count_t ID; //ID of the host
  const unsigned time; //time of infection
  hash<string>& hs_sim_;

};
}

#endif

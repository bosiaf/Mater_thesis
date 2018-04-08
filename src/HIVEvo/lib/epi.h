#ifndef EPIDEMICS_HPP
#define EPIDEMICS_HPP
#include<iostream>
#include<vector>
#include<random>
#include<string>
#include<cmath>
#include<memory>
#include<list>
#include"host.h"
#include"strain.h"
#include"global_fct.h"
#include"par_class.h"
#include "TransitionMatrix.h"


using namespace std;

namespace epi {

class epidemics {
 public:

  hash<string> hs_sim;
  bool is_over;//bool to see if the simulation is over
  unsigned epi_time;//the global simulation time
  const unsigned diseaseFreePopulation;
  vector<unique_ptr<host> > h_vec; //vector collecting all the infected hosts

  const TransitionMatrix TransitionMatrix_;

  explicit epidemics(unsigned diseaseFreePop);

  ~epidemics();
  void new_host_infection(mt19937& rng, const string& path, const par& par);
  unsigned int next_inf_time(mt19937& rng);
  void add_host(host::count_t hc, const vector<strain>& A, double k_spread);
  void delete_host(unsigned int ind);
 private:
  unsigned host_nr;//save number of hosts
};

}

#endif

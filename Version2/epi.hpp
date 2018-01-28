#ifndef EPIDEMICS_HPP
#define EPIDEMICS_HPP
#include<iostream>
#include<array>
#include<vector>
#include<random>
#include<string>
#include<cmath>
#include<memory>
#include<list>
#include"host.hpp"
#include"strain.hpp"
#include"global_fct.hpp"
#include"par_class.hpp"


using namespace std;

namespace epi
{
  class epidemics
  {
    public:

      bool is_over;//bool to see if the simulation is over
      unsigned epi_time;//the global simulation time
      const unsigned S_df;
      vector<unique_ptr<host> > h_vec; //vector collecting all the infected hosts

      //Add and delete hosts
      void add_host(const host::count_t hc, const vector<strain> & A, 
                    const double k_spread)
      {
        h_vec.push_back( make_unique<host>(hc, A, k_spread, epi_time) );
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

      //CALCULATE NEXT INFECTION TIME
      unsigned next_inf_time(mt19937 & rng)
      {
        double rate = 0;
        //calculate overall rate
        for (unsigned i = 0; i < h_vec.size(); ++i)
        {
          rate += h_vec[i]->k_spread;
        }
        exponential_distribution<> expd(rate*S_df);
        double t = expd(rng);
        cout << "Next infection time is in " << t << " days" << endl;
        return t; 
      }

      //Infect new host
      void new_host_infection(mt19937 rng, const string path, const par & par)
      {
        //Select a host
        //calculate weights
        vector<double> h_w(h_vec.size());
        double tot = 0;
        for (unsigned i = 0; i < h_vec.size(); ++i)
        {
          h_w.push_back(h_vec->k_spread);
        }
        //generate discrete distribution and choose one
        discrete_distribution<unsigned> dd(h_w.begin(), h_w.end());
        unsigned h_n = dd(rng); //index of infecting host
        
        //Sample a strain from previous host "h_n"
        //calculate propensities
        const unsigned str_n = h_vec[h_n]->V.size();
        //times of each strain from beginning
        vector<unsigned> times (str_n);
        /*tr_fit will contain the strain transmission fitnesses.
        * They are calculated as 
        * N_v * exp(-k_fit*p)
        * N_v: virion number of strain v
        * k_fit: tuning parameter for exponential dependency
        * (how bad the fitness decays with time)
        * p: hamming distance, approximated with the JC69 model
        * p = 0.75 - 0.75 * exp(-t*k_mut)
        * k_mut: mutation rate
        */
        vector<double> tr_fit (str_n);
        double k_mut = -log(1-par.pmut);
        for (unsigned i = 0; i < str_n; ++i)
        {
          times[i] = h_vec[h_n]->get_time() - h_vec[h_n]->get_time();
          tr_fit[i] = h_vec[h_n]->V[i].vir/par.vol * exp(-0.75 * par.k_fit *
                      (1 - exp(-times[i]*k_mut)));
        }
        //sample a strain
        discrete_distribution<unsigned> sampler (tr_fit.begin(), tr_fit.end());
        const unsigned str = sampler(rng);
        cout << "Next infecting strain's ID is " << h_vec[h_n]->V[str].get_ID();
        cout << " from host number " << h_vec[h_n]->get_ID() << endl;

        //Construct new strain and host

        //get new spreading rate
        const vector<double> ss_w = {par.kbtw/4, par.kbtw*4};
        discrete_distribution<unsigned> h_sp {4,1};
        //construct new strain vector
        vector<strain> V;
        //construct the host
        e.add_host(par.h0, V, ss_w[h_sp(rng)]);
        //get new virion number
        unsigned n = static_cast<unsigned>(rnorm(v0, v0/3, rng));
        if (n < 1) n = 1;
        e.h_vec.back()->add_strain(n, 0, 0, 0, 1, h_vec[h_n]->V[str].get_sequence(),
                                   epi_time, hs_sim(h_vec[h_n]->V[str].get_sequence()));

      //WHAT TO DO WITH WRITING??
      }

      //Constructor
      epidemics(const vector<unique_ptr<host> > & h, const unsigned S_df)
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

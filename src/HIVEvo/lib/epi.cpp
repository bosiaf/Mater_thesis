//
// Created by iotn_ on 04.04.2018.
//

#include <random>
#include <vector>
#include "par_class.h"
#include "epi.h"
#include "TransitionMatrix.h"

//CALCULATE NEXT INFECTION TIME
unsigned epi::epidemics::next_inf_time(mt19937& rng) {
  double rate = 0;
  //calculate overall rate
  for (auto& i : h_vec) {
    rate += i->k_spread;
  }
  exponential_distribution<> expd(rate * diseaseFreePopulation);
  auto t = static_cast<unsigned>(expd(rng));
  cout << "Next infection time is in " << t << " days" << endl;
  return t;
}

//Infect new host
void epi::epidemics::new_host_infection(std::mt19937& rng, const std::string& path, const par& par) {
  //Select a host
  //calculate weights
  vector<double> h_w(h_vec.size());
  double tot = 0;
  for (auto& i : h_vec) {
    h_w.push_back(i->k_spread);
  }
  //generate discrete distribution and choose one
  discrete_distribution<unsigned> dd(h_w.begin(), h_w.end());
  unsigned h_n = dd(rng); //index of infecting host

  //Sample a strain from previous host "h_n"
  //calculate propensities
  const auto str_n = h_vec[h_n]->V.size();
  //times of each strain from beginning
  vector<unsigned> strainAges(str_n);
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
  vector<double> tr_fit(str_n);
  double k_mut = -log(1 - par.pmut);
  for (unsigned i = 0; i < str_n; ++i) {
    strainAges[i] = h_vec[h_n]->V[i].get_time() - h_vec[h_n]->get_time();
    tr_fit[i] = h_vec[h_n]->V[i].vir / par.vol * exp(-0.75 * par.k_fit *
      (1 - exp(-strainAges[i] * k_mut)));
  }
  //sample a strain
  discrete_distribution<unsigned> sampler(tr_fit.begin(), tr_fit.end());
  const unsigned str = sampler(rng);
  cout << "Next infecting strain's ID is " << h_vec[h_n]->V[str].get_ID();
  cout << " from host number " << h_vec[h_n]->get_ID() << endl;

  //Construct new strain and host

  //get new spreading rate
  const vector<double> ss_w = {par.kbtw / 4, par.kbtw * 4};
  discrete_distribution<unsigned> h_sp{4, 1};
  //construct new strain vector
  vector<strain> V;
  //construct the host
  add_host(par.h0, V, ss_w[h_sp(rng)]);
  //get new virion number
  auto n = static_cast<unsigned>(rnorm(par.v0, par.v0 / 3., rng));
  if (n < 1) n = 1;
  h_vec.back()->add_strain(n, 0, 0, 0, 1, h_vec[h_n]->V[str].get_sequence(),
                           epi_time, hs_sim(h_vec[h_n]->V[str].get_sequence()));

  //WHAT TO DO WITH WRITING??
}

//Add and delete hosts
void epi::epidemics::add_host(const host::count_t hc, const vector<strain>& A,
                              const double k_spread) {
  h_vec.push_back(make_unique<host>(hc, A, k_spread, epi_time, hs_sim));
  ++host_nr;
}

void epi::epidemics::delete_host(const unsigned ind) {
  swap(h_vec[ind], h_vec[h_vec.size() - 1]);
  h_vec.pop_back();
  --host_nr;
}

//Constructor
epi::epidemics::epidemics(const unsigned diseaseFreePop)
  : is_over(false), epi_time(0), h_vec( std::vector<unique_ptr<host>>() ),
    host_nr(0), diseaseFreePopulation(diseaseFreePop),
    TransitionMatrix_(epi::TransitionMatrix()) {}

epi::epidemics::~epidemics() {
  h_vec.clear();//make sure memory is released
#ifdef NDEBUG
  cout << "host vector is clear? " << h_vec.empty() << endl;
#endif
}
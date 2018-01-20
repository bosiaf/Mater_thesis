#include<iostream>
#include<vector>
#include<random>
#include"host.hpp"


bool wi_host_dyn_std(mt19937 & rng, const epi::par * par)
{
  bool result = false;
  bool outcome = false;
  double prob = 0.0;
  const count_t nstr = nr_strains;
  vector<count_t> weight(nstr);
  vector<double> w_norm(nstr);
  vector<double> fit(nstr);
  vector<double> fi_vi(nstr);

  double expk = exp(par.kinf);

  uniform_real_distribution<> u01(0.0,1.0);

  //Vectors for Vose sampler
  vector<double> probs;
  vector<unsigned> alias;
 
  //Variables containing sums of virions and fitness*virions
  count_t v_tot = 0;  
  double sum_fi_vi = .0;

  const count_t hc = healthy_cells;

  //initialize stuff:
  // initialize fit (1 per strain), precalculated exp(k_inf * fit_i) 
  // (1 per strain), number of virions weight (1 per strain) and total number of virions v_sum,
  // as well as the sum S_i^{nstr} (fit_i * n_vir_i). w_norm is a normed weight.
  //fi_vi is for infective fitness.

  for (count_t i = 0; i < nstr; ++i)
  {
    fit[i] = V[i]->fitness;
    weight[i] = V[i]->vir;
    sum_v += weight[i];
    //just calculate infection fitness variable if it is necessary
    if (par.inf_fit_dep)
    {
      eKF[i] = exp(par.kinf*fit[i]);
      fi_vi[i] = weight[i]*fit[i];
      sum_fi_vi += fi_vi[i];
    }
  } 
  //calculate normed weight
  for (count_t i = 0; i < nstr; ++i) w_norm[i] = static_cast<double>(weight[i])/v_sum;
   
  //INFECTION

  //precompute exponential part of the probability
  //it is exp(-K_inf * V_tot) if no dependency to the fitness is given
  // it is exp(-K_inf * Sum(f_i * V_i)) if dependant on fitness.
  double exp_prob = exp(-par.kinf * sum_v);
  if (par.inf_fit_dep) exp_prob = exp(-par.kinf * sum_fi_vi);

  //initialize Vose sampler
  Vose_smpl_init(w_norm, probs, alias, nstr);

  cout << "Starting main infection loop with " << hc;
  cout << " healthy cells and " << nstr << " strains." endl;


  /*Here infection begins with a loop over all healthy cells:
   *A bernoulli trial is done, and a strain is sampled if it is successful
   *A new temporary infected cell is added to that strain (to see if it will evolve)
   */
  for (count_t i = 0; i < hc; ++i)
  {
    prob = 1.0 - exp_prob;
    bernoulli_distribution d(prob);
    
    if (d(rng))
    {
      unsigned i_str = Vose_smpl(probs, alias, nstr, weight, rng, u01);
      V[i_str]->--vir;
      V[i_str]->++temp_cell;
      --healthy_cells;
      --weight[i_str];
      --v_sum;
      //update probability of infection:
      //one less virion
      if (par.inf_fit_dep) exp_prob *= eKF[i_str];
      else exp_prob *= expK;     
    }
    //if the virions reach 0, break from loop
    if (v_sum == 0) break;

    //reset Vose sampler after given steps interval
    if (par.nr_chunks != 1)
    {
      //if the step is the correct one
      if (i % (hc / par.nr_chunks) == 0)
      {
        //recalculate probabilities for each strain
        for (count_t j = 0; j < nstr; ++j) w_norm[j] = static_cast<double>(weight[j])/v_sum;              //and reinitialize the Vose sampler
        Vose_smpl_init(w_norm, probs, alia, nstr);
      }
    }
  }
  
//DYNAMICS OF SPAWNING AND REMOVAL
  
  //Immunocompetence
  count_t tot_c = hc;
  //calculate number of cells
  for (count_t i = 0; i < nstr; ++i) tot_c += V[i].inf_cell;
  //immunocompetence icomp
  double icomp = static_cast<double>(tot_c)/par.h0;
  
  //Loop over strains to propagate dynamics
  for (count_t i = 0; i < nstr; ++i)
  {
    count_t norm_burst = 0;

    //LEFT HERE!!!!
    
    
  }
}


bool wi_host_dyn_lat(mt19937 & rng,)
{}


void evolve()
{}

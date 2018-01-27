#include<iostream>
#include<vector>
#include<random>
#include<functional>
#include<array>
#include"host.hpp"


bool host::wi_host_dyn_dyn(mt19937 & rng, const epi::par & par)
{
  bool result = false;
  double prob = 0.0;
  const count_t nstr = nr_strains;
  vector<count_t> weight(nstr);
  vector<double> w_norm(nstr);
  vector<double> fit(nstr);
  vector<double> fi_vi(nstr);
  vector<count_c> to_elim;

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
    fit[i] = V[i].fitness;
    weight[i] = V[i].vir;
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
      --V[i_str].vir;
      ++V[i_str].temp_cell;
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

  //precalculate mean burst of a single strain
  vector<double> b_norm = rnorm<double>(nstr, par.b_size, par.b_size/4, rng)
  //precalculate latent cell probability of proliferation
  //calculate how many cells are born, and extract from a
  //discrete distribution of what kind they are.
  // (Just execute if lat_max is not 0)
  if (par.lat_max != 0)
  {
    vector<count_t> lat_weigths(nstr);
    vector<count_t> lat_to_add(nstr);
    count_t lat_tot = 0;
    for (count_t i = 0; i < nstr; ++i) 
    {
      lat_weights[i] = V[i].lat_cell;
      lat_tot += V[i].lat_cell;
    }
 
    discrete_distribution<> lat_d(lat_weights.begin(), lat_weights.end());
    const double prob_lat = rtp(par.lat_prol * (1 - lat_tot/par.lat_max));
    const count_t how_many_lat = rbinom(lat_tot, prob_lat, rng);
    for (count_t i = 0; i < how_many_lat; ++i) ++lat_to_add[lat_d(rng)];
  }

  for (count_t i = 0; i < nstr; ++i)
  {
    //check that burst is not negative
    if (b_norm[i] <= 0) b_norm[i] = 0;

    count_t norm_burst = 0;

    //BURST SIZE
    if (V[i].inf_cell)//calculate just if there are infected cells
    {
      //calculate burst with or without fitness dependency
      if (par.burst_fit_dep)
      {
        norm_burst = static_cast<count_t>(b_norm[i]*fit[i]*V[i].inf_cell);
      }
      else
      {
        norm_burst = static_cast<count_t>(b_norm[i]*V[i].inf_cell);
      }
    }
    //INFECTED CELLS REMOVAL
    //VIRIONS REMOVAL
    if (par.ad_im_sys)
    {
      if (par.dic_fit_dep)
      {
        V[i].inf_cell -= rbinom(V[i].inf_cell, rtp(par.dic*icomp/fit[i]), rng);
      }
      else
      {
        V[i].inf_cell -= rbinom(V[i].inf_cell, rtp(par.dic*icomp), rng); 
      }
      if (par.dv_fit_dep)
      {
        V[i].vir -= rbinom(V[i].vir, rtp(par.dv*icomp/fit[i]), rng);
      }
      else
      {
        V[i].vir -= rbinom(V[i].vir, rtp(par.dv*icomp), rng);
      }
    }
    else
    {   
      if (par.dic_fit_dep)
      {
        V[i].inf_cell -= rbinom(V[i].inf_cell, rtp(par.dic/fit[i]), rng);
      }
      else
      {
        V[i].inf_cell -= rbinom(V[i].inf_cell, rtp(par.dic), rng); 
      }
      if (par.dv_fit_dep)
      {
        V[i].vir -= rbinom(V[i].vir, rtp(par.dv/fit[i]), rng);
      }
      else 
      {
        V[i].vir -= rbinom(V[i].vir, rtp(par.dv), rng);
      }
    }   
    //LATENT CELLS
    if (V[i].lat_cell != 0)
    {
      // Proliferate
      V[i].lat_cell += lat_to_add[i];
      // Death
      V[i].lat_cell -= rbinom(V[i].lat_cell, rtp(par.dl, rng)); 
    }
    //spawn new virions
    V[i].vir += norm_burst;

    //if there are no virions and infected cells of a strain,
    //eliminate that strain
    if (!(V[i].vir && V[i].inf_cell && V[i].lat)) to_elim.push_back(i);
    else result = true;
  }
  for (unsigned i = to_elim.size(); i > 0; --i) delete_strain(to_elim[i-1]);

  //HEALTHY CELLS
  healthy_cells -= rbinom(healthy_cells, rtp(par.dhc), rng);
  healthy_cells += floor(rnorm(hc_renew + 0.5, hc_renew/5.0, rng);

  //returns 0 if no viable strain is present
  return result;
}


void host::evolve(mt19937 & rng, const array<const array<const double> > & tmat, 
            const vector<const unsigned> & SNPs_list, const double p_mut, 
            const unsigned time, const par & par)
{
  const unsigned s_sz = V[0].get_size();
  uniform_real_distribution<> ud(0.0, 1.0);
  //get number of mutation per sequence
  //and get where on the single sequence these mutations happen
  //rbinom(seq_length, p.mut, rng)
  //runif_int(0, seq_length - 1, rng)
  //if the index sampled is present in SNPs_list, then write function "get_new_fitness()"
  //to generate the new fitness, taking into account what kind of mutation it is.
  //If it is not a SNP, fitness can either stay the same or decrease, 
  //if it is a SNP, increase
  //fitness through escape in a first approximation (better models will follow)
  //new host infection fitness decreases with mutations.

  //use constant variable to prevent down scaling of loop
  const count_t nstr = V.size();
  count_t hs = 0;
  for (count_t str = 0; str < nstr; ++str)
  {
    //using this variable prevents down-scaling of for loop
    const count_t temp = V[str].temp_cell;
    for (int v = 0; v < temp; ++v)
    {
      --V[str].temp_cell;
      //copy sequence in order not to change the strain 
      //sequence (would affect each virion
      //in the same strain)
      string sq = V[str].get_sequence();
      const count_t n_mut = rbinom(s_sz, p_mut, rng);
      //if there are no mutations...
      if (n_mut == 0)
      {
        //...directly add an infected cell and you're done
        //or a latent cell
        if (ud(rng) > par.inf_to_lat) ++V[str].inf_cell;
        else ++V[str].lat_cell;
      }//If there is a mutation...
      else if (n_mut == 1)
      {
      //...sample a position where this happens
        const unsigned ind = floor(s_sz*ud(rng));
      //record old nucleotide at that position
        const char o_nt = sq[ind];
      //find new nucleotide
        const char subst = evo_nt(tmat, ud(rng), o_nt);
      //substitute the old with the new
        sq[ind] = subst;

      //calculate hash of sequence to look for it
        hs = epidemics::hs_sim(sq);

      //now look for the sequence in the already available sequences
        bool found = false;
        for (count_t i = 0; i < nstr; ++i)
        {
      //How to do string comparison? Maybe hash table would be the best
      //if not use Rabin-Karp. But a lookup-table would be nice
          if (V[i].get_hash() == hs && V[i].get_sequence() == sq)
          {
            if (ud(rng) > par.inf_to_lat) ++V[str].inf_cell;
            else ++V[str].lat_cell;
            found = true;
            break;
          }
        }
        //if the new sequence is not found btw all the strains, create a new one

        if (!found)
        {
          //get a fitness
          double f = V[str].fitness + get_new_fitness(ind, SNPs_list, 
                     o_nt, subst, rng, par);
          if (f < fit_low_cap) f = fit_low_cap;
          //Do not allow for fitness to increase indefinitely
          if (f > fit_high_cap) f = fit_high_cap;
          //instantiate new sequence, strain classes and add a line to host::V.
          if (ud(rng) > par.inf_to_lat) 
          {
            strain * st = new strain(0, 1, 0, 0, f, sq, tot_strains, epidemics::epi_time);
          }
          else  
          {
            strain * st = new strain(0, 0, 0, 1, f, sq, tot_strains, epidemics::epi_time);
          }
          add_strain(st);
        }
      }
      else if (n_mut > 1)
      {
        vector<unsigned> ind = sample_int_wo_repl(s_sz, n_mut, rng);
        
        vector<char> o_nt(n_mut);
        
        for (unsigned m = 0; m < n_mut; ++m)
        {
          o_nt[m] = sq[ind[m]];
          char subst = evo_nt(tmat, ud(rng), sq[ind[m]]);
          sq[ind[m]] = subst;
        }
	
        hs = epidemics::hs_sim(sq);

        bool found = false;

        for (unsigned i = 0; i < V.size(); ++i)
        {
          if (V[i].get_hash() == hs && V[i].get_sequence() == sq)
          {
            if (ud(rng) > par.inf_to_lat) ++V[str].inf_cell;
            else ++V[str].lat_cell;
            found = true;
            break;
          }
        }
        
        //if the new sequence is not found btw all the strains, create a new one
        if (!found)
        {
          //get a fitness
          double f = V[str].fitness;
          for (unsigned m = 0; m < n_mut; ++m)
          {
            f += get_new_fitness(ind[m], SNPs_list, o_nt[m], 
                 sq[ind[m]], rng, par);
          }
          if (f < fit_low_cap) f = fit_low_cap;
          if (f > fit_high_cap) f = fit_high_cap;
          //instantiate new sequence, strain classes and add a line to host::V.
          if (ud(rng) > par.inf_to_lat) 
          {
            strain * st = new strain(0, 1, 0, 0, f, sq, tot_strains, epidemics::epi_time);
          }
          else  
          {
            strain * st = new strain(0, 0, 0, 1, f, sq, tot_strains, epidemics::epi_time);
          }
          add_strain(st);
        }
      }
    }
  }
}

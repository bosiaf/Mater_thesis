//IMPLEMENT PRINT TO FILE AND READ PAPERS ON FITNESS

#include<iostream>
#include<iomanip>
#include<random>
#include<vector>
#include<string>
#include"header.hpp"
#include<omp.h>

using namespace std;
using namespace personal;
unsigned personal::MAX_TSTEP;
int personal::v0, personal::h0, personal::hc_renew, personal::burst_size, personal::nr_chunks;
unsigned personal::seed;
string personal::path_to_tmat, personal::path_output_dyn, personal::path_output_seq, personal::seq0;
double personal::d_hc, personal::d_ic, personal::d_v, personal::k_inf, personal::S_df, personal::k_btw, 
personal::k_mut, personal::fit_snp, personal::fit_change, personal::k_fit, personal::fit_low_cap;
vector<unsigned> personal::SNPs;
vector<long int> personal::weight_not_snp;
vector<double> personal::fit_not_snp;
bool personal::dic_fit_dep, personal::dv_fit_dep, personal::inf_fit_dep, personal::ad_imm_sys;
mt19937 personal::gen(42);
vector<mt19937*> personal::gens;

unsigned host::total = 0;



int main(int argc, char * argv[])
{

	uniform_int_distribution<> u(1,100000);
	for (int i = 0; i < omp_get_max_threads(); ++i)
	{
		gens.push_back(new mt19937(u(gen) * i));
	}
	if(omp_get_max_threads() > 1)
	{
		cout << "Using parallelization OpenMP with " << omp_get_max_threads() << " threads" << endl;
	}
	else
	{
		cout << "No parallelization with OpenMP" << endl;
	}
	//if number of arguments provided is correct
	if (argc == 2){
	//read in parameters, in mm^{-3} day{-1}
		read_pars(argv[1], MAX_TSTEP, path_to_tmat, path_output_dyn, path_output_seq, seq0, SNPs, v0, h0, hc_renew, d_hc, d_ic, burst_size, d_v, k_inf, S_df, k_btw, k_mut, fit_snp, fit_not_snp, weight_not_snp, dic_fit_dep, dv_fit_dep, inf_fit_dep, k_fit, ad_imm_sys, fit_change, fit_low_cap, seed, nr_chunks);
	}
	else //else print a statement with the correct usage
	{
		cout << "Usage: ./CEvo.o path/to/parameters.dat" << endl;
	}
	
	//get probability from rate
	double p_mut = rtp(k_mut);
	//reseed with given seed
	gen.seed(seed);

	vector<vector<double> > tmat;

	if (argc == 2)
	{
		//name/path to the file, container for the transistion matrix, header = T/F
		read_in(path_to_tmat, tmat, true);
	}
	else
	{
		cout << "usage: ./progr path/to/parameters.dat" << endl;
		return 1;
	}
	
	if (nr_chunks < 1)
	{
		cout << "Number of chunks to reinitialize sampler cannot be smaller than 1!!" << endl;
	}
	//Discrete time algorithm
	
	//set the initial strain values
	sequences * s0 = new sequences(seq0, 1);
	
	//initial sequence, one virion, 0 infected cells
	strain * st = new strain(s0, v0, 0, 0, 1, 0);

	//initialize a host specific RNG
	mt19937 g1(u(gen) * (host::total));

	//instantiate a class instance with 1000 healthy cells, the previously defined strain
	//and the local RNG
	host * h1 = new host(h0, st, g1);

	//initialize epidemics with first host and number (constant) of uninfected population
	epidemics e = epidemics(h1, S_df);
	//initialize a first host infection time
	long int t_next_inf = e.next_inf_time(k_btw, gen);
	cout << "Next infection time is " << t_next_inf << endl;

	//begin the epidemics, set a checker for premature end of epidemics
	bool over = 0;
	for (unsigned tstep = 0; tstep < MAX_TSTEP; ++tstep)
	{		
		//check if it is time to infect another host
		while (t_next_inf == 0)
		{
			cout << "New host infected at time " << e.get_time() << "!" << endl;
			e.new_host_infection(gen, path_output_dyn);
			cout << "Virion number of new strain is " << e.get_hosts().back()->get_V().back()->get_vir() << endl;
			t_next_inf = e.next_inf_time(k_btw, gen);
			cout << "New infection time is " << t_next_inf << endl;
		}
		//vector to contain indexes of empty hosts (that will be eliminated)
		vector<unsigned> host_to_drop;
		for(unsigned i = 0 ; i < e.get_hosts().size(); ++i)
		{
			bool is_host_full;
			//check if virions numbers are correct:
			//new healthy cells are created
			e.get_hosts()[i]->set_hc(floor(rnorm(1,hc_renew + 0.5,hc_renew/5.0, gen).back()));

			//infection occurs and death takes its toll
			//with virions and infected cells (strain specific)
			is_host_full = e.get_hosts()[i]->wi_host_inf_death(k_inf, d_ic, d_v, burst_size, gen, tmat, SNPs, p_mut, e.get_time());
			//death comes to take everybody:
			//also target healthy cells
			e.get_hosts()[i]->set_hc(-rbinom(1, e.get_hosts()[i]->get_hc(), rtp(d_hc), gen).back());
			//In the beginning death of healthy cells happened prior to infection, now it is opposite 
			//for no particular reason.
			//is_host_full = e.get_hosts()[i]->wi_host_inf_death(k_inf, d_ic, d_v, burst_size, gen, tmat);
			
			//I created the bool is_host_full to see if I have to eliminate a healed host
			//or a host where the infection did not take place.
			if (!is_host_full)
			{
				//here I gather the indexes of the hosts that became empty
				//to be eliminated afterwards (not now because it could influence the 
				//epidemics in non-model-like ways)
				host_to_drop.push_back(i);
			}	

		}
		//here I eliminate the hosts that healed. No garbage around. 
		//indexing decreases to prevent shrinking aberration in index.
		for (unsigned i = host_to_drop.size(); i > 0; --i)
		{
			e.delete_host(host_to_drop[i-1]);
		}
		//print the epidemics
		//e.print_epidemics();
		e.print_epidemics(path_output_dyn);
		//e.print_seq_epidemics();
		e.print_seq_epidemics(path_output_seq);
		cout << "Time = " << e.get_time() << endl;
		//increment the time by 1
		e.set_time(1);
		//decrement the time to the next infection of a host
		--t_next_inf;

		//decrease the fitness of the strains due to immune system
		e.change_fitness(gen);

		//see if the epidemics is over (everybody healed)
		//and, if so, end the simulation.
		over = e.is_it_over();
		if (over)
		{
			break;
		}
	}

	//Clean up allocated memory
	for (int i = omp_get_max_threads(); i > 0; --i) delete gens[i];
	for (unsigned i = e.get_hosts().size(); i > 0; --i)
	{
		for (unsigned j = e.get_hosts()[i - 1]->get_V().size(); j > 0; --j)
		{
			e.get_hosts()[i - 1]->delete_strain(j - 1);
		}
		e.delete_host(i - 1);
	}
	return 0;
}

#include<iostream>
#include<random>
#include<string>
#include<memory>
#include<functional>
#include<cstdlib>
#include <epi.h>


using namespace std;
using namespace epi;

//Initialize static class members
host::count_t host::total_hosts = 0;//current host number in epidemics
unsigned strain::s_size = 0;

int main(int argc, char * argv[])
{
  //read in files from parameter file if it is present
  if (argc < 2)
  {
    cout << "Please provide parameter file name\nAborting." << endl;
    exit(EXIT_FAILURE);
  }
  
  const string parameter = argv[1]; 
  
  const par a = read_pars(parameter);
  
  if (a.nr_chunks < 1)
  {
    cout << "Number of chunks after which to reinitialize sampler cannot be smaller than 1!" << endl;
    cout << "Exiting program." << endl;
    exit(EXIT_FAILURE);
  }

  a.print_par();//print parameters

  //Initialize static members that need parameters
  strain::set_size( static_cast<unsigned int>(a.seq.size()) );//Strain class sequence size
  //Random number generator
  mt19937 rng(a.seed);

  //Distribution for the host propensities
  //according to the 20/80 rule
  const vector<double> ss_w = {a.kbtw/4, a.kbtw*4};
  discrete_distribution<unsigned> host_spread {4,1}; 

  //initialize host container
  vector<unique_ptr<host> > h;
  vector<strain> V;
  //instantiate epidemics
  epidemics e(a.sdf);//constructor called
  //add the first host to the epidemics
  e.add_host(a.h0, V, ss_w[host_spread(rng)]);
  //add the first strain to the first host
  e.h_vec[0]->add_strain(a.v0, 0, 0, 0, 1, a.seq, e.epi_time, e.hs_sim(a.seq));

  //calculate next infection time
  unsigned t_next_inf = e.next_inf_time(rng);

  //begin the epidemics, set a checker for premature end
  bool is_over = 0;
  for (unsigned tstep = 0; tstep < a.max_tstep; ++tstep)
  {
    //check if new host is to be infected
    while (t_next_inf == 0)
    {
      cout << "New host infected at time " << e.epi_time << "!" << endl;//print current time
      e.new_host_infection(rng, a.path_output_dyn, a);//infect a new host
      t_next_inf = e.next_inf_time(rng);//produce a new infection time
    }
  }



  
  

  return 0;
}

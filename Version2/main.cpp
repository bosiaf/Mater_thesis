#include<iostream>
#include<random>
#include<string>


using namespace std;
using namespace epi;

//Initialize static class members
host::count_t host::total_hosts = 0;//current host number in epidemics


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
  
  a.print_par();//print parameters

  //Initialize static members that need parameters
  unsigned strain::s_size = a.seq.size();//Strain class sequence size
  //Random number generator
  mt19937 rng(a.seed);

  
  

  return 0;
}

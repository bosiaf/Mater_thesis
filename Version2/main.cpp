#include<iostream>
#include<random>
#include<string>


using namespace std;
using namespace epi;

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

  //Random number generator
  mt19937 rng(a.seed);

  return 0;
}

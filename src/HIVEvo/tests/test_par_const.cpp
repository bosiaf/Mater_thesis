#include<iostream>
#include<exception>
#include"../lib/par_class.h"

using namespace std;
using namespace epi;


int main()
{
  const string filename = "parameters_test.dat";
  const epi::par a = read_pars(filename);
  
  try 
  {
    a.v0 = 2;
  }
  catch (exception& e) 
  {
    cout << "Exeption " << e.what() << " happened" << endl;
    exit(EXIT_SUCCESS);
  }
  
  a.print_par();
  return 0;
}

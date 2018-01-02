#include<iostream>
#include"../par_class.hpp"

using namespace epi;


int main()
{
  const string filename = "parameters_test.dat";
  const epi::par a = read_pars(filename);

  a.print_par();
  exit(EXIT_SUCCESS);
  return 0;
}

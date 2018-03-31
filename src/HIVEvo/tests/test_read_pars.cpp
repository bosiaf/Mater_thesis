#include<iostream>
#include"../lib/par_class.h"

using namespace epi;


int main()
{
  const string filename = "parameters_test.dat";
  const epi::par a = read_pars(filename);

  a.print_par();
  exit(EXIT_SUCCESS);
  return 0;
}

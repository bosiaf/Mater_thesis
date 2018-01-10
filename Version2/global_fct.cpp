#include<iostream>
#include"global_fct.hpp"

using namespace std;

int nt_to_i (const char nt)
{
  if (nt == 'A') return 0;
  else if (nt == 'G') return 1;
  else if (nt == 'C') return 2;
  else if (nt == 'T') return 3;
  else
  {
    cout << "Only nucleotides supported are A, C, G and T" << endl;
    cout << "Found " << nt << "." << endl;
    return -1;
  }
}

char i_to_nt (const int nt)
{
  if (nt == 0) return 'A';
  else if (nt == 1) return 'G';
  else if (nt == 2) return 'C';
  else if (nt == 3) return 'T';
  else
  {
    cout << "Only nucleotide codes supported are 0, 1, 2, 3. " << endl;
    cout << "Found " << nt << endl;
    return 'F';
  }
}

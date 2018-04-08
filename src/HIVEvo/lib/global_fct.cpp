#include<iostream>
#include"global_fct.h"

using namespace std;

int epi::nt_to_i (const char nt)
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

char epi::i_to_nt (const int nt)
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

//takes a double transition matrix, a random number and an old nucleotide and
//gives a new nucleotide.
char epi::evo_nt(const double ran_nr, const char old_nt) {
  int nt = nt_to_i(old_nt);
  int new_nt = 0;
  auto oldNucleotideRowCumulativeSum = epi::TransitionMatrix().matrixOfCumulativeSums[nt];
  for (unsigned i = 0; i < 4; ++i) {
    if (ran_nr < oldNucleotideRowCumulativeSum[i]) {
      new_nt = i;
      break;
    }
  }
  return i_to_nt(new_nt);
}

//Vose sampler
void epi::Vose_smpl_init(vector<double> & p, vector<double> & probs,
                    vector<unsigned> & alias, const unsigned size)
{
  unsigned l = 0, s = 0;
  vector<unsigned> large(size);
  vector<unsigned> small(size);

  for (unsigned j = 0; j < size; ++j)
  {
    if (p[j] > (1.0/size))
    {
      large[l] = j;
      ++l;
    }
    else
    {
      small[s] = j;
      ++s;
    }
  }
  while (s != 0 && l != 0)
  {
    --s;
    --l;
    unsigned j = small[s];
    unsigned k = large[l];
    probs[j] = size * p[j];
    alias[j] = k;
    p[k] += p[j] - 1.0 / size;

    if (p[k] > 1.0 / size)
    {
      large[l] = k;
      ++l;
    }
    else
    {
      small[s] = k;
      ++s;
    }
  }
  while (s > 0)
  {
    --s;
    probs[s] = 1;
  }
  while (l > 0)
  {
    --l;
    probs[l] = 1;
  }
}


unsigned epi::Vose_smpl (const vector<double> & probs, const vector<unsigned> & alias,
                    unsigned size, const vector<unsigned> & not_empty, 
                    mt19937 & gen, uniform_real_distribution<> & d)
{
  while(true)
  {
    double u = d(gen)*size;
    unsigned j = floor(u);
    //if the random number lies in the correct bin, the right height
    //and the virion number is not 0, return j or its alias
    if ( (u - j) <= probs[j] )
    {
      if (not_empty[j]) return j;
    }
    else
    {
      if (not_empty[alias[j]]) return alias[j];
    }
  }
}

std::istream& epi::operator>> (std::istream& fin, std::vector<std::vector<double> >& table) {
//clear everything that might be in the table
  table.clear();
  std::vector<double> row;

  while (fin >> row) {
    table.push_back(row);
  }
//return new input
  return fin;
}
//RANDOM DISTRIBUTIONS


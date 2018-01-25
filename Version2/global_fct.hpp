#ifndef GLOBAL_FCT_HPP
#define GLOBAL_FCT_HPP
#include<iostream>
#include<string>
#include<vector>
#include<cmath>
#include<algorithm>

template <class T>
istream& operator>>(istream& fin, vector<T> & row)
{
  //clean the row if it wasn't
  row.clear();
  
  //declare a string buffer for the line
  string l;
  
  //read line of input
  getline(fin, l);

  //this line has now a character as first element, and subsequent elements
  //are numbers apart from header, so I will treat them separately
  stringstream ss(l);
  string element;

  //read row name
  getline(ss, element, ',');
  stringstream field(element);
  string rname;
  field >> rname;

  //read in csv fields without row name
  while (getline (ss, element, ',')
  {
    stringstream field(element);
    //setup variable
    T fld;
    field >> t_prob;
  
    //and now put the read in variable in the row vector
    row.push_back(t_prob);
  }

#ifdef NDEBUG
  for (unsigned i = 0; i < row.size(); ++i)
  {
    std::cout << row[i] << "\t";
  }
#endif

  return fin;
}

template<>
istream& operator>><vector<double> > (istream& fin, vector< vector<> > & table)
{
  //clear everything that might be in the table
  table.clear();
  vector<double> row;
  
  while (fin >> row)
  {
    table.push_back(row);
  }
  //return new input
  return fin;
}


template<class T>
double rtp(const T rate, const double t = 1.)
{
  return double(1 - exp(-t * rate));
}


int nt_to_i (const char nt);
char i_to_n (const int nt);

//takes a double transition matrix, a random number and an old nucleotide and
//gives a new nucleotide.
template<class TM>
char evo_nt (const TM tmat, const double ran_nr, const char old_nt)
{
  char result = '';
  int nt = nt_to_i(old_nt);
  int new_nt = 0;
  for (unsigned i = 0; i < 4; ++i)
  {
    if (random_nr < csum[i])
    {
      new_nt = i;
      break;
    }
  }
  
  return i_to_nt(new_nt);
}


//Vose sampler
void Vose_smpl_init(const vector<double> & p, vector<double> & probs, 
                    vector<unsigned> & alias, const unsigned size);
unsigned Vose_smpl (const vector<double> & probs, const vector<unsigned> & alias, 
                    unsigned size, const vector<unsigned> & not_empty, 
                    mt19937 & gen, uniform_real_distribution<> & d);

//function to calculate new fitness
template<class T>
double get_new_fitness(const T position, const vector<unsigned> & SNPs_vec, 
                       const char old_nt, const char new_nt, mt19937 & rng,
                       const par & par)
{
  //search for the mutation position in the SNP vector.
  //If found, return true
  
  //generate a discrete distribution for getting new fitness.
  dicrete_distribution<> d(par.weight_not_snp.begin(), par.weight_not_snp.end());
  //search
  if (binary_search(SNPs_vec.begin(), SNPs_vec.end(), position)) return par.fit_snp;
  else return par.fit_not_snp[d(rng)];
}

//RANDOM DISTRIBUTIONS
template<class T>
T rnorm(const double mean, const double sd, mt19937 & rng)
{
  normal_distribution<> norm(mean, sd);
  return T(norm(rng));
}

template<class T>
vector<T> rnorm(const unsigned n, const double mean, const double sd, mt19937 & rng)
{
  normal_distribution<> norm(mean,sd);
  vector<T> result(n);

  for (unsigned i = 0; i < n; ++i) result[i] = norm(rng);

  return result;
}

template<class T>
vector<T> rbinom(const unsigned n, const unsigned drows, const double p, mt19937 & rng)
{
  binomial_distribution<> binom(drows, p);
  vector<T> result(n);
  
  if (drows == 0) return result;
  for (unsigned i = 0; i < n; ++i) result[i] = binom(rng);

  return result;
}

template<class T>
T rbinom(const unsigned drows, const double p, mt19937 & rng)
{
  if (drows == 0) return 0;
  binomial_distribution<> binom(drows, p);
  return T(binom(rng));
}

template<class T>
vector<T> sample_int_wo_repl(const T n, const unsigned size, mt19937 & rng)
{
  vector<T> result(size);
  uniform_int_distribution<T> d(0,n-1);
  if (size == 1)
  {
    result[0] = d(rng);
  }
  else
  {
    unsigned i = 0;
    while (i < size)
    {
      T attempt = d(rng);
      bool found = false;
      for (unsigned j = 0; j < result.size(); ++j)
      {
        if (attempt == result[j]) 
        {
          found = true;
          break;
        }
      }
      
      if (!found)
      {
        result[i] = attempt;
        ++i;
      }
    }
  }
  return result;
}

}
#endif

#ifndef GLOBAL_FCT_HPP
#define GLOBAL_FCT_HPP
#include<iostream>
#include<string>
#include<vector>
#include<cmath>

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


template<class T, class C>
double rtp(const T rate, const C t = 1.)
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

#endif

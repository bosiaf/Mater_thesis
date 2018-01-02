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
double rtp(T rate, C t = 1.)
{
  return double(1 - exp(-t * rate));
}


int nt_to_i (char nt)
{
  if (nt == 'A') return 0;
  else if (nt == 'G') return 1;
  else if (nt == 'C') return 2;
  else if (nt == 'T') return 3;
  else
  {
    cout << "Only nucelotides supported are A, G, C, T" << endl;
    cout << "Found " << nt << "." << endl;
    return -1;
  }
}

//takes a double transition matrix, a random number and an old nucleotide and
//gives a new nucleotide.
//TODO:
//tmat hardcoding
//cszums hardcoded
//end evo_nt code
char evo_nt (vector<vector<double> > tmat, double ran_nr, char old_nt)
{
  char result = '';
  int nt = nt_to_i(old_nt);
  double csum[4];
}

#endif

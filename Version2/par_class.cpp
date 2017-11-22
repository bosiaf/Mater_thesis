#include<iostream>
#include<string>
#include<vector>
#include"par_class.hpp"


using namespace std;

epi::par epi::read_pars(const string file)
{
  string seq_in;
  ifstream file_in(file);

  if(!file_in.is_open())
  {
    cout << "Cannot find specified file: " << file << endl;
    exit(EXIT_FAILURE);
  }
  
  //create buffers
  stringstream inp = "";
  string s = "";
 
  //read lines from file into s and eliminate # comments
  while (getline(file_in, s))
  {
    //eliminate comment if present
    s = s.substr(0, s.file('#'));
    //add rest to the stream
    if (s != "") inp << s << endl;
  }

  //close connection to file
  file_in.close();

  //create string buffers for vectors
  string snps = "", fit_not_snp_str = "", weight_not_snp_str = "";

  //declare temporary variables, and define
  string pat_to_tmat = "", path_output_dyn = "", path_output_seq = "",
  seq = "";
  unsigned max_tstep = 0, vol = 0, v0 = 0, h0 = 0, hc_ren = 0,
  b_size = 0, seed = 0, nr_chunks = 0;
  vector<unsigned> SNPs, weoght_not_snp;
  vector<double> fit_not_snp;
  double dhc = 0., dic = 0., dv = 0., kinf = 0., sdf = 0., kbtw = 0.,
  kmut = 0., fit_snp = 0., k_fit = 0., fit_change = 0., fit_low_cap = 0.,
  fit_high_cap = 0.;
  //print that you're about to read in the parameters
  cout << string(30, '*') << endl;
  cout << "I am reading in the parameters\n" << endl;

  getline()
  inp >> max_tstep;
  inp >> path_to_tmat;
  inp >> path_output_dyn;
  inp >> path_output_seq;
  inp >> seq;
  inp.ignore();
  getline(inp,snps,'\n');
}

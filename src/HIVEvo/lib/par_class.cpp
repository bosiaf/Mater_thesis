#include<iostream>
#include<string>
#include<vector>
#include<algorithm>
#include<fstream>
#include<sstream>
#include<iterator>
#include"par_class.h"
#include"global_fct.h"


using namespace std;

const epi::par epi::read_pars(const string& filename)
{
  ifstream file_in(filename);

  if(!file_in.is_open())
  {
    cout << "Cannot find specified file: " << filename << endl;
    exit(EXIT_FAILURE);
  }
  
  //create buffers
  stringstream inp("");
  string s;
 
  //read lines from file into s and eliminate # comments
  while (getline(file_in, s))
  {
    //eliminate comment if present
    s = s.substr(0, s.find('#'));
    //add rest to the stream
    if (s != "") inp << s << endl;
  }

  //close connection to file
  file_in.close();

  //create string buffers for vectors
  string snps, fit_not_snp_str, weight_not_snp_str;

  //declare temporary variables, and define
  string path_to_tmat, path_output_dyn, path_output_seq,
    seq, seq_in;
  unsigned max_tstep = 0, v0 = 0, h0 = 0, hc_ren = 0,
  b_size = 0, lat_max = 0, seed = 0, nr_chunks = 0, sdf = 0;
  vector<unsigned> SNPs, weight_not_snp;
  vector<double> fit_not_snp;
  double vol = 0., dhc = 0., dic = 0., dv = 0., dl = 0.,  
  kinf = 0., kbtw = 0., inf_to_lat = 0., lat_act = 0., k_fit = 0.,
  lat_prol = 0., kmut = 0., fit_snp = 0., fit_change = 0., fit_low_cap = 0.,
  fit_high_cap = 0.;
  bool dic_fit_dep = false, dv_fit_dep = false, burst_fit_dep = false, 
       inf_fit_dep = false, ad_imm_sys = false, parallel = false, 
       seq_per_time = false, seq_print = false;
  
  //print that you're about to read in the parameters
  cout << string(50, '*') << endl;
  cout << "I am reading in the parameters\n" << endl;

  //use s as temporary buffer
  s = "";

  getline(inp, s);
  max_tstep = stoul(s);

  //getline(inp, path_to_tmat); Now Hardcoded
  getline(inp, path_output_dyn);
  getline(inp, path_output_seq);
  getline(inp, seq_in);
  getline(inp, snps);
  getline(inp, s);
  vol = stod(s);

  getline(inp, s);
  v0 = stoul(s);

  getline(inp, s);
  h0 = stoul(s);
 
  getline(inp, s);
  hc_ren = stoul(s);

  getline(inp, s);
  dhc = stod(s);
 
  getline(inp, s);
  dic = stod(s);

  getline(inp, s);
  dl = stod(s);

  getline(inp, s);
  b_size = stoul(s);

  getline(inp, s);
  dv = stod(s);

  getline(inp, s);
  kinf = stod(s);

  getline(inp, s);
  inf_to_lat = stod(s);

  getline(inp, s);
  lat_act = stod(s);

  getline(inp, s);
  lat_max = stoul(s);

  getline(inp, s);
  lat_prol = stod(s);

  getline(inp, s);
  sdf = stoul(s);

  getline(inp, s);
  kbtw = stod(s);

  getline(inp, s);
  kmut = stod(s);

  getline(inp, s);
  fit_snp = stod(s);

  getline(inp, fit_not_snp_str);
  getline(inp, weight_not_snp_str);

  getline(inp, s);
  dic_fit_dep = stoi(s);
  getline(inp, s);
  dv_fit_dep = stoi(s);
  getline(inp, s);
  burst_fit_dep = stoi(s);
  getline(inp, s);
  inf_fit_dep = stoi(s);
  getline(inp, s);
  ad_imm_sys = stoi(s);
  getline(inp, s);
  fit_change = stod(s);
  getline(inp, s);
  fit_low_cap = stod(s);
  getline(inp, s);
  fit_high_cap = stod(s);
  getline(inp, s);
  k_fit = stod(s);
  getline(inp, s);
  seed = stoul(s);
  getline(inp, s);
  nr_chunks = stoul(s);
  getline(inp, s);
  parallel = stoi(s);
  getline(inp, s);
  seq_print = stoi(s);
  getline(inp, s);
  seq_per_time = stoi(s);

  if(getline(inp, s)) cout << s << endl;
  else cout << "Buffer read in successfully and completely." << endl;

  cout << string(50, '*') << endl;
  //first check format of snps, looking for the dash
  unsigned long found = snps.find('-');
  //if it's found, then generate the range, 
  //if not, then put the desired positions.
  if (found != string::npos)
  {
    //read in start and end of range delimited by dash
    unsigned long start = 0, end = 0;
    string element;
    stringstream read_range(snps);
    getline(read_range, element, '-');
    stringstream field(element);
    field >> start;
    getline(read_range, element);
    stringstream field2(element);
    field2 >> end;

    //generate the value range
    for (unsigned i = start; i <= end; ++i)
    {
      SNPs.push_back(i);
    }
  }
  else //if the input is not given with the dash
  {
    istringstream buffer(snps);
    copy(istream_iterator<unsigned> (buffer),
         istream_iterator<unsigned> (), back_inserter(SNPs));
  }

  //read in other vectors
  istringstream buffer2(fit_not_snp_str);
  copy(istream_iterator<double> (buffer2), 
       istream_iterator<double> (), back_inserter(fit_not_snp));
  istringstream buffer3(weight_not_snp_str);
  copy(istream_iterator<unsigned> (buffer3),
       istream_iterator<unsigned> (), back_inserter(weight_not_snp));

  cout << "Check input: Sequence" << endl;
  //read in the sequence from a file
  ifstream s_in(seq_in);
  if(!s_in.is_open())
  {
    cout << "Cannot find specified sequence file: " << seq_in << endl;
    exit(EXIT_FAILURE);
  }

  s_in >> seq;

  s_in.close();
  //transform to have all uppercase characters
  transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
  
  //CHECK INPUT
  //sequence validity:
  const string valid_chars("ACGT");
  found = seq.find_last_not_of(valid_chars);
  if (found != string::npos)
  {
    cout << "Invalid initial sequence at position " << found;
    cout << "!! Exiting program." << endl;
    cout << "Invalid character is " << seq[found] << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    cout << "Sequence was accepted." << endl;
  }

  //Now transform the constants with the volume
  h0 *= vol;
  hc_ren *= vol;
  kinf /= vol;
  lat_max *= vol;
  const double pmut = rtp(kmut);
  

  const par a(path_output_dyn, path_output_seq,
                     seq, max_tstep, v0, h0, hc_ren, b_size, lat_max,
                     seed, nr_chunks, sdf, SNPs, weight_not_snp, fit_not_snp,
                     vol, dhc, dic, dv, dl, inf_to_lat, lat_act, 
                     lat_prol, k_fit, kinf, kbtw, pmut, fit_snp,
                     fit_change, fit_low_cap, fit_high_cap,
                     dic_fit_dep, dv_fit_dep, inf_fit_dep, ad_imm_sys,
                     parallel, seq_per_time, seq_print);

  return a;
}


void epi::par::print_par() const
{
  cout << string(50, '*') << endl;
  
  cout << "Time steps to simulate: " << max_tstep << endl;
  //cout << "Path to transition matrix " << path_to_tmat << endl;
  cout << "Path to output folder of dynamics files: " << path_output_dyn << endl;
  cout << "Path to output folder of sequence files: " << path_output_seq << endl;
  cout << "Initial sequence: \n" << seq << endl;
  cout << "Sequence size is: " << seq.size() << endl;
  cout << "Location of SNPs: ";
  for (unsigned int SNP : SNPs) cout << SNP << " ";
  cout << endl;
  cout << "The simulation volume is " << vol << " mm^3" << endl;
  cout << "Initial virions: " << v0 << endl;
  cout << "Initial healthy cells: " << h0 << endl;
  cout << "Renewal rate of healthy cells: " << hc_ren << endl;
  cout << "Death rate of a single healthy cell: " << dhc << endl;
  cout << "Death rate of a single infected cell: " << dic << endl;
  cout << "Burst size: " << b_size << endl;
  cout << "Death rate of a single virion: " << dv << endl;
  cout << "Infection rate constant: " << kinf << endl;
  cout << "Disease free susceptibles: " << sdf << endl;
  cout << "Btw. host infection rate constant: " << kbtw << endl;
  cout << "Change of infection fitness with age of strain exponential parameter: " << k_fit << endl;
  cout << "Mutation probability of a nt after a day: " << pmut << endl;
  cout << "Fitness increase for SNP mutation: " << fit_snp << endl;
  cout << "Fitness changes for non-SNP mutation: ";
  for (double i : fit_not_snp) cout << i << " ";
  cout << endl;
  cout << "Probability weight of change in fitness: ";
  for (unsigned int i : weight_not_snp) cout << i << " ";
  cout << endl;
}

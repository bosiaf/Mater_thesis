#include<iostream>
#include<iomanip>
#include<vector>
#include<list>
#include<string>
#include<cmath>
#include<fstream>
#include<sstream>
#include"header.hpp"
#include<random>
#include<algorithm>
#include<iterator>
#include<utility>
#include<chrono>
#include<thread>
#include<omp.h>
#include<sys/stat.h>


using namespace std;
using namespace personal;

//Overloaded operators
istream& personal::operator >> (istream& fin, vector<double> & row)
{
	//if there was something in row, eliminate it
	row.clear();

	//declare a string buffer to contain a line
	string l;

	//read a line of the input
	getline(fin, l);

	//this line has now a character as first element, and subsequent elements are numbers
	//apart from header, so i will treat them separately
	stringstream ss(l);
	string element;

	//read row name
	getline(ss, element, ',');
	stringstream field(element);
	string rname;
	field >> rname;

	//read in the comma separated field without the row name
	while (getline(ss, element, ','))
	{
		stringstream field(element);
		//set up the variable for storing
		double t_prob = 0.0;
		field >> t_prob;

		//and put the read in variable in the row vector
		row.push_back(t_prob);
	}

	/*for (unsigned i = 0; i < row.size(); ++i)
	{
		cout << row[i] << "\t";
	}
	cout << endl;
	*/
	//return the rest of the input from file
	return fin;

}

istream& personal::operator >> (istream& fin, vector<string> & row)
{
	//if there was something in row, eliminate it
	row.clear();

	//declare a string buffer to contain a line
	string l;

	//read a line of the input
	getline(fin, l);

	//this line has now a character as first element, and subsequent elements are numbers
	//apart from header, so i will treat them separately
	stringstream ss(l);
	string element;

	//read in the comma separated field
	while (getline(ss, element, ','))
	{
		stringstream field(element);
		//set up the variable for storing
		string cname;
		field >> cname;

		//and put the read in variable in the row vector
		row.push_back(cname);
	}
	/*cout << "Header is ";
	for (unsigned i = 0; i < row.size(); ++i)
	{
		cout << row[i] << "\t";
	}
	cout << endl;
	*/
	//return the rest of the input from file
	return fin;

}


istream& personal::operator >> (istream& fin, vector<vector<double> > & table)
{
	//clear everything that might be in the table
	table.clear();
	vector<double> row;

	while (fin >> row)
	{
		table.push_back(row);
	}
	//return rest of input
	return fin;
}



//GLOBAL FUNCTIONS

int personal::nt_to_i(char nt)
{
	if (nt == 'A') return 0;
	else if (nt == 'G') return 1;
	else if (nt == 'C') return 2;
	else if (nt == 'T')	return 3;
	else
	{
		cout << "Only nucleotides supported are A, G, C, T" << endl;
		return -1;
	}
}

char personal::i_to_nt(int nt)
{
	if (nt == 0) return 'A';
	else if (nt == 1) return 'G';
	else if (nt == 2) return 'C';
	else if (nt == 3) return 'T';
	else
	{
		cout << "Only nucleotides supported are A, G, C, T" << endl;
		return 'F';
	}
}

//functions to implement random number generators, samplers,...
vector<double> personal::rexp(int n, double lam, mt19937 & genr)
{
	exponential_distribution<> exp(lam);
	vector<double> result(n);

	if (lam >= 0)
	{
		for (unsigned i = 0; i < n; ++i)
		{
			result[i] = exp(genr);
		}
	}
	else
	{
		for (unsigned i = 0; i < n; ++i)
		{
			result[i] = 0;
		}
	}
	return result;
}

vector<long int> personal::rbinom(int n, unsigned drows, double p, mt19937 & genr)
{
	binomial_distribution<> binom(drows, p);
	vector<long int> result(n);

	for (unsigned i = 0; i < n; ++i)
	{
		result[i] = binom(genr);
	}
	return result;
}

vector<double> personal::rnorm(int n, double mean, double sd, mt19937 & genr)
{
	normal_distribution<> norm(mean, sd);
	vector<double> result(n);

	for (unsigned i = 0; i < n; ++i)
	{
		result[i] = norm(genr);
	}
	return result;
}

vector<long int> personal::runif_int(int n, unsigned from, unsigned to, mt19937 & genr)
{
	uniform_int_distribution<unsigned> unif(from, to);
	vector<long int> result(n);

	for (unsigned i = 0; i < n; ++i)
	{
		result[i] = unif(genr);
	}
	return result;
}

vector<bool> personal::rber(int n, double p, mt19937 & genr)
{
	bernoulli_distribution ber(p);
	vector<bool> result(n);

	for (unsigned i = 0; i < n; ++i)
	{
		result[i] = ber(genr);
	}
	return result;
}

long int personal::smpl_weight(vector<long int> w, mt19937 & gen, uniform_real_distribution<> d)
{
	size_t w_sz = w.size();
	vector<double> w_csum(w_sz);
	double r = d(gen);
	w_csum[0] = w[0];
	for (unsigned i = 1; i < w_sz; ++i)	w_csum[i] = w_csum[i - 1] + w[i];
	for (int i = 0; i < w_sz; ++i)
	{
		w_csum[i] /= w_csum.back();
	}
	for (unsigned i = 0; i < w_sz; ++i)
	{
		if (w_csum[i] > r) return i;
	}
	cout << "not found" << endl;
	return -1;
}

long int personal::smpl_weight(vector<double> w, mt19937 & gen, uniform_real_distribution<> d)
{
	//discrete_distribution<long int> dd(w.begin(), w.end());
	//return dd(gen);
	size_t w_sz = w.size();
	vector<double> w_csum(w_sz);
	double r = d(gen);
	w_csum[0] = w[0];

	for (unsigned i = 1; i < w_sz; ++i)	w_csum[i] = w_csum[i - 1] + w[i];
	//vectorize this
	for (int i = 0; i < w_sz; ++i)
	{
		w_csum[i] /= w_csum.back();
	}
	for (unsigned i = 0; i < w_sz; ++i)
	{
		if (w_csum[i] > r) return i;
	}
	cout << "not found" << endl;
	return -1;
}

//rand and init functions from Vose alias method
//"A linear algorithm for generating random numbers
//with a given distribution"

//Init function
void personal::Vose_smpl_init(vector<double> p, vector<double> & probs, vector<int> & alias, int size)
{
	int l = 0, s = 0;
	vector<int> large(size);
	vector<int> small(size);

	for (int j = 0; j < size; ++j)
	{
		if (p[j] > (1.0 / size))
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
		int j = small[s];
		int k = large[l];
		probs[j] = size*p[j];
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

//rand function as described in the abovementioned paper
int personal::Vose_smpl(vector<double> probs, vector<int> alias, int size, vector<long> not_empty, mt19937 & gen, uniform_real_distribution<> d)
{
	while (true)
	{
		double u = d(gen)*size;
		int j = floor(u);
		//if the random number lies in the correct bin, the right height and 
		//the virion number is not 0, return j or its alias
		if ((u - j) <= probs[j])
		{
			if (not_empty[j]) return j;
		}
		else
		{
			if (not_empty[alias[j]]) return alias[j];
		}
	}
}

vector<long int> personal::sample_int_wo_repl(long int n, unsigned size, mt19937 & gen)
{
	vector<long int> result;
	uniform_int_distribution<> d(0, n - 1);
	if (size == 1)
	{
		result.push_back(d(gen));
	}
	else
	{
		unsigned i = 0;
		while (i < size)
		{
			long int attempt = d(gen);
			bool found = false;
			for (unsigned j = 0; j < result.size(); ++j)
			{
				if (attempt == result[j]) found = true;
			}

			if (!found)
			{
				result.push_back(attempt);
				++i;
			}
		}
	}
	return result;
}

//evolution of a nucleotide (random evolution), takes as argument a random number from 0-1 unif dist
char personal::evo_nt(vector<vector<double> > tmat, double random_nr, char old_nt)
{
	char result;
	int nt = nt_to_i(old_nt);
	double csum[4];
	csum[0] = tmat[nt][0];
	//initialize to 3 to prevent rounding errors
	int new_nt = 2;

	for (unsigned i = 1; i < 4; ++i)
	{
		csum[i] = csum[i - 1] + tmat[nt][i];
	}
	for (unsigned i = 0; i < 4; ++i)
	{
		if (random_nr < csum[i])
		{
			new_nt = i;
			break;
		}
	}

	result = i_to_nt(new_nt);

	return result;
}

//function to read in the transmission matrix from a file. Matrix must 
//have first line with the names of the 4 bases as Character (header). The following rows
//are composed of a column of character (rowname) and 4 doubles.
//The transition matrix will just contain the numbers, whereas
//tmat[1][2] indicates the first row and second column (as R).
void personal::read_in(string file, vector<vector<double> > &trans_mat, bool header)
{
	vector<string> h;
	ifstream file_in(file);
	if (!file_in.is_open())
	{
		cout << "Cannot find specified file: " << file << endl;
		return;
	}

	//if there's a header, read it in	
	file_in >> h;

	//then read the rest of the matrix:
	file_in >> trans_mat;

	//close the connection to the file
	file_in.close();

	cout << "\n****************************************" << endl;
	cout << "Reading Transition Matrix:" << endl;
	for (unsigned i = 0; i < h.size(); ++i)
	{
		cout << "\t" << h[i];
	}
	cout << endl;
	for (unsigned i = 0; i < trans_mat.size(); ++i)
	{
		cout << h[i] << "\t";
		for (unsigned j = 0; j < trans_mat.size(); ++j)
		{
			cout << trans_mat[i][j] << "\t";
		}
		cout << endl;
	}
	cout << "****************************************\n" << endl;

	if (trans_mat[0][0] != 0)
	{
		cout << "Matrix has diagonal elements different from 0" << endl;
		cout << "Do not worry, I take care of that automaticaly" << endl;
		for (unsigned i = 0; i < trans_mat.size(); ++i)
		{
			//set diagonal elements to 0
			double sum = 0.0;
			trans_mat[i][i] = 0;
			//renormalize matrix so that sum(row) = 1
			for (unsigned j = 0; j < trans_mat.size(); ++j)
			{
				//calculate sum of row
				sum += trans_mat[i][j];
			}
			for (unsigned j = 0; j < trans_mat.size(); ++j)
			{
				//divide probability by the sum of  the row to get a total prob of 1
				trans_mat[i][j] /= sum;
			}
		}
	}

	cout << "\n****************************************" << endl;
	cout << "Normalized Transition Matrix:" << endl;
	for (unsigned i = 0; i < h.size(); ++i)
	{
		cout << "\t\t" << h[i];
	}
	cout << endl;
	for (unsigned i = 0; i < trans_mat.size(); ++i)
	{
		cout << h[i] << "\t\t";
		for (unsigned j = 0; j < trans_mat.size(); ++j)
		{
			cout << trans_mat[i][j];
			if (j == i) cout << "\t\t";
			else cout << "\t";

		}
		cout << endl;
	}
	cout << "****************************************\n" << endl;

}

void personal::read_pars(string file, unsigned & max_tstep, string & path_to_tmat, string & path_output_dyn, string & path_output_seq, string & seq, vector<unsigned> & SNPs, int & v0, int & h0, int & hc_ren, double & dhc, double & dic, int & b_size, double & dv, double & kinf, double & sdf, double & kbtw, double & kmut, double & fit_snp, vector<double> & fit_not_snp, vector<long int> & weight_not_snp, bool & dic_fit_dep, bool & dv_fit_dep, bool & inf_fit_dep, double & k_fit, bool & ad_imm_sys, double & fit_change, double & fit_l_c, unsigned & seed, int & nr_chunks)
{
	ifstream file_in(file);
	if (!file_in.is_open())
	{
		cout << "Cannot find specified file: " << file << endl;
		return;
	}

	//create buffers
	stringstream inp;
	string s;

	//read lines from file to the string s
	while (getline(file_in, s))
	{
		//find the comment character and eliminate it
		s = s.substr(0, s.find('#'));
		//if it was not commented, then add it to the stream
		if (s != "")
		{
			inp << s << endl;
		}
	}

	//close the connection to the file
	file_in.close();

	//create string buffers for vectors
	string snps, fit_not_snp_str, weight_not_snp_str;

	//Print that you're readin in the parameters
	cout << "\n***********************************" << endl;
	cout << "I am reading in the parameters.\n" << endl;
	//read in the file line per line, vectors are first read in the string buffer
	inp >> max_tstep;
	inp >> path_to_tmat;
	inp >> path_output_dyn;
	inp >> path_output_seq;
	inp >> seq;
	inp.ignore();
	inp.ignore();
	getline(inp, snps, '\n');
	inp >> v0;
	inp >> h0;
	inp >> hc_ren;
	inp >> dhc;
	inp >> dic;
	inp >> b_size;
	inp >> dv;
	inp >> kinf;
	inp >> sdf;
	inp >> kbtw;
	inp >> kmut;
	inp >> fit_snp;
	inp.ignore();
	inp.ignore();
	getline(inp, fit_not_snp_str, '\n');
	getline(inp, weight_not_snp_str, '\n');
	inp >> dic_fit_dep;
	inp >> dv_fit_dep;
	inp >> inf_fit_dep;
	inp >> k_fit;
	inp >> ad_imm_sys;
	inp >> fit_change;
	inp >> fit_l_c;
	inp >> seed;
	inp >> nr_chunks;

	//first look which one of the two versions of SNPs was given, 
	//checking if the dash is there.
	size_t found = snps.find('-');
	//if it is not there
	if (found != string::npos)
	{
		//read in start and end by splitting the line at the dash
		int start, end;
		string element;
		stringstream read_range(snps);
		getline(read_range, element, '-');
		stringstream field(element);
		field >> start;
		getline(read_range, element, '-');
		stringstream field2(element);
		field2 >> end;
		//put in SNPs all loci from a to b
		for (unsigned i = start; i <= end; ++i)
		{
			SNPs.push_back(i);
		}
	}
	else //if not
	{
		//read all the elements on the line in SNPs
		istringstream buffer(snps);
		copy(istream_iterator<unsigned>(buffer), istream_iterator<unsigned>(), back_inserter(SNPs));
	}

	//read in the vectors in the right vectors
	istringstream buffer1(fit_not_snp_str);
	copy(istream_iterator<double>(buffer1), istream_iterator<double>(), back_inserter(fit_not_snp));
	istringstream buffer2(weight_not_snp_str);
	copy(istream_iterator<long int>(buffer2), istream_iterator<long int>(), back_inserter(weight_not_snp));

	cout << "Time steps to simulate: " << max_tstep << endl;
	cout << "Path to transistion matrix: " << path_to_tmat << endl;
	cout << "Path to output folder for dynamics files: " << path_output_dyn << endl;
	cout << "Path to output folder for sequence files: " << path_output_seq << endl;
	cout << "Initial sequence: " << seq << endl;
	cout << "Sequence size is: " << seq.size() << endl;
	cout << "Location of SNPs: ";
	for (unsigned i = 0; i < SNPs.size(); ++i) cout << SNPs[i] << " ";
	cout << endl;
	cout << "Initial healthy cells (hc): " << h0 << endl;
	cout << "Renewal rate of hc: " << hc_ren << endl;
	cout << "Death rate of a hc: " << dhc << endl;
	cout << "Death rate of a infected cell: " << dic << endl;
	cout << "Burst size: " << b_size << endl;
	cout << "Death rate of virions: " << dv << endl;
	cout << "Infection rate constant: " << kinf << endl;
	cout << "Disease free susceptibles: " << sdf << endl;
	cout << "Btw. host infection rate constant: " << kbtw << endl;
	cout << "Mutation rate " << kmut << endl;
	cout << "Fitness increase for SNP mutation: " << fit_snp << endl;
	cout << "Fitness changes for non-SNP mutation: ";
	for (unsigned i = 0; i < fit_not_snp.size(); ++i) cout << fit_not_snp[i] << " ";
	cout << endl;
	cout << "Probability weight of change in fitness: ";
	for (unsigned i = 0; i < weight_not_snp.size(); ++i) cout << weight_not_snp[i] << " ";
	cout << endl;
	cout << "Infected Cells fitness dependency switch set to: " << dic_fit_dep << endl;
	cout << "Virion death rate fitness dependency switch set to: " << dv_fit_dep << endl;
	cout << "Infection rate fitness dependency switch set to: " << inf_fit_dep << endl;
	cout << "Burst size fitness multiplicative dependency: " << k_fit << endl;
	cout << "Adaptive immune system switch set to: " << ad_imm_sys << endl;
	cout << "Change of fitness over time is: " << fit_change << endl;
	cout << "Fitness low cap is: " << fit_l_c << endl;
	cout << "Seed for random number generator: " << seed << endl;
	cout << "Number of chunks in the parallel infection region: " << nr_chunks << endl;

	cout << "************************************************************" << endl;

}


//for description, see header file
double personal::rtp(double rate)
{
	return (1 - exp((-1)*rate));
}

bool personal::fileExists(const string & file)
{
	struct stat buf;

	return stat(file.c_str(), &buf) != -1;
}

//sequences class member functions
sequences::sequences(string s_i, double fit)
{
	sequence = s_i;
	fitness = fit;
	s_size = s_i.size();
}

string sequences::get_sequence()
{
	return sequence;
}

void sequences::set_sequence(string s)
{
	sequence = s;
}

void sequences::set_sequence(char nt, unsigned ind)
{
	sequence[ind] = nt;
}

double sequences::get_fitness()
{
	return fitness;
}

void sequences::set_fitness(double fit)
{
	fitness = fit;
}

unsigned sequences::get_size()
{
	return s_size;
}

//virion class member functions
virion::virion(sequences * s)
{
	pars = s;
}

//strain class member functions
strain::strain(sequences * s, unsigned id, unsigned t)
{
	seq = s;
	vir = 0;
	icell = 1;
	tcell = 0;
	ID = id;
	time = t;
}

strain::strain(sequences * s, long int v, long int ic, long int tc, unsigned id, unsigned t)
{
	vir = v;
	icell = ic;
	seq = s;
	tcell = tc;
	ID = id;
	time = t;

}

sequences * strain::get_sequence()
{
	return seq;
}

void strain::change_fitness(double diff)
{
	seq->set_fitness(seq->get_fitness() + diff);
}

long int strain::get_vir()
{
	return vir;
}

void strain::set_vir(int diff)
{
	vir += diff;
	if (vir < 0)
	{
		vir = 0;
	}
}

long int strain::get_icell()
{
	return icell;
}

void strain::set_icell(int diff)
{
	icell += diff;
	if (icell < 0)
	{
		icell = 0;
	}
}

unsigned strain::get_ID()
{
	return ID;
}

long int strain::get_tcell()
{
	return tcell;
}

void strain::set_tcell(int diff)
{
	tcell += diff;
	if (tcell < 0) tcell = 0;
}

unsigned strain::get_time()
{
	return time;
}

//host member functions

host::host(long int hc, strain * s, mt19937 gen_local)
{
	healthy_cells = hc;
	V.push_back(s);
	ID = total;
	gen_l = gen_local;
	++total;
	l_time = 0;
	n_str = 1;

}

long int host::get_ltime()
{
	return l_time;
}

void host::set_ltime(int diff)
{
	l_time += diff;
}

void host::new_str()
{
	++n_str;
}

mt19937 host::get_rng()
{
	return gen_l;
}

//returns 0 if host is empty, 1 if host is not empty (at least 1 strain is still present)
//rewritten function. Previously it was implemented with a double for-loop over all strains 
//and virions, and they would try to infect healthy cells. A bias towards the first strains in the
//container was introduced, as they have more targerts at disposal than the last in line.
//So I reversed the loops, in a first attemp (over all virions, then over all strains), 
//but it was overly complicated and did not eliminate the bias completely.
//This version loops over the healthy cells as they try to become infected and then, if
//this happens, the infected strain is determined in a Gillespie similar way. Elegant, simple.
//Previous version kept below in a huge comment for historic reasons.
bool host::wi_host_inf_death(double k, double d_infc, double d_vir, long int burst, mt19937 & gen, vector<vector<double> > tmat, vector<unsigned> SNPs_list, double p_mut, unsigned t)
{
	//number of strains in host
	unsigned N_STR = V.size();
	//bool outcome;
	bool result = 0;
	double prob;
	vector<double> fit(N_STR);
	//precomputed exp(k*fit_i)
	vector<double> eKF(N_STR);
	vector<long> weight(N_STR);
	long int v_sum = 0;
	double sum_fi_vi = 0;
	vector<double> cumsum_v(N_STR + 1);
	vector<double> weight_prob(N_STR);
	uniform_real_distribution<> u01(0.0, 1.0);
	vector<unsigned> to_elim;
	vector<vector<unsigned> > parallel_elim(omp_get_max_threads());
	vector<double> probs(N_STR);
	vector<int> alias(N_STR);
	vector< vector <long> > weights_par(omp_get_max_threads());
	//make a copy of healthy_cells to run the loop on, if not it will change the
	//loop variable and analyze less than the total of healthy cells when
	//--healthy_cells is called.
	long int hc = healthy_cells;

	//initialize stuff
	// initialize fit (1 per strain), precaclulated exp(k*fit_i) (1 per strain)
	// number of virions weight ( 1 per strain) and total number of virion v_sum, as well as
	// the sum S_i^{N_STR}(fit_i * n_vir_i) Weight_prob is the samew as weight but normed at 1
	//CAN BE MADE PARALLEL if slow
	for (int i = 0; i < N_STR; ++i)
	{
		fit[i] = V[i]->get_sequence()->get_fitness();
		eKF[i] = exp(k*fit[i]);
		weight[i] = V[i]->get_vir();// *fit[i];
		v_sum += weight[i];
		sum_fi_vi += fit[i] * weight[i];
	}
	for (int i = 0; i < N_STR; ++i)
	{
		weight_prob[i] = static_cast<double>(weight[i])/v_sum;
	}
	
	//initialize the parallel weights vector as a 2D vector of [number_threads]*[N_STRAINS]
	//so that each thread has its own copy to play with. in the end it will be updated globally.
	//wêights_par causes false sharing
	for (size_t i = 0; i < omp_get_max_threads(); ++i)
	{
		weights_par[i].resize(N_STR);
	}			

	double exp_prob = exp(-k * sum_fi_vi);
	cout << "a" << endl;
	Vose_smpl_init(weight_prob, probs, alias, N_STR);
	cout << "I am before healthy cell infection loop with " << hc << " healthy cells and " << N_STR << " strains." << endl;
	//BOTTLENECK, MUST FIND WAY OF SPEEDING IT UP
	//TRY DIVIDING IT IN CHUNCKS
	// Set weight as private and initialized (firstprivate?)
	//Idea: 
	//1) Find a sensible size determination for the chunks (fixed, dynamic?)
	//2) Find a way to update the properties after chunk is done, for example
	//store the used up virions in a thread-private vector and then sum up all the entries to
	//get an overall update.
	//3) Beware of race conditions!
	//4) Random number generator
	//5) empty vector must be synchronized btw all threads
	
//Chunk size determination
	int chunk_size = hc / nr_chunks;

	for (int chunk = 0; chunk < nr_chunks; ++chunk)
	{
		//reset the temporary weight change vector and the hc count
		for (int i = 0; i < omp_get_max_threads(); ++i)
		{
			fill(weights_par[i].begin(), weights_par[i].end(), 0);
		}
		//update the probability
		prob = 1.0 - exp_prob;

		int nvir = rbinom(1, chunk_size, prob, gen).back();
		int i_str;

		#pragma omp parallel for
		for (int i = 0; i < nvir; ++i)
		{
			//here a strain index is sampled according to the abundance of the relative
			//virions frequencies. This is an approximate algorithm: weights are not changed
			//within the chunk to allow seamless parallelization. This is not a big deal,though.
			//weight: vector containing the number of virions at the beginning of the chunk.
			//probs, alias: used for Vose sampler, no need to know that: read the source code
			//of the sampler.
			//gens: vector containing many mt19937 RNG
			//u01: uniform 0-1 distribution
			//weights_par: 2D vector of dimension [#threads]*[#strains] to store the number of
			//virions that infected for updating the probabilities after chunk end.
			//
			//A virion infects, and a temporarly infected cell is created
			//The weight vector is NOT yet adaptes, to allow easier parallelization
			//at cost of some accuracy.
			//
			int i_str = Vose_smpl(probs, alias, N_STR, weight, *gens[omp_get_thread_num()], u01);
			//Now update the temporary dimensions.
			++weights_par[omp_get_thread_num()][i_str];	
		}
			//update the quantities
		for (int str = 0; str < N_STR; ++str)
		{
			long temp_v = 0;
			//get number of virions of strain "str"
			for (int j = 0; j < omp_get_max_threads(); ++j)
			{
				temp_v += weights_par[j][str];
			}

			if (temp_v > weight[str]) temp_v = weight[str];
			V[str]->set_vir(-temp_v);
			V[str]->set_tcell(temp_v);
			healthy_cells -= temp_v;
			weight[str] -= temp_v;
			v_sum -= temp_v;
			for (int k = 0; k < temp_v; ++k)
			{
				exp_prob *= eKF[str]; 
			}				
		}
		//if the number of virions reaches 0, break.
		//since the highest number of virions attainable each chunk is smaller equal
		//v_sum, this never becomes negative.
		//temp_v is always <= sum(weight) = sum_v
		if (v_sum == 0) break;
		//for (int i = 0; i < N_STR; ++i) weight_prob[i] = static_cast<double>(weight[i])/v_sum;
		//Vose_smpl_init(weight_prob, probs, alias, N_STR);			
	}
	
	//execute further only if virions are not depleted
	if (v_sum != 0)
	{
		//infect the last cells (from the last chunk to the end)
		prob = 1.0 - exp_prob;
		
		int n_vir = rbinom(1, hc - (nr_chunks * chunk_size), prob, gen).back();

		for (int i = 0; i < omp_get_max_threads(); ++i)
		{
			fill(weights_par[i].begin(), weights_par[i].end(), 0);
		}
		int i_str;
		#pragma omp parallel for
		for (int i = nr_chunks*chunk_size; i < hc; ++i)
		{
				int i_str = Vose_smpl(probs, alias, N_STR, weight, *gens[omp_get_thread_num()], u01);
				++weights_par[omp_get_thread_num()][i_str];
		}
		for(int str = 0; str < N_STR; ++str)
		{
			long temp_v = 0;
			for (int j = 0; j < omp_get_max_threads(); ++j) temp_v += weights_par[j][str];
			if (temp_v > weight[str]) temp_v = weight[str];
			V[str]->set_vir(-temp_v);
			V[str]->set_icell(temp_v);
			healthy_cells -= temp_v;
			weight[str] -= temp_v;
			v_sum -= temp_v;
		}
	}

/*
*	#pragma omp parallel for firstprivate(weight)
*	for (size_t i = 0; i < hc; ++i)
*	{
*		prob = 1.0 - exp_prob;
*		outcome = rber(1, prob, gen).back();
*
*		if (outcome)
*		{
*			//here a strain index is sampled according to the abundance of the relative virions
*			//frequencies
*			//weights: vector containing the cumulative sum of the frequency in [0,1] of the virions
*			//cumsum_v: vector containing the cumulative sum of the absolute abundance of the virions
*			int i_str;
*			//Vose_smpl_init(weight_prob, probs, alias, N_STR);
*			i_str = Vose_smpl(probs, alias, N_STR, weight_loc, gen, u01);
*			//i_str = smpl_weight(weight, gen, u01);
*			//one virion infects, and a temporarly infected cell is created.
*			//The weight vector is adapted, the number of healthy cells decreases of one,
*			//and the total number of virions v_sum too.
*			//Also change the cumulative sum
*			
*			//let weight be updated globally only if a part reaches 0. 
*			//Do atomic updates: first let the thread that reached 0 update the global 
*			//copy, then let the others update the local copies with the global one.
*			//But is the global copy even necessary?
*			//
*			//if (!(--weight_loc[i_str]))
*			//{
*			//	#pragma omp atomic
*			//	weight[i_str] = 0;
*			//} 
*			//
*			//if ( !(--wheight_loc[i_str]) )
*			//{
*			//	for (int i = 0; i < omp_max_num_threads(); ++i)
*			//	{
*			//		if(omp_thread_id() == i) weight_loc[i_str] = 0;
*			//	}
*			//}
*			V[i_str]->set_vir(-1);
*			V[i_str]->set_tcell(1);
*			--healthy_cells;
*			--weight[i_str];
*			--v_sum;
*			exp_prob *= eKF[i_str];
*			*weight_prob = weight;
* #pragma omp parallel for
*			for (int i = 0; i < N_STR; ++i)
*			{
*				weight_prob[i] /= v_sum;
*			}
*		}
*		if (v_sum == 0)	break;
*	}
*/
	/*for (size_t i = 0; i < hc; ++i)
	{
		prob = 1.0 - exp_prob;
		outcome = rber(1, prob, gen).back();

		if (outcome)
		{
			vector<double> probs(N_STR);
			vector<int> alias(N_STR);
			//here a strain index is sampled according to the abundance of the relative virions
			//frequencies
			//weights: vector containing the cumulative sum of the frequency in [0,1] of the virions
			//cumsum_v: vector containing the cumulative sum of the absolute abundance of the virions
			int i_str;
			Vose_smpl_init(weight_prob, probs, alias, N_STR);
			i_str = Vose_smpl(probs, alias, N_STR, gen, u01);
			//i_str = smpl_weight(weight, gen, u01);
				//one virion infects, and a temporarly infected cell is created.
				//The weight vector is adapted, the number of healthy cells decreases of one,
				//and the total number of virions v_sum too.
				//Also change the cumulative sum

			V[i_str]->set_vir(-1);
			V[i_str]->set_tcell(1);
			--healthy_cells;
			--weight[i_str]; //-= fit[i_str];
			--v_sum;
			exp_prob *= eKF[i_str];
			weight_prob = weight;
#pragma omp parallel for
			for (int i = 0; i < N_STR; ++i)
			{
				weight_prob[i] /= v_sum;
			}
		}
		if (v_sum == 0)	break;
	}*/
	//Now evolution can take place on the temporarly infected cells  tcell 
	//(infected, but not yet sure from which sequence yet)
	evolve(gen, tmat, SNPs_list, p_mut, t);

	//death comes to gather everyone (infected cells are destroyed, then virion destroyed,
	//then virions are created.
	//infected cells
	//
	//Get total number of cells tot_c
	long int tot_c = hc;
	for (int i = 0; i < V.size(); ++i)
	{
		tot_c += V[i]->get_icell();
	}
	double immunocompetence = static_cast<double>(tot_c) / h0;


	//PARALLEL, NEED TO THINK ABOUT RANDOM NUMBER GENERATOR
	//MIGHT USE VECTOR OF RNG FOR NOT TOO MANY NUMBERS
	//CHANGED GEN TO GENL
	//Loop over old strains ( not the newest just created ones.)
#pragma omp parallel for
	for (int i = 0; i < N_STR; ++i)
	{
		mt19937 genl = *(gens[omp_get_thread_num()]);
		int norm_burst = 0;
		//Calculate the burst size with k_fit as fitness dependency factor
		if (V[i]->get_icell() != 0)
		{
			norm_burst = static_cast<int>(rnorm(1, k_fit * fit[i] * burst*V[i]->get_icell(), k_fit * fit[i] * burst*V[i]->get_icell() / 3.0, genl).back());
		}
		else
		{
			norm_burst = 0;
		}
		if (norm_burst < 0) norm_burst = 0;

		//rtp() generates probabilities from rates for a time of 1 day
		//switch on the option selected for in the parameter file
		if (dic_fit_dep)
		{
			if (ad_imm_sys)
			{
				V[i]->set_icell(-rbinom(1, V[i]->get_icell(), rtp(d_infc*immunocompetence / fit[i]), genl).back());
			}
			else
			{
				V[i]->set_icell(-rbinom(1, V[i]->get_icell(), rtp(d_infc / fit[i]), genl).back());
			}
		}
		else
		{
			if (ad_imm_sys)
			{
				V[i]->set_icell(-rbinom(1, V[i]->get_icell(), rtp(d_infc*immunocompetence), genl).back());
			}
			else
			{
				V[i]->set_icell(-rbinom(1, V[i]->get_icell(), rtp(d_infc), genl).back());
			}
		}

		//virions that die and are born, fitness does its magic here.
		//switch on the option selected for in the parameter file
		if (dv_fit_dep)
		{
			if (ad_imm_sys)
			{
				V[i]->set_vir(-rbinom(1, V[i]->get_vir(), rtp(d_vir*immunocompetence / fit[i]), genl).back());
			}
			else
			{
				V[i]->set_vir(-rbinom(1, V[i]->get_vir(), rtp(d_vir / fit[i]), genl).back());
			}
		}
		else
		{
			if (ad_imm_sys)
			{
				V[i]->set_vir(-rbinom(1, V[i]->get_vir(), rtp(d_vir*immunocompetence), genl).back());
			}
			else
			{
				V[i]->set_vir(-rbinom(1, V[i]->get_vir(), rtp(d_vir), genl).back());
			}
		}
		//burst virions
		V[i]->set_vir(norm_burst);
		//if both virion number and infected cell number of a strain are = 0
		//eliminate the strain.
		if (!(V[i]->get_vir() && V[i]->get_icell()))
		{
			//to_elim.push_back(i);
			parallel_elim[omp_get_thread_num()].push_back(i);
		}
		else
		{
#if _OPENMP >= 200505
#pragma omp critical
			{
				result = 1;
			}
#else
#pragma omp atomic
			result = 1;
#endif
		}
	}

	//strains stored for elimination are here deleted.
	//Go through the array once again and look for the empty strains to eliminate.
	for (int i = 0; i < parallel_elim.size(); ++i)
	{
		for (int j = 0; j < parallel_elim[i].size(); ++j)
		{
			to_elim.push_back(parallel_elim[i][j]);
		}
	}

	sort(to_elim.begin(), to_elim.end());

	for (unsigned i = to_elim.size(); i > 0; --i)
	{
		delete_strain(to_elim[i - 1]);
	}

	return result;
}

//THIS FUNCTION IS BIASED AND SHOULD BE CHANGED: FIRST STRAINS HAVE MORE HEALTHY CELLS AT DISPOSAL FOR INFECTION THAN LAST STRAINS, version w/ virion loop and not healthy cells loop
/*bool host::wi_host_inf_death(double k, double d_infc, double d_vir, long int burst, mt19937 & gen)
{
	bool outcome;
	bool result = 0;
	double prob;
	double fit;
	//vector containing the strains to be eliminated
	vector<unsigned> to_elim;
	//rate = k1*T*fit -> if units are day^{-1}, P(t < 1 day) = 1 - e^{-k1*T*fit}
	for (unsigned i = 0; i < V.size(); ++i)
	{
		fit = V[i]->get_sequence()->get_fitness();
		prob = 1 - exp(-k*healthy_cells*fit);
		unsigned V_amount = V[i]->get_vir();

		for (unsigned j = 0; j < V_amount; ++j)//break if healthy cells are =0
		{
			outcome = rber(1, prob, gen).back();
			//if the infection succedes (decided for each virion by a bernoulli experiment)
			if (outcome)
			{
				//decrement healthy cell and virus by 1 and increase infected cell of the same
				//strain i by 1
				--healthy_cells;
				V[i]->set_vir(-1);
				V[i]->set_icell(1);
				//and modify the probability
				prob = 1 - exp(-k*healthy_cells*fit);
			}
			if(prob == 0)
			{
				break;
			}
		}
		//death comes to gather everyone
		//infected cells
		V[i]->set_icell(-rbinom(1, V[i]->get_icell(), rtp(d_infc), gen).back());
		//virions, that die and are born, fitness does its magic here.
		V[i]->set_vir(-rbinom(1, V[i]->get_vir(), rtp(d_vir/fit), gen).back());
		V[i]->set_vir(rnorm(1, fit*burst*V[i]->get_icell(), fit*burst*V[i]->get_icell()/5.0, gen).back());
		if (!(V[i]->get_vir() && V[i]->get_icell()))
		{
			to_elim.push_back(i);
		}
		else
		{
			result = 1;
		}
	}
	for (unsigned i = 0; i < to_elim.size(); ++i)
	{
		delete_strain(to_elim[i]);
	}
	return result;
}
*/


unsigned host::get_ID()
{
	return ID;
}

void host::add_line(strain * s)
{
	V.push_back(s);
}

void host::delete_strain(unsigned index)
{
	vector<strain*>::iterator it = V.begin();
	if (index >= V.size())
	{
		cout << "Index is out of bounds, boy " << endl;
	}
	else
	{
		delete(*(it + index));
		V.erase(it + index);
	}
}

long int host::get_hc()
{
	return healthy_cells;
}

void host::set_hc(long int diff)
{
	healthy_cells += diff;
	if (healthy_cells < 0)
	{
		healthy_cells = 0;
	}
}

double host::get_new_fitness(unsigned position, vector<unsigned> SNPs_list, char old_nt, char new_nt, mt19937 & gen)
{
	//search for the mutation position in the SNP list. If found, return true
	bool is_SNP = binary_search(SNPs_list.begin(), SNPs_list.end(), position);
	//generate a discrete distribution with weights read in from the parameters
	discrete_distribution<> d(weight_not_snp.begin(), weight_not_snp.end());
	//separate the cases of being in a SNP position or being in a neutral position.
	if (is_SNP)
	{
		return (fit_snp);
	}
	else
	{
		return (fit_not_snp[d(gen)]);
	}
}
//MAYBE RECORD TIME OF NEW EVOLUTION
void host::evolve(mt19937 & gen, vector<vector<double> > tmat, vector<unsigned> SNPs_list, double p_mut, unsigned time)
{
	unsigned s_sz = V[0]->get_sequence()->get_size();
	uniform_real_distribution<> ud(0.0, 1.0);
	//get number of mutation per sequence
	//and get where on the single sequence these mutations happen
	//rbinom(1, seq_length, p.mut).back()
	//runif_int(1, 0, seq_length - 1).back()
	//if the index sampled is present in SNPs_list, then write function "get_new_fitness()"
	//to generate the new fitness, taking into account what kind of mutation it is.
	//If it is not a SNP, fitness can either stay the same or decrease, if it is a SNP, increase
	//fitness through escape in a first approximation (better models will follow)
	//use constant variable to prevent down scaling of loop
	int N_STR = V.size();
	for (int str = 0; str < N_STR; ++str)
	{
		//using this variable prevents down-scaling of for loop
		long int temp = V[str]->get_tcell();
		for (int v = 0; v < temp; ++v)
		{
			V[str]->set_tcell(-1);
			//copy sequence in order not to change the strain sequence (would affect each virion
			//in the same strain)
			string sq = V[str]->get_sequence()->get_sequence();
			int n_mut = rbinom(1, s_sz, p_mut, *(gens[omp_get_thread_num()])).back();
			//if there are no mutations...
			if (n_mut == 0)
			{
				//...directly add an infected cell and you're done
				V[str]->set_icell(1);
			}//If there is a mutation...
			else if (n_mut == 1)
			{
				//...sample a position where this happens
				long int ind = sample_int_wo_repl(s_sz, 1, *(gens[omp_get_thread_num()])).back();
				//long int ind = sample(s_sz, 1, gen, ud).back();
				//record old nucleotide at that position
				char o_nt = sq[ind];
				//find new nucleotide
				char subst = evo_nt(tmat, ud(*(gens[omp_get_thread_num()])), sq[ind]);
				//substitute the old with the new
				sq[ind] = subst;
				//now look for the sequence in the already available sequences
				bool found = false;
				for (int i = 0; i < V.size(); ++i)
				{
					//How to do string comparison? Maybe hash table would be the best
					//if not use Rabin-Karp. But a lookup-table would be nice
					//FIX THIS
					if (V[i]->get_sequence()->get_sequence() == sq)
					{
						V[i]->set_icell(1);
						found = true;
						break;
					}
				}
				//if the new sequence is not found btw all the strains, create a new one
				if (!found)
				{
					//get a fitness
					double f = V[str]->get_sequence()->get_fitness() + get_new_fitness(ind, SNPs_list, o_nt, subst, *(gens[omp_get_thread_num()]));
					if (f < fit_low_cap) f = fit_low_cap;
					//instantiate new sequence, strain classes and add a line to host::V.
					new_str();
					sequences * s0 = new sequences(sq, f);
					strain * st = new strain(s0, n_str, time);
					add_line(st);
				}
			}
			else if (n_mut > 1)
			{
				vector<long int> ind = sample_int_wo_repl(s_sz, n_mut, *(gens[omp_get_thread_num()]));
				//vector<long int> ind = sample(s_sz, n_mut, gen, ud);
				vector<char> o_nt(n_mut);
				for (unsigned m = 0; m < n_mut; ++m)
				{
					o_nt[m] = sq[ind[m]];
					char subst = evo_nt(tmat, ud(*(gens[omp_get_thread_num()])), sq[ind[m]]);
					sq[ind[m]] = subst;
				}

				bool found = false;

				for (unsigned i = 0; i < V.size(); ++i)
				{
					//How to do string comparison? Maybe hash table would be the best
					//if not use Rabin-Karp. But a lookup-table would be nice
					//FIX THIS
					if (V[i]->get_sequence()->get_sequence() == sq)
					{
						V[i]->set_icell(1);
						found = true;
						break;
					}
				}
				//if the new sequence is not found btw all the strains, create a new one
				if (!found)
				{
					//get a fitness
					double f = V[str]->get_sequence()->get_fitness();
					for (unsigned m = 0; m < n_mut; ++m)
					{
						f += get_new_fitness(ind[m], SNPs_list, o_nt[m], sq[ind[m]], *(gens[omp_get_thread_num()]));
					}
					if (f < fit_low_cap) f = fit_low_cap;
					//instantiate new sequence, strain classes and add a line to host::V.
					new_str();
					sequences * s0 = new sequences(sq, f);
					strain * st = new strain(s0, n_str, time);
					add_line(st);
				}
			}
		}
	}
}

vector<strain*> host::get_V()
{
	return V;
}

void host::print()
{
	cout << "Total number of healthy cells in the host: " << healthy_cells << endl;
	cout << "Strain\tVirion\tInfected cell" << endl;

	for (unsigned i = 0; i < V.size(); ++i)
	{
		cout << V[i]->get_ID() << "\t" << V[i]->get_vir() << "\t" << V[i]->get_icell() << endl;
	}
}

void host::print_seq()
{
	for (unsigned i = 0; i < V.size(); ++i)
	{
		cout << "Strain\t" << V[i]->get_ID() << endl;
		cout << "Sequence\t" << (V[i]->get_sequence())->get_sequence() << endl;
		cout << "Fitness\t" << (V[i]->get_sequence())->get_fitness() << endl;
	}
}

//epidemics class member functions

epidemics::epidemics(host * h, long int s)
{
	time = 0;
	hosts.push_back(h);
	S = s;
}

epidemics::epidemics()
{
	time = 0;
}

long int epidemics::get_time()
{
	return time;
}

void epidemics::set_time(long int t_diff)
{
	time += t_diff;
}

unsigned epidemics::get_S()
{
	return S;
}

void epidemics::set_S(long int diff)
{
	S += diff;
	if (S < 0)
	{
		S = 0;
	}
}

//samples an exponentially distributed time for the next host infection
long int epidemics::next_inf_time(double k, mt19937 & gen)
{
	if (k > 0)
	{
		return static_cast<long>(rexp(1, k*S, gen).back());
	}
	else
	{
		return -1;
	}
}

//samples a host for infection and a strain therein to infect a new host
void epidemics::new_host_infection(mt19937 & gen, string path)
{
	vector<host*>::iterator it = hosts.begin();
	int run = runif_int(1, 0, (hosts.size() - 1), gen).back();
	host * inf = *(it + run);

	vector<strain*> inf_strains = inf->get_V();
	//vector<strain*>::iterator it_s = inf_strains.begin();
	//probs is a container for the weigths of the viruses
	vector<long int> probs;
	for (unsigned i = 0; i < inf_strains.size(); ++i)
	{
		//if you want to change the weighting, here it is happening
		probs.push_back(inf_strains[i]->get_vir());
	}
	//now sample the sequences in the selected host, weighting them
	//by the number of virions (and the fitness?)
	discrete_distribution<long int> sampler(probs.begin(), probs.end());
	long int index = sampler(gen);
	cout << "Index of next infected strain is " << index << endl;
	strain * strain_infecting = inf_strains[index];
	//now create the new host (must have new everything, a pointer to the
	//already existing sequences will not suffice!)
	//A pointer to an old sequence would make the two sequences impossible to
	//evolve independently
	sequences * a = new sequences(strain_infecting->get_sequence()->get_sequence(), 1);
	//the new host starts with v0 virions, but could actually start with a random number of virions.
	int n = static_cast<int>(rnorm(1, v0, v0 / 3, gen).back());
	if (n < 0) n = 1;
	strain * b = new strain(a, n, 0, 0, 1, time);
	mt19937 gen2(host::total + 1);
	host * h_new = new host(h0, b, gen2);
	add_host(h_new);

	//create the stream to print to the file
	string filename = string(path + "Infection_history");
	ofstream fout;

	if (!fileExists(filename + ".dat"))
	{
		fout.open(filename + ".dat");

		if (!fout.is_open())
		{
			cout << "Could not establish connection with file " << filename + ".dat" << endl;
			return;
		}
		//add a nice informative header just the first timestep
		fout << "Time\tInfecting\tInfected\tStrain" << endl;
	}
	else
	{
		fout.open(filename + ".dat", ios::app); //append if the file already exists

		if (!fout.is_open())
		{
			cout << "Could not establish connection with file " << filename + ".dat" << endl;
			return;
		}
	}

	//print the data
	fout << time << "\t" << run << "\t" << hosts.size() - 1 << "\t" << index << endl;


	//close the connection
	fout.close();
}


void epidemics::add_host(host * h)
{
	hosts.push_back(h);
}

void epidemics::delete_host(unsigned index)
{
	vector<host*>::iterator it = hosts.begin();
	if (index >= hosts.size())
	{
		cout << "Index out of bounds, 'tard" << endl;
	}
	else
	{
		delete(*(it + index));
		hosts.erase(it + index);
	}
}

vector<host*> epidemics::get_hosts()
{
	return hosts;
}

void epidemics::print_epidemics()
{
	cout << "\n****************************************\n" << endl;
	cout << "Epidemics has gone on for " << time << " days" << endl;
	vector<host*>::iterator it_h = hosts.begin();
	for (; it_h != hosts.end(); ++it_h)
	{
		cout << "Host number " << (*it_h)->get_ID() << endl;
		(*it_h)->print();
	}
}

void epidemics::print_epidemics(string path)
{

	vector<host*>::iterator it_h = hosts.begin();
	for (; it_h != hosts.end(); ++it_h)
	{
		//instantiate a stream for the output of V and healthy cells
		string filename = string(path + "host_" + to_string((*it_h)->get_ID()));
		ofstream fout, fout_hc;

		if (!fileExists(filename + ".dat"))
		{
			fout.open(filename + ".dat");
			fout_hc.open(filename + "_healthy_cells.dat");

			if (!fout.is_open())
			{
				cout << "Could not establish connection with file " << filename + ".dat" << endl;
				return;
			}
			if (!fout_hc.is_open())
			{
				cout << "Could not establish connection with file " << filename + "_healthy_cells.dat" << endl;
				return;
			}
			//add a nice informative header just the first timestep
			fout << "Time\tStrain\tVirion\tInfected cell" << endl;
			fout_hc << "Time\tHealthy cells\tTotal cells\tTotal virions" << endl;
		}
		else
		{
			fout.open(filename + ".dat", ios::app);
			fout_hc.open(filename + "_healthy_cells.dat", ios::app);

			if (!fout.is_open())
			{
				cout << "Could not establish connection with file " << filename + ".dat" << endl;
				return;
			}

			if (!fout_hc.is_open())
			{
				cout << "Could not establish connection with file " << filename + "_healthy_cells.dat" << endl;
				return;
			}
		}



		//initialize counters for virions and infected cells
		long int tot_vir = 0, tot_icell = 0;
		//and now print the data
		for (unsigned i = 0; i < (*it_h)->get_V().size(); ++i)
		{
			fout << time << "\t" << (*it_h)->get_V()[i]->get_ID() << "\t";
			fout << (*it_h)->get_V()[i]->get_vir() << "\t";
			fout << (*it_h)->get_V()[i]->get_icell() << endl;
			tot_vir += (*it_h)->get_V()[i]->get_vir();
			tot_icell += (*it_h)->get_V()[i]->get_icell();
		}


		fout_hc << time << "\t" << (*it_h)->get_hc() << "\t" << (*it_h)->get_hc() + tot_icell << "\t" << tot_vir << endl;
		//close the connections
		fout.close();
		fout_hc.close();
	}
}

void epidemics::print_seq_epidemics()
{
	vector<host*>::iterator it_h = hosts.begin();
	for (; it_h != hosts.end(); ++it_h)
	{
		cout << "Host number " << (*it_h)->get_ID() << endl;
		(*it_h)->print_seq();
	}
}

void epidemics::print_seq_epidemics(string path)
{
	vector<host*>::iterator it_h = hosts.begin();
	for (; it_h != hosts.end(); ++it_h)
	{
		string filename = path + "host_" + to_string((*it_h)->get_ID()) + "_seq.dat";
		ofstream fout;

		if (!fileExists(filename))
		{
			fout.open(filename);

			//check that fout is correctly opened
			if (!fout.is_open())
			{
				cout << "Could not establish connection with file " << filename << endl;
				return;
			}

			//print nice informative header the first time ever ever
			fout << "Time\tStrain\tSequence\tFitness" << endl;
		}
		else
		{
			fout.open(filename, ios::app);

			//check that fout is correctly opened
			if (!fout.is_open())
			{
				cout << "Could not establish connection with file " << filename << endl;
				return;
			}
		}

		//print data to file
		for (unsigned i = 0; i < (*it_h)->get_V().size(); ++i)
		{
			fout << time << "\t" << (*it_h)->get_V()[i]->get_ID() << "\t";
			fout << (*it_h)->get_V()[i]->get_sequence()->get_sequence() << "\t";
			fout << (*it_h)->get_V()[i]->get_sequence()->get_fitness() << endl;
		}
		//close connection
		fout.close();
	}
}

void epidemics::change_fitness(mt19937 & gen)
{
	for (unsigned i = 0; i < hosts.size(); ++i)
	{
		for (unsigned j = 0; j < hosts[i]->get_V().size(); ++j)
		{
			if (hosts[i]->get_V()[j]->get_sequence()->get_fitness() > fit_low_cap)
			{
				hosts[i]->get_V()[j]->change_fitness(rnorm(1, fit_change, fit_change / 3, gen).back());
			}
			if (hosts[i]->get_V()[j]->get_sequence()->get_fitness() < fit_low_cap)
			{
				hosts[i]->get_V()[j]->change_fitness(fit_low_cap - hosts[i]->get_V()[j]->get_sequence()->get_fitness());
			}
		}
	}
}

bool epidemics::is_it_over()
{
	if (hosts.size() == 0)
	{
		cout << "\n**************************************" << endl;
		cout << "EPIDEMICS IS OVER" << endl;
		return 1;
	}
	return 0;
}

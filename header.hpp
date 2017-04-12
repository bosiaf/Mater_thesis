#include<iomanip>
#include<iostream>
#include<random>
#include<vector>
#include<sstream>
#include<string>

using namespace std;

namespace personal
{
	//global variables
	//
	//Why not using an Hash table to store the viral sequences? Easy to see if new seq is new or not.
	extern unsigned MAX_TSTEP;
	extern int v0, h0, hc_renew, burst_size;
	extern unsigned seed;
	extern string path_to_tmat, path_output_dyn, path_output_seq, seq0;
	extern double d_hc, d_ic, d_v, k_inf, S_df, k_btw, k_mut, fit_snp, fit_change, k_fit;
	extern vector<unsigned> SNPs;
	extern vector<long int> weight_not_snp;
	extern vector<double> fit_not_snp;
	extern mt19937 gen;
	extern vector<mt19937*> gens;
	extern bool dic_fit_dep, dv_fit_dep, inf_fit_dep, ad_imm_sys;


	//classes
	
	//sequence with the actual string and its fitness
	class sequences
	{
		private:
			string sequence;
			double fitness;
			unsigned s_size;
		public:
			sequences(string s_i, double fit);

			string get_sequence();
			void set_sequence(string s);
			void set_sequence(char nt, unsigned ind);

			double get_fitness();
			void set_fitness(double fit);
			unsigned get_size();

	};

	//virion class, don't know right now if I will use it, but it is there for 
	//single virion representation. it just contains an istance of a sequence class
	class virion
	{
		private:
			sequences * pars;	
		public:
			virion(sequences * s);
	};
	
	//a row of the R Viral.Array: it contains the pointer to the sequence class, 
	//the number of virions and the number of infected cells. 
	//Double for storage reasons.
	//You can also change the fitness from there.
	class strain
	{
		private:
			//which strain?
			sequences * seq;
			//number of virions
			long int vir;
			//number of infected cells
			long int icell;
			//cells that are infected and await evolution
			long int tcell;
			//give it an ID for recognition
			unsigned ID;
			//time of creation
			unsigned time;
		public:

			strain(sequences * s, unsigned id, unsigned t);
			strain(sequences * s, long int v, long int ic, long int tc, unsigned id, unsigned t);

			sequences * get_sequence();

			void change_fitness(double diff);

			unsigned get_ID();

			long int get_vir();
			void set_vir(int diff);

			long int get_icell();
			void set_icell(int diff);

			long int get_tcell();
			void set_tcell(int diff);

			unsigned get_time();
	};

	//a host class with the healthy cell number, a personal MT rng,
	//and a complete vector of all the strains in the host 
	//(see class strain)
	class host
	{
		private:
			vector<strain *> V;
			unsigned ID;
			mt19937 gen_l;
			long int healthy_cells;
			long int l_time;
			unsigned n_str;

		public:
			static unsigned total;
			
			host(long int hc, strain * s, mt19937 gen_local);

			long int get_ltime();
			void set_ltime(int diff);

			void new_str();

			unsigned get_ID();

			mt19937 get_rng();

			bool wi_host_inf_death(double k, double d_infc, double d_vir, long int burst, mt19937 & gen, vector<vector<double> > tmat, vector<unsigned> SNPs_listn, double p_mut, unsigned t);

			void add_line(strain * s);
			void delete_strain(unsigned index);

			long int get_hc();
			void set_hc(long int diff);

			double get_new_fitness(unsigned position, vector<unsigned> SNPs_list, char old_nt, char new_nt, mt19937 & gen);

			void evolve(mt19937 & gen, vector<vector<double> > tmat, vector<unsigned> SNPs_list, double p_mut, unsigned time);

			vector<strain *> get_V();

			void print();
			void print_seq();	
	};

	//master class containing all the host in play, it is
	//the total system. It contains the master time as well as
	//the population parameters like S, the susceptible population.
	class epidemics
	{
		private:
			vector<host*> hosts;
			long int time;
			long int S;

		public:
			epidemics();
			epidemics(host * h, long int s);

			long int get_time();
			void set_time(long int t);

			unsigned get_S();
			void set_S(long int diff);

			long int next_inf_time(double k, mt19937 & gen);
			void new_host_infection(mt19937 & gen);

			void add_host(host * h);
			void delete_host(unsigned index);

			vector<host*> get_hosts();

			void print_epidemics();
			void print_epidemics(string path);
			void print_seq_epidemics();
			void print_seq_epidemics(string path);

			void change_fitness(mt19937 & gen);

			bool is_it_over();
	};


	//overloaded operators
	istream& operator >> (istream& fin, vector<double> & row);
	istream& operator >> (istream& fin, vector<string> & row);
	istream& operator >> (istream& fin, vector<vector<double> > & table);

	//global functions
	
	//function to fo from nucleotide to integer and vice versa
	int nt_to_i(char nt);	
	
	char i_to_nt(int nt);
	//rexp generates a random exponential realization using the random library
	vector<double> rexp(int n, double lam, mt19937 & genr);

	//rbinom generates a random binomial realization using the random library
	vector<long int> rbinom(int n, unsigned drows, double p, mt19937 & genr);

	//rnorm generates a random normal realization using the random library
	vector<double> rnorm(int n, double mean, double sd, mt19937 & genr);

	//runif_int generates a uniform integer random realization using random library
	vector<long int> runif_int(int n, unsigned from, unsigned to, mt19937 & genr);

	//rber generates a random bernoulli distribution using random library
	vector<bool> rber(int n, double p, mt19937 & genr);

	//smpl_probs samples an index using an uniform real distribution(0,1) and a probability vector
	long int smpl_weight(vector<long int> w, mt19937 & gen, uniform_real_distribution<> d);
	long int smpl_weight(vector<double> w, mt19937 & gen, uniform_real_distribution<> d);			

	//sample samples size elements from 0 to (n-1) without replacement.
	vector<long int> sample(long int n, unsigned size, mt19937 & gen, uniform_real_distribution<> d);
	vector<long int> sample_int_wo_repl(long int n, unsigned size, mt19937 & gen);

	//rtp (rate to probability) converts ODE rates into probabilities through
	//1 - exp(-rate), since it is tuned to 1 day
	double rtp(double rate);

	//Gives as return a pair <old_nt, new_nt> if given a transition prob.matrix,
	//a random number from a uniform distribution in 0-1, and the old nucleotide.
	char evo_nt(vector<vector<double> > tmat, double random_nr, char old_nt);

	//read_in (might overload it) reads the transition matrix in
	//from a file connection
	void read_in(string file, vector<vector <double> > &trans_mat, bool header);
	void read_pars(string file, unsigned & max_tstep, string & path_to_tmat, string & path_output_dyn, string & path_output_seq, string & seq, vector<unsigned> & SNPs, int & v0, int & h0, int & hc_ren, double & dhc, double & dic, int & b_size, double & dv, double & kinf, double & sdf, double & kbtw, double & kmut, double & fit_snp, vector<double> & fit_not_snp, vector<long int> & weight_not_snp, bool & dic_fit_dep, bool & dv_fit_dep, bool & inf_fit_dep, double & k_fit, bool & ad_imm_sys, double & fit_change, unsigned & seed);
	bool fileExists(const std::string& file);
}

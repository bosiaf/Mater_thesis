#ifndef PAR_CLASS_HPP
#define PAR_CLASS_HPP
#include<iostream>
#include<vector>
#include<string>

using namespace std;

namespace epi
{
  struct parameters
  {
    const string path_to_tmat, path_output_dyn, path_output_seq, seq;
    const unsigned max_tstep, vol, v0, h0, hc_ren, b_size, seed, nr_chunks;
    const vector<unsigned> SNPs, weight_not_snp;
    const vector<double> fit_not_snp;
    const double dhc, dic, dv, kinf, sdf, kbtw, kmut, fit_snp, k_fit,
                 fit_change, fit_low_cap, fit_high_cap;
  

    parameters(const string file, p_tmat, p_out_dyn, p_out_seq, s0,
               const unsigned MAX_TSTEP, VOL, V0, H0, HC_REN, B_SIZE, SEED, NR_CHUNKS,
               const vector<unsigned> SNPs, weight_not_snp,
               const vector<double> fit_not_snp,
               const double DHC, DIC, DV, KINF, SDF, KBTW, KMUT, FIT_SNP, K_FIT,
                            FIT_CHANGE, FIT_LOW_CAP, FIT_HIGH_CAP);
  };

}

#endif

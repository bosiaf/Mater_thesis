#ifndef PAR_CLASS_HPP
#define PAR_CLASS_HPP
#include<iostream>
#include<vector>
#include<string>

using namespace std;

namespace epi
{
  struct par
  {
    const string path_to_tmat, path_output_dyn, path_output_seq, seq;
    const unsigned max_tstep, vol, v0, h0, hc_ren, b_size, seed, nr_chunks;
    const vector<unsigned> SNPs, weight_not_snp;
    const vector<double> fit_not_snp;
    const double dhc, dic, dv, kinf, sdf, kbtw, kmut, fit_snp, k_fit,
                 fit_change, fit_low_cap, fit_high_cap;
  
    //copy constructor by compiler
    parameters(const string p_tmat, p_out_dyn, p_out_seq, s0,
               const unsigned MAX_TSTEP, VOL, V0, H0, HC_REN, B_SIZE, SEED, NR_CHUNKS,
               const vector<unsigned> SNPs, weight_not_snp,
               const vector<double> fit_not_snp,
               const double DHC, DIC, DV, KINF, SDF, KBTW, KMUT, FIT_SNP, K_FIT,
                            FIT_CHANGE, FIT_LOW_CAP, FIT_HIGH_CAP)
              : path_to_tmat(p_tmat),
                path_output_dyn(p_out_dyn),
                path_output_seq(p_out_seq),
                seq(s0),
                max_tstep(MAX_TSTEP),
                vol(VOL), v0(V0), h0(H0), hc_ren(HC_REN),
                b_size(B_SIZE), seed(SEED), nr_chunks(NR_CHUNKS),
                SNPs(SNPs), weight_not_snp(weight_not_snp),
                fit_not_snp(fit_not_snp), 
                dhc(DHC), dic(DIC), dv(DV), kinf(KINF), sdf(SDF),
                kbtw(KBTW), kmut(KMUT), fit_snp(FIT_SNP), kfit(KFIT),
                fit_change(FIT_CHANGE), fit_low_cap(FIT_LOW_CAP),
                fit_high_cap(FIT_HIGH_CAP)
    { }

  };
  
  par read_pars(const string file);
}

#endif

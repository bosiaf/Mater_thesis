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
    const unsigned max_tstep, v0, h0, hc_ren, b_size, seed, nr_chunks;
    const vector<unsigned> SNPs, weight_not_snp;
    const vector<double> fit_not_snp;
    const double vol, dhc, dic, dv, kinf, sdf, dl, inf_to_lat, lat_act, 
                 kbtw, kmut, fit_snp, k_fit,
                 fit_change, fit_low_cap, fit_high_cap;
    const bool dic_fit_dep, dv_fit_dep, inf_fit_dep, ad_imm_sys, parallel,
               seq_per_time, seq_print;
    void print_par() const;
    //copy constructor by compiler
    par(const string p_out_dyn, const string p_out_seq, const string s0, 
               const unsigned max_tstep, const unsigned v0, 
               const unsigned h0, const unsigned hc_ren, const unsigned b_size, 
               const unsigned seed, const unsigned nr_chunks,
               const vector<unsigned> SNPs, const vector<unsigned> weight_not_snp,
               const vector<double> fit_not_snp,
               const double vol, const double dhc, const double dic, const double dv,
               const double dl, const double inf_to_lat, const double lat_act, 
               const double kinf, const double sdf, const double kbtw, 
               const double kmut, const double fit_snp, const double k_fit,
               const double fit_change, const double fit_low_cap,
               const double fit_high_cap, const bool dic_fit_dep,
               const bool dv_fit_dep, const bool inf_fit_dep, 
               const bool ad_imm_sys, const bool parallel, const bool seq_per_time,
               const bool seq_print)
              : path_output_dyn(p_out_dyn),
                path_output_seq(p_out_seq),
                seq(s0),
                max_tstep(max_tstep),
                v0(v0), h0(h0), hc_ren(hc_ren),
                b_size(b_size), seed(seed), nr_chunks(nr_chunks),
                SNPs(SNPs), weight_not_snp(weight_not_snp),
                fit_not_snp(fit_not_snp), vol(vol), 
                dhc(dhc), dic(dic), dv(dv), kinf(kinf), sdf(sdf),
                dl(dl), inf_to_lat(inf_to_lat), lat_act(lat_act),
                kbtw(kbtw), kmut(kmut), fit_snp(fit_snp), k_fit(k_fit),
                fit_change(fit_change), fit_low_cap(fit_low_cap),
                fit_high_cap(fit_high_cap),
                dic_fit_dep(dic_fit_dep), dv_fit_dep(dv_fit_dep),
                inf_fit_dep(inf_fit_dep), ad_imm_sys(ad_imm_sys),
                parallel(parallel), seq_per_time(seq_per_time),
                seq_print(seq_print)
    { }

  };
  
  const par read_pars(const string file);
}

#endif

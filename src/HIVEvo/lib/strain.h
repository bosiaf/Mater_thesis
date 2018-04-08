#ifndef STRAIN_HPP
#define STRAIN_HPP
#include<iostream>
#include<string>
#include <utility>
using namespace std;

namespace epi {
class strain {
 public:
  typedef unsigned count_t;

  //Public members variables
  count_t vir; // number of virions
  count_t inf_cell; // number of infected cells
  count_t temp_cell; // temporarly infected cells for evolution
  count_t lat_cell;//latently infected cells
  double fitness; //fitness of the strain

  //getters
  string get_sequence() const { return sequence; }
  double get_fitness() const { return fitness; }
  static unsigned get_size() { return s_size; }
  count_t get_ID() const { return ID; }
  unsigned get_time() const { return time; }
  unsigned get_hash() const { return hash; }

  static void set_size(unsigned newSize) { s_size = newSize; };

  strain(const count_t v, const count_t i_c, const count_t t_c,
         const count_t l_c,
         const double fit, const string& s, const unsigned id,
         const unsigned t, const unsigned hash)
    : vir(v), inf_cell(i_c), temp_cell(t_c),
      lat_cell(l_c), fitness(fit),
      sequence(s), ID(id), time(t), hash(hash) {}

 private:
  string sequence; //sequence of the strain
  static unsigned s_size; //sequence size
  count_t ID; //ID of the strain
  unsigned time; //time of creation
  //hash for sequence search speedup
  unsigned hash;
};
}
#endif

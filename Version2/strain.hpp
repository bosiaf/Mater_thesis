#ifndef STRAIN_HPP
#define STRAIN_HPP
#include<iostream>
#include<string>

using namespace std;

namespace epi
{
  class strain
  {
    public:
      typedef long count_t; 

      //Public members variables
      count_t vir; // number of virions
      count_t inf_cell; // number of infected cells
      count_t temp_cell; // temporarly infected cells for evolution
      count_t lat_cell;//latently infected cells
      double fitness; //fitness of the strain
      
      //getters
      string get_sequence() const
      { return sequence; }
      double get_fitness() const
      { return fitness; }
      static const unsigned get_size()
      { return s_size; }
      const count_t get_ID() const
      { return ID; }
      const unsigned get_time() const
      { return time; }
      const unsigned get_hash() const
      { return hash; } 

      strain (const count_t v, const count_t i_c, const count_t t_c, 
              const count_t l_c, 
              const double fit, const string & s, const unsigned id,
              const unsigned t, const unsigned hash)
      : vir(v), inf_cell(i_c), temp_cell(t_c),
        lat_cell(l_c), fitness(fit),
        sequence(s), ID(id), time(t), hash(hash)
      {}
      
    private:
      const string sequence; //sequence of the strain
      static unsigned s_size; //sequence size
      const count_t ID; //ID of the strain
      const unsigned time; //time of creation
      //hash for sequence search speedup
      const unsigned hash;
  };    
}

#endif

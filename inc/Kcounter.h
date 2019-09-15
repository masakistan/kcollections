#pragma once

#include <pybind11/pybind11.h>
#include "Kcontainer.h"

namespace py = pybind11;

class Kcounter
{
private:
  Kcontainer* kc;
  int m_k;
  std::function<int(int, int)> merge_func = [] (int prev_val, int new_val)->int{ return prev_val + new_val; };
public:
  Kcounter( const int k );
  ~Kcounter();
  void insert( char* kmer, count_dtype count );
  bool contains( char* kmer );
  void clear();
  uint64_t size();
  void remove( char* kmer );
  void add_seq(const char* seq, uint32_t length);
  count_dtype get( char* kmer );
  int get_k() { return m_k; }
  Kcontainer* get_kc() { return kc; }
  void parallel_add_init(int threads) {
    parallel_kcontainer_add_init(kc, threads, merge_func);
  };
  void parallel_add(const char* kmer) {
    parallel_kcontainer_add(kc, kmer);
  }
  void parallel_add_join() {
    parallel_kcontainer_add_join(kc);
  }

  void parallel_add_seq(const char* seq, uint32_t length);

  std::string get_uc_kmer( Vertex<int>* v, int k, int idx )
  {
    UC<int>* uc = v->get_uc();
    char* kmer = deserialize_kmer(k, uc->get_suffix(k, idx));
    std::string skmer(kmer);
    free(kmer);
    return skmer;
  }

  int get_uc_size( Vertex<int>* v )
  {
    return v->get_uc()->get_size();
  }

  Vertex<int>* get_root() { return kc->v; }
  int get_vs_size( Vertex<int>* v ){ return v->get_vs_size(); }
  Vertex<int>* get_child_vertex( Vertex<int>* v, int idx )
  {
    return v->get_vs()[idx];
  }
  std::string get_child_suffix( Vertex<int>* v, int idx )
  {
    return kcontainer_get_child_suffix(v, idx);
  }
};

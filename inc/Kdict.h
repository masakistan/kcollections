#pragma once

#include <pybind11/pybind11.h>
#include "Kcontainer.h"

namespace py = pybind11;

class Kdict
{
private:
  Kcontainer* kc;
  int m_k;
  std::function<int(int, int)> merge_func = [] (int prev_val, int new_val)->int{ return new_val; };
public:
  Kdict( const int k );
  ~Kdict();
  void add( char* kmer, int obj );
  bool contains( char* kmer );
  void clear();
  uint64_t size();
  void remove( char* kmer );
  int get( char* kmer );
  int get_k() { return m_k; }
  Kcontainer* get_kc() { return kc; }
  void add_seq(const char* seq, uint32_t length, py::iterable values, std::function<int(int, int)> &f);

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
  
  int get_vs_size( Vertex<int>* v ){ return v->get_vs_size(); }
  
  Vertex<int>* get_child_vertex( Vertex<int>* v, int idx )
  {
    return v->get_vs()[idx];
  }
  
  std::string get_child_suffix( Vertex<int>* v, int idx )
  {
    return kcontainer_get_child_suffix(v, idx);
  }
  
  Vertex<int>* get_root() { return kc->v; }

  void parallel_add_init(int threads, const std::function<int(int, int)> &f)  {
    parallel_kcontainer_add_init(kc, threads, f);
  }
  
  void parallel_add(const char* kmer, int value) {
    parallel_kcontainer_add(kc, kmer, value);
  }
  
  void parallel_add_join() {
    parallel_kcontainer_add_join(kc);
  }

  void parallel_add_seq(const char* seq, uint32_t length, py::iterable values) {
    parallel_kcontainer_add_seq(kc, seq, length, values);
  }
};

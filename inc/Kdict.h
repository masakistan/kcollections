#pragma once

#include <pybind11/pybind11.h>
#include "Kcontainer.h"

namespace py = pybind11;

class Kdict
{
private:
  Kcontainer* kc;
  int m_k;
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
  void add_seq( char* seq, uint32_t length, py::iterable values, std::function<int(int, int)> &f);

  std::string get_uc_kmer( Vertex* v, int k, int idx )
  {
    UC<int>* uc = v->uc;
    char* kmer = deserialize_kmer(k, uc->get_suffix(k, idx));
    std::string skmer(kmer);
    free(kmer);
    return skmer;
  }

  int get_uc_size( Vertex* v )
  {
    return v->uc->get_size();
  }
  
  int get_vs_size( Vertex* v ){ return v->vs_size; }
  
  Vertex* get_child_vertex( Vertex* v, int idx )
  {
    return &v->vs[idx];
  }
  
  std::string get_child_suffix( Vertex* v, int idx )
  {
    return kcontainer_get_child_suffix(v, idx);
  }
  
  Vertex* get_root() { return &kc->v; }

  void parallel_add_init(int threads, const std::function<int(int, int)> &f)  {
    parallel_kcontainer_add_init(kc, threads, f);
  }
  
  void parallel_add(const char* kmer, int value) {
    parallel_kcontainer_add(kc, kmer, value);
  }
  
  void parallel_add_join() {
    parallel_kcontainer_add_join(kc);
  }

  void parallel_add_seq(char* seq, uint32_t length, py::iterable values) {
    parallel_kcontainer_add_seq(kc, seq, length, values);
  }
};

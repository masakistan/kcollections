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
  void add( char* kmer, py::handle* obj );
  bool contains( char* kmer );
  void clear();
  uint64_t size();
  void remove( char* kmer );
  py::handle* get( char* kmer );
  int get_k() { return m_k; }
  Kcontainer* get_kc() { return kc; }

  char* get_uc_kmer( Vertex* v, int k, int idx )
  {
    UC* uc = &v->uc;
    int bk = calc_bk( k );
    int suffix_idx = bk * idx;
    return deserialize_kmer( k, bk, &uc->suffixes[ suffix_idx ] );
  }

  int get_uc_size( Vertex* v )
  {
    return v->uc.size;
  }
  
  int get_vs_size( Vertex* v ){ return v->vs_size; }
  
  Vertex* get_child_vertex( Vertex* v, int idx )
  {
    return &v->vs[idx];
  }
  
  char* get_child_suffix( Vertex* v, int idx )
  {
    return kcontainer_get_child_suffix(v, idx);
  }
  
  Vertex* get_root() { return &kc->v; }

  void parallel_add_init(int threads) {
    parallel_kcontainer_add_init(kc, threads);
  }
  
  void parallel_add(const char* kmer, py::handle* value) {
    parallel_kcontainer_add(kc, kmer, value);
  }
  
  void parallel_add_join() {
    parallel_kcontainer_add_join(kc);
  }

  void parallel_add_seq(char* seq, uint32_t length, py::iterable* values) {
    parallel_kcontainer_add_seq(kc, seq, length, values);
  }
};

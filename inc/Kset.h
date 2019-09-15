#pragma once

#include "uint256_t.h"
#include <pybind11/pybind11.h>
#include "Kcontainer.h"

namespace py = pybind11;

class Kset
{
private:
  Kcontainer* kc;
  int m_k;
public:
  Kset( const int k );
  ~Kset();
  void add( const char* kmer );
  bool contains( const char* kmer );
  void clear();
  uint64_t size();
  void remove( const char* kmer );
  int get_k() { return m_k; }
  Kcontainer* get_kc() { return kc; }

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
  void add_seq(const char* seq, uint32_t length);

  void parallel_add_init(int threads) {
    parallel_kcontainer_add_init(kc, threads);
  };
  void parallel_add(const char* kmer) {
    parallel_kcontainer_add(kc, kmer);
  }
  void parallel_add_join() {
    parallel_kcontainer_add_join(kc);
  }

  void parallel_add_seq(const char* seq, uint32_t length);
};

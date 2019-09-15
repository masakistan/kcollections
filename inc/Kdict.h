#pragma once

#include <pybind11/pybind11.h>
#include "Kcontainer.h"

namespace py = pybind11;

template <class T>
class Kdict
{
private:
  Kcontainer<T>* kc;
  int m_k;
  std::function<T(T, T)> merge_func = [] (T prev_val, T new_val)->T{ return new_val; };
public:
  Kdict(const int k) {
    kc = new Kcontainer<T>(k);
    m_k = k;
  }
  
  ~Kdict() {
    delete kc;
  }
  
  void add(char* kmer, T obj) {
    CHECK_KMER_LENGTH(kmer, m_k, "Kdict");
    kc->kcontainer_add(kmer, obj, merge_func);
  }
  
  bool contains(char* kmer) {
    CHECK_KMER_LENGTH(kmer, m_k, "Kdict");
    return kc->kcontainer_contains(kmer);
  }
  
  void clear() {
    delete kc;
    kc = new Kcontainer<T>(m_k);
  }

  uint64_t size() {
    return kc->kcontainer_size();
  }
  
  void remove(char* kmer) {
    CHECK_KMER_LENGTH(kmer, m_k, "Kdict");
    kc->kcontainer_remove(kmer);
  }

  T get(char* kmer) {
    CHECK_KMER_LENGTH(kmer, m_k, "Kdict");
    return kc->kcontainer_get(kmer);
  }
  
  int get_k() { return m_k; }
  Kcontainer<T>* get_kc() { return kc; }
  
  void add_seq(const char* seq, uint32_t length, py::iterable values, std::function<T(T, T)> &f) {
    kc->kcontainer_add_seq(seq, length, values, f);
  }

  std::string get_uc_kmer( Vertex<T>* v, int k, int idx )
  {
    UC<T>* uc = v->get_uc();
    char* kmer = deserialize_kmer(k, uc->get_suffix(k, idx));
    std::string skmer(kmer);
    free(kmer);
    return skmer;
  }

  int get_uc_size( Vertex<T>* v )
  {
    return v->get_uc()->get_size();
  }
  
  int get_vs_size( Vertex<T>* v ){ return v->get_vs_size(); }
  
  Vertex<T>* get_child_vertex( Vertex<T>* v, int idx )
  {
    return v->get_vs()[idx];
  }
  
  std::string get_child_suffix( Vertex<T>* v, int idx )
  {
    return kc->kcontainer_get_child_suffix(v, idx);
  }
  
  Vertex<T>* get_root() { return kc->get_v(); }

  void parallel_add_init(int threads, const std::function<T(T, T)> &f)  {
    kc->parallel_kcontainer_add_init(threads, f);
  }
  
  void parallel_add(const char* kmer, T value) {
    kc->parallel_kcontainer_add(kmer, value);
  }
  
  void parallel_add_join() {
    kc->parallel_kcontainer_add_join();
  }

  void parallel_add_seq(const char* seq, uint32_t length, py::iterable values) {
    kc->parallel_kcontainer_add_seq(seq, length, values);
  }
};
template<class T>
ThreadGlobals<T>* Kcontainer<T>::tg = (ThreadGlobals<T>*) calloc(1, sizeof(ThreadGlobals<T>));

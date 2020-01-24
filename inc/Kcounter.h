#pragma once

#include <functional>
#include "Kcontainer.h"

class Kcounter
{
private:
  Kcontainer<int>* kc;
  int m_k;
  std::function<int(int&, int&)> merge_func = [] (int& prev_val, int& new_val)->int&{
						prev_val += new_val;
						return prev_val;
					       };
  std::function<int(int&, int&)>overwrite_merge_func = [] (int& prev_val, int& new_val)->int&{ return new_val;};
public:
  Kcounter( const int k );
  ~Kcounter();
  void insert( char* kmer, count_dtype count );
  bool contains( char* kmer );
  void clear();
  uint64_t size();
  void remove( char* kmer );
  void add_seq(const char* seq);
  count_dtype get( char* kmer );
  int get_k() { return m_k; }
  Kcontainer<int>* get_kc() { return kc; }
  void parallel_add_init(int threads) {
    kc->parallel_kcontainer_add_init(threads, merge_func);
  };
  void parallel_add(const char* kmer) {
    kc->parallel_kcontainer_add(kmer);
  }
  void parallel_add_join() {
    kc->parallel_kcontainer_add_join();
  }

  void parallel_add_seq(const char* seq);

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

  Vertex<int>* get_root() { return kc->get_v(); }
  int get_vs_size( Vertex<int>* v ){ return v->get_vs_size(); }
  Vertex<int>* get_child_vertex( Vertex<int>* v, int idx )
  {
    return &v->get_vs()[idx];
  }
  std::string get_child_suffix( Vertex<int>* v, int idx )
  {
    return kc->kcontainer_get_child_suffix(v, idx);
  }

  Kcontainer<int>::iterator begin() {
    return kc->begin();
  }

  Kcontainer<int>::iterator end() {
    return kc->end();
  }
};
template<class T>
ThreadGlobals<T>* Kcontainer<T>::tg = (ThreadGlobals<T>*) calloc(1, sizeof(ThreadGlobals<T>));

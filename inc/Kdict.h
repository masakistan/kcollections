#pragma once

#include <fstream>
#include <functional>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include "Kcontainer.h"

template <class T>
class Kdict
{
private:
  Kcontainer<T>* kc;
  int m_k;
  std::function<T(T&, T&)> merge_func;
  std::function<T(T&, T&)> overwrite_merge_func;
public:
  Kdict(){
    merge_func = [] (T& prev_val, T& new_val)->T&{ return new_val;};
    overwrite_merge_func = [] (T& prev_val, T& new_val)->T&{ return new_val;};
  }
  Kdict(const int k) : m_k(k) {
    kc = new Kcontainer<T>(k);
    merge_func = [] (T& prev_val, T& new_val)->T&{ return new_val;};
    overwrite_merge_func = [] (T& prev_val, T& new_val)->T&{ return new_val;};
  }
  
  ~Kdict() {
    delete kc;
  }

  template<class Archive>
  void save(Archive& ar, const unsigned int version) const {
    ar & m_k;
    ar & *kc;
  }

  template<class Archive>
  void load(Archive& ar, const unsigned int version) {
    ar & m_k;
    CDEPTH = calc_bk(m_k);

    kc = new Kcontainer<T>(m_k);
    ar & *kc;
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()

  void write(const char* opath) {
    CDEPTH = calc_bk(get_k());
    std::ofstream ofs(opath);

    boost::archive::binary_oarchive oa(ofs);
    oa << *this;
    CDEPTH = -1;
  }

  void read(const char* ipath) {
    std::ifstream ifs(ipath);

    boost::archive::binary_iarchive ia(ifs);
    ia >> *this;

    CDEPTH = -1;
  }

  
  void set_merge_func(std::function<T(T&, T&)> merge_func) {
    this->merge_func = merge_func;
  }     
  
  void add(char* kmer, T& obj) {
    CHECK_KMER_LENGTH(kmer, m_k, "Kdict");
    kc->kcontainer_add(kmer, obj, overwrite_merge_func);
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
  
#if defined(PYTHON)
  void add_seq(const char* seq, py::iterable& values)
#else
    template <typename Iter>
  void add_seq(const char* seq, Iter& values)
#endif
  {
    kc->kcontainer_add_seq(seq, strlen(seq), values, merge_func);
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
    return &v->get_vs()[idx];
  }
  
  std::string get_child_suffix( Vertex<T>* v, int idx )
  {
    return kc->kcontainer_get_child_suffix(v, idx);
  }
  
  Vertex<T>* get_root() { return kc->get_v(); }

  void parallel_add_init(int threads)  {
    kc->parallel_kcontainer_add_init(threads, merge_func);
  }
  
  void parallel_add(const char* kmer, T& value) {
    kc->parallel_kcontainer_add(kmer, value);
  }
  
  void parallel_add_join() {
    kc->parallel_kcontainer_add_join();
  }

#if defined(PYTHON)
  void parallel_add_seq(const char* seq, py::iterable& values)
#else
    template <typename Iter>
  void parallel_add_seq(const char* seq, Iter& values)
#endif
  {
    //std::cout << "parallel adding seq" << std::endl;
    kc->parallel_kcontainer_add_seq(seq, strlen(seq), values);
  }

  typedef typename Kcontainer<T>::iterator iterator;
  iterator begin() {
    return kc->begin();
  }

  iterator end() {
    return kc->end();
  }
};


#pragma once

#include <fstream>
#include <functional>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

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
  Kcounter(){};
  Kcounter( const int k );
  ~Kcounter();

  template<class Archive>
  void save(Archive& ar, const unsigned int version) const {
    ar & m_k;
    ar & *kc;
  }

  template<class Archive>
  void load(Archive& ar, const unsigned int version) {
    ar & m_k;
    CDEPTH = calc_bk(m_k);
    kc = new Kcontainer<int>(m_k);
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
  
  void insert( const char* kmer, count_dtype count );
  bool contains( const char* kmer );
  void clear();
  uint64_t size();
  void remove( const char* kmer );
  void add_seq(const char* seq);
  count_dtype get( const char* kmer );
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


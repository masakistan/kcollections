#pragma once

#include <fstream>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include "uint256_t.h"
#include "Kcontainer.h"
#include "globals.h"
#include "helper.h"

class Kset
{
private:
  Kcontainer* kc;
  int m_k;
public:
  Kset(){};
  Kset( const int k );
  ~Kset();

  template<class Archive>
  void save(Archive& ar, const unsigned int version) const {
    ar & m_k;
    ar & *kc;
  }

  template<class Archive>
  void load(Archive& ar, const unsigned int version) {
    ar & m_k;
    CDEPTH = calc_bk(m_k);
    kc = new Kcontainer(m_k);
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
  
  void add( const char* kmer );
  bool contains( const char* kmer );
  void clear();
  uint64_t size();
  void remove( const char* kmer );
  int get_k() { return m_k; }
  Kcontainer* get_kc() { return kc; }

  std::string get_uc_kmer( Vertex* v, int k, int idx )
  {
    UC* uc = v->get_uc();
    char* kmer = deserialize_kmer(k, uc->get_suffix(k, idx));
    std::string skmer(kmer);
    free(kmer);
    return skmer;
}

  int get_uc_size( Vertex* v )
  {
    return v->get_uc()->get_size();
  }

  Vertex* get_root() { return kc->get_v(); }
  int get_vs_size( Vertex* v ){ return v->get_vs_size(); }
  Vertex* get_child_vertex( Vertex* v, int idx )
  {
    return &v->get_vs()[idx];
  }

  std::string get_child_suffix(Vertex* v, int idx)
  {
    return kc->kcontainer_get_child_suffix(v, idx);
  }
  void add_seq(const char* seq);

  void parallel_add_init(int threads) {
    kc->parallel_kcontainer_add_init(threads);
  }
  
  void parallel_add(const char* kmer) {
    kc->parallel_kcontainer_add(kmer);
  }
  void parallel_add_join() {
    kc->parallel_kcontainer_add_join();
  }

  void parallel_add_seq(const char* seq);

  Kcontainer::iterator begin() {
    return kc->begin();
  }

  Kcontainer::iterator end() {
    return kc->end();
  }
};


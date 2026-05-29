#pragma once

#include <fstream>
#include <functional>

#include "kc/kcounter_core.h"

class Kcounter
{
private:
  KcCounterContainer<int>* kc;
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
    kc->save(ar, version);
  }

  template<class Archive>
  void load(Archive& ar, const unsigned int version) {
    ar & m_k;
    CDEPTH = calc_bk(m_k);
    kc = new KcCounterContainer<int>(m_k);
    kc->load(ar, version);
  }

  void write(const char* opath) {
    CDEPTH = calc_bk(get_k());
    std::ofstream ofs(opath, std::ios::binary);
    kc_io::BinaryOutArchive ar(ofs);
    kc_io::write_file_header(ar, kc_io::ContainerKind::KIND_COUNTER);
    save(ar, 0);
    CDEPTH = -1;
  }

  void read(const char* ipath) {
    std::ifstream ifs(ipath, std::ios::binary);
    kc_io::BinaryInArchive ar(ifs);
    kc_io::read_file_header(ar, kc_io::ContainerKind::KIND_COUNTER);
    load(ar, 0);
    CDEPTH = -1;
  }
  
  void insert( const char* kmer, count_dtype count );
  bool contains( const char* kmer );
  void clear();
  uint64_t size();
  void remove( const char* kmer );
  void add_seq(const char* seq);
  void add_seq(const char* seq, size_t length);
  count_dtype get( const char* kmer );
  int get_k() { return m_k; }
  KcCounterContainer<int>* get_kc() { return kc; }
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

  std::string get_uc_kmer( KcCounterVertex<int>* v, int k, int idx )
  {
    KcCounterUC<int>* uc = v->get_uc();
    char* kmer = deserialize_kmer(k, uc->get_suffix(k, idx));
    std::string skmer(kmer);
    free(kmer);
    return skmer;
  }

  int get_uc_size( KcCounterVertex<int>* v )
  {
    return v->get_uc()->get_size();
  }

  KcCounterVertex<int>* get_root() { return kc->get_v(); }
  int get_vs_size( KcCounterVertex<int>* v ){ return v->get_vs_size(); }
  KcCounterVertex<int>* get_child_vertex( KcCounterVertex<int>* v, int idx )
  {
    return &v->get_vs()[idx];
  }
  std::string get_child_suffix( KcCounterVertex<int>* v, int idx )
  {
    return kc->kcontainer_get_child_suffix(v, idx);
  }

  KcCounterContainer<int>::iterator begin() {
    return kc->begin();
  }

  KcCounterContainer<int>::iterator end() {
    return kc->end();
  }
};


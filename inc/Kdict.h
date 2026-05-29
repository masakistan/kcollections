#pragma once

#include <fstream>
#include <functional>

#include "kc/kdict_core.h"
#if defined(PYTHON)
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif
#include "opaque.h"

template <class T>
class Kdict
{
private:
  KcDictContainer<T>* kc;
  int m_k;
  std::function<T(T&, T&)> merge_func;
  std::function<T(T&, T&)> overwrite_merge_func;
public:
  Kdict(){
    merge_func = [] (T& prev_val, T& new_val)->T&{ return new_val;};
    overwrite_merge_func = [] (T& prev_val, T& new_val)->T&{ return new_val;};
  }
  Kdict(const int k) : m_k(k) {
    kc = new KcDictContainer<T>(k);
    merge_func = [] (T& prev_val, T& new_val)->T&{ return new_val;};
    overwrite_merge_func = [] (T& prev_val, T& new_val)->T&{ return new_val;};
  }
  
  ~Kdict() {
    delete kc;
  }

  template<class Archive>
  void save(Archive& ar, const unsigned int version) const {
    ar & m_k;
    kc->save(ar, version);
  }

  template<class Archive>
  void load(Archive& ar, const unsigned int version) {
    ar & m_k;
    CDEPTH = calc_bk(m_k);

    kc = new KcDictContainer<T>(m_k);
    kc->load(ar, version);
  }

  void write(const char* opath) {
    CDEPTH = calc_bk(get_k());
    std::ofstream ofs(opath, std::ios::binary);
    kc_io::BinaryOutArchive ar(ofs);
    kc_io::write_file_header(ar, kc_io::ContainerKind::KIND_DICT);
    save(ar, 0);
    CDEPTH = -1;
  }

  void read(const char* ipath) {
    std::ifstream ifs(ipath, std::ios::binary);
    const kc_io::FileHeader hdr =
        kc_io::read_file_header_raw(ifs, kc_io::ContainerKind::KIND_DICT);
    kc_io::require_little_endian_for_v1(hdr.version);
    if (hdr.version == kc_io::FILE_VERSION) {
      kc_io::BinaryInArchive ar(ifs);
      load(ar, 0);
    } else {
      kc_io::BinaryInArchiveNative ar(ifs);
      load(ar, 0);
    }
    CDEPTH = -1;
  }

  
  void set_merge_func(std::function<T(T&, T&)> merge_func) {
    this->merge_func = merge_func;
  }     
  
  void add(const char* kmer, T obj) {
    CHECK_KMER_LENGTH(kmer, m_k, "Kdict");
    kc->kcontainer_add(kmer, obj, overwrite_merge_func);
  }
  
  bool contains(const char* kmer) {
    CHECK_KMER_LENGTH(kmer, m_k, "Kdict");
    return kc->kcontainer_contains(kmer);
  }
  
  void clear() {
    delete kc;
    kc = new KcDictContainer<T>(m_k);
  }

  uint64_t size() {
    return kc->kcontainer_size();
  }
  
  void remove(const char* kmer) {
    CHECK_KMER_LENGTH(kmer, m_k, "Kdict");
    kc->kcontainer_remove(kmer);
  }

  T& get(const char* kmer) {
    CHECK_KMER_LENGTH(kmer, m_k, "Kdict");
    return kc->kcontainer_get(kmer);
  }
  
  int get_k() { return m_k; }
  KcDictContainer<T>* get_kc() { return kc; }
  
#if defined(PYTHON)
  void add_seq(const char* seq, py::iterable& values)
#else
    template <typename Iter>
  void add_seq(const char* seq, Iter& values)
#endif
  {
    kc->kcontainer_add_seq(seq, strlen(seq), values, merge_func);
  }

  std::string get_uc_kmer( KcDictVertex<T>* v, int k, int idx )
  {
    KcDictUC<T>* uc = v->get_uc();
    char* kmer = deserialize_kmer(k, uc->get_suffix(k, idx));
    std::string skmer(kmer);
    free(kmer);
    return skmer;
  }

  int get_uc_size( KcDictVertex<T>* v )
  {
    return v->get_uc()->get_size();
  }
  
  int get_vs_size( KcDictVertex<T>* v ){ return v->get_vs_size(); }
  
  KcDictVertex<T>* get_child_vertex( KcDictVertex<T>* v, int idx )
  {
    return &v->get_vs()[idx];
  }
  
  std::string get_child_suffix( KcDictVertex<T>* v, int idx )
  {
    return kc->kcontainer_get_child_suffix(v, idx);
  }
  
  KcDictVertex<T>* get_root() { return kc->get_v(); }

  void parallel_add_init(int threads)  {
    kc->parallel_kcontainer_add_init(threads, merge_func);
  }
  
  void parallel_add(const char* kmer, T value) {
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

  typedef typename KcDictContainer<T>::iterator iterator;
  iterator begin() {
    return kc->begin();
  }

  iterator end() {
    return kc->end();
  }
};


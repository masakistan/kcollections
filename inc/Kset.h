#pragma once

#include <fstream>

#include "kc/kset_core.h"
#include "helper.h"

class Kset {
 private:
  KcSetContainer* kc;
  int m_k;

 public:
  Kset() {}
  Kset(const int k);
  ~Kset();

  template <class Archive>
  void save(Archive& ar, const unsigned int version) const {
    ar& m_k;
    kc->save(ar, version);
  }

  template <class Archive>
  void load(Archive& ar, const unsigned int version) {
    ar& m_k;
    CDEPTH = calc_bk(m_k);
    kc = new KcSetContainer(m_k);
    kc->load(ar, version);
  }

  void write(const char* opath) {
    CDEPTH = calc_bk(get_k());
    std::ofstream ofs(opath, std::ios::binary);
    kc_io::BinaryOutArchive ar(ofs);
    kc_io::write_file_header(ar, kc_io::ContainerKind::KIND_SET);
    save(ar, 0);
    CDEPTH = -1;
  }

  void read(const char* ipath) {
    std::ifstream ifs(ipath, std::ios::binary);
    kc_io::BinaryInArchive ar(ifs);
    kc_io::read_file_header(ar, kc_io::ContainerKind::KIND_SET);
    load(ar, 0);
    CDEPTH = -1;
  }

  void add(const char* kmer);
  bool contains(const char* kmer);
  void clear();
  uint64_t size();
  void remove(const char* kmer);
  int get_k() { return m_k; }
  KcSetContainer* get_kc() { return kc; }

  std::string get_uc_kmer(KcSetVertex* v, int k, int idx) {
    KcSetUC* uc = v->get_uc();
    char* kmer = deserialize_kmer(k, uc->get_suffix(k, idx));
    std::string skmer(kmer);
    free(kmer);
    return skmer;
  }

  int get_uc_size(KcSetVertex* v) { return v->get_uc()->get_size(); }

  KcSetVertex* get_root() { return kc->get_v(); }
  int get_vs_size(KcSetVertex* v) { return v->get_vs_size(); }
  KcSetVertex* get_child_vertex(KcSetVertex* v, int idx) { return &v->get_vs()[idx]; }

  std::string get_child_suffix(KcSetVertex* v, int idx) {
    return kc->kcontainer_get_child_suffix(v, idx);
  }
  void add_seq(const char* seq);
  void add_seq(const char* seq, size_t length);

  void parallel_add_init(int threads) { kc->parallel_kcontainer_add_init(threads); }

  void parallel_add(const char* kmer) { kc->parallel_kcontainer_add(kmer); }
  void parallel_add_join() { kc->parallel_kcontainer_add_join(); }

  void parallel_add_seq(const char* seq);

  KcSetContainer::iterator begin() { return kc->begin(); }

  KcSetContainer::iterator end() { return kc->end(); }
};

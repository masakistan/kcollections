#pragma once

#if defined(PYTHON)
#include <pybind11/pybind11.h>
#endif

#include <functional>

#include "globals.h"
#include "kc_io.h"

#include "UContainer.h"
//#include <jemalloc/jemalloc.h>
#include "pref_mask.h"
#include "kc/class_names.h"

#if defined(PYTHON)
namespace py = pybind11;
#endif

#if defined(KDICT) || defined(KCOUNTER)
template <class T>
#endif
class KC_VERTEX {
private:
#if defined(KDICT) || defined(KCOUNTER)
  KC_VERTEX<T>* vs;
#else
  KC_VERTEX* vs;
#endif
  PrefMask pref_pres;

#if defined(KDICT) || defined(KCOUNTER)
  KC_UC<T> uc;
#else
  KC_UC uc;
#endif

  uint16_t vs_size;
public:
  KC_VERTEX() : vs(NULL), vs_size(0) {
  }

  ~KC_VERTEX() {
    clear();
  }

  template<class Archive>
  void save(Archive& ar, const unsigned int version) const {
    ar & vs_size;
    ar & pref_pres;
    ar & uc;

    CDEPTH -= 1;
    for(size_t i = 0; i < vs_size; i++) {
      ar & vs[i];
    }
    CDEPTH += 1;
  }

 
  template<class Archive>
  void load(Archive& ar, const unsigned int version) {
    ar & vs_size;
    ar & pref_pres;
    ar & uc;

#if defined(KDICT) || defined(KCOUNTER)
    vs = new KC_VERTEX<T>[vs_size];
#else
    vs = new KC_VERTEX[vs_size];
#endif

    CDEPTH -= 1;
    for(size_t i = 0; i < vs_size; i++) {
      ar & vs[i];
    }
    CDEPTH += 1;
  }

  KC_VERTEX& operator=(KC_VERTEX&& o) {
    uc = std::move(o.uc);
    
    vs = o.vs;
    o.vs = NULL;
    
    std::swap(o.vs_size, vs_size);
    std::swap(o.pref_pres, pref_pres);
    return *this;
  }

  void clear() {
    //std::cout << "vertex clear" << std::endl;
    //delete uc;
    pref_pres = PrefMask();
    uc.clear();

    if(vs != NULL) {
      for(uint16_t i = 0; i < vs_size; i++) {
	vs[i].clear();
      }
      delete[] vs;
      vs = NULL;
      vs_size = 0;
    }
  }

#if defined(KDICT) || defined(KCOUNTER)
  KC_VERTEX<T>* get_vs() { return vs; }
#else
  KC_VERTEX* get_vs() { return vs; }
#endif

#if defined(KDICT) || defined(KCOUNTER)
  KC_UC<T>* get_uc() { return &uc; }
#else
  KC_UC* get_uc() { return &uc; }
#endif

  PrefMask get_pref_pres() { return pref_pres; }
  void set_pref_pres(const PrefMask& p) { pref_pres = p; }
  uint16_t get_vs_size() { return vs_size; }
  void set_vs_size(uint16_t vs_size) { this->vs_size = vs_size; }

#if defined(KDICT) || defined(KCOUNTER)
  void set_vs(KC_VERTEX<T>* vs) { this->vs = vs; }
#else
  void set_vs(KC_VERTEX* vs) { this->vs = vs; }
#endif

  void vertex_remove(uint8_t* bseq, int k) {
    uint8_t prefix = bseq[0];

    if((pref_pres >> (unsigned)prefix).test_bit(0)) {
      int vidx = PrefMask::popcount_vidx(pref_pres, prefix);
      vs[vidx].vertex_remove(&bseq[1], k - 4);
    }

    std::pair< bool, int > sres = uc.uc_find(k, bseq);
    int uc_idx = sres.second;
    if(sres.first) {
      uc.uc_remove(calc_bk(k), uc_idx);
      return;
    }

#if defined(PYTHON)
    throw pybind11::key_error( "Key not found!" );
#endif
  }

  uint64_t get_vertex_size()
  {
    uint64_t c = uc.get_size();
    for(int i = 0; i < vs_size; i++) {
      c += vs[i].get_vertex_size();
    }

    return c;
  }

  bool vertex_contains(uint8_t* bseq, int k)
  {
    uint8_t prefix = bseq[0];
    if((pref_pres >> (unsigned)prefix).test_bit(0)) {
      int vidx = PrefMask::popcount_vidx(pref_pres, prefix);
      return vs[vidx].vertex_contains(&bseq[1], k - 4);
    }

    std::pair<bool, int> sres = uc.uc_find(k, bseq);
    if( sres.first ) {
      return true;
    }

    return false;
  }

#if defined(KDICT) || defined(KCOUNTER)
#if defined(KDICT)
  T& vertex_get(uint8_t* bseq, int k) {
#elif defined(KCOUNTER)
  T vertex_get(uint8_t* bseq, int k) {
#endif
    uint8_t prefix = bseq[0];

    if((pref_pres >> (unsigned)prefix).test_bit(0)) {
      int vidx = PrefMask::popcount_vidx(pref_pres, prefix);
      return vs[vidx].vertex_get(&bseq[1], k - 4);
    }

    std::pair< bool, int> sres = uc.uc_find(k, bseq);
    int uc_idx = sres.second;
    if(sres.first) {
      return uc.get_obj(uc_idx);
    }

    // TODO: fix for kcounter, catch exception

#if defined(PYTHON)
#if defined(KCOUNTER)
    //int default_count = 0;
    //return default_count;
    //this->vertex_insert(bseq[0], k, default_count);
    //return this->vertex_get(bseq, k);
    return 0;
#elif defined(KDICT)
   
    throw pybind11::key_error( "Key not in dictionary!" );
#endif
#endif
  }
#endif

#if defined(KDICT) || defined(KCOUNTER)
  void vertex_insert(uint8_t* bseq, int k, T obj, std::function<T(T&, T&)>& merge_func)
#elif defined(KSET)
  void vertex_insert(uint8_t* bseq, int k)
#endif
  {
    uint8_t prefix = bseq[ 0 ];
    if((pref_pres >> (unsigned)prefix).test_bit(0)) {
      int vidx = PrefMask::popcount_vidx(pref_pres, prefix);
#if defined(KDICT) || defined(KCOUNTER)
      vs[vidx].vertex_insert(&bseq[1], k - 4, obj, merge_func);
#elif defined(KSET)
      vs[vidx].vertex_insert(&bseq[1], k - 4);
#endif
      return;
    }

    // NOTE: if the key already exisrts, we need to update it somehow
    std::pair< bool, int > sres = uc.uc_find(k, bseq);
    int uc_idx = sres.second;
    if( sres.first ) {
      // replace object here
#if defined(KDICT) || defined(KCOUNTER)
      if(merge_func != NULL) {
	//std::cout << "merging objects!" << std::endl;
	T merged_obj = merge_func(uc.get_obj(uc_idx), obj);
	uc.set_obj(uc_idx, merged_obj);
      } else {
	uc.set_obj(uc_idx, obj);
      }
#endif
      return;
    }

#if defined(KDICT) || defined(KCOUNTER)
    uc.uc_insert(bseq, k, uc_idx, obj);
#elif defined(KSET)
    uc.uc_insert(bseq, k, uc_idx);
#endif

    if(uc.get_size() == CAPACITY) {
#if defined(KDICT) || defined(KCOUNTER)
      burst_uc(k, merge_func);
#else
      burst_uc(k);
#endif
    }
  }

  void realloc_vertex_array(uint16_t insert_at = 0)
  {
#if defined(KDICT) || defined(KCOUNTER)
    KC_VERTEX<T>* vs_temp = new KC_VERTEX<T>[vs_size + 1];
#else
    KC_VERTEX* vs_temp = new KC_VERTEX[vs_size + 1];
#endif
    uint16_t insert_idx = 0, orig_idx = 0;
    for(; orig_idx < vs_size; insert_idx++, orig_idx++) {
      if(orig_idx == insert_at) {
	//std::cout << "skipping: " << orig_idx << std::endl;
	insert_idx++;
      }
      
      vs_temp[insert_idx] = std::move(vs[orig_idx]);
    }

    if(vs != NULL) {
      delete[] vs;
    }

    vs = vs_temp;
    vs_size++;
  }

#if defined(KDICT) || defined(KCOUNTER)
  void burst_uc(int k, std::function<T(T&, T&)>& merge_func)
#else
  void burst_uc(int k)
#endif
  {
    //std::cout << "burst" << std::endl;
    int suffix_size = calc_bk( k );

    uint8_t* suffixes = uc.get_suffixes();
#if defined(KDICT) || defined(KCOUNTER)
    std::vector<T> objs = uc.get_objs();
#endif
    int idx;
    for(size_t i = 0; i < uc.get_size(); i++) {
      idx = i * suffix_size;

      uint8_t* bseq = &suffixes[ idx ];
      uint8_t prefix = bseq[ 0 ];
      uint8_t* suffix = &bseq[ 1 ];
      uint8_t bits_to_shift = (unsigned) prefix;
      uint16_t vidx = PrefMask::popcount_vidx(pref_pres, prefix);

      // check if there is already a vertex that represents this prefix
      if(!(pref_pres >> (unsigned)bits_to_shift).test_bit(0)) {
	realloc_vertex_array(vidx);
	pref_pres.set_bit(bits_to_shift);
      } else {
	//std::cout << "\t" << vidx << "\tusing existing vertex" << std::endl;
      }

#if defined(KDICT) || defined(KCOUNTER)
      vs[vidx].vertex_insert(suffix, k - 4, objs[ i ], merge_func);
#elif defined(KSET)
      vs[vidx].vertex_insert(suffix, k - 4);
#endif
    }

    uc.clear();
  }
};

namespace kc_io {

#if defined(KDICT) || defined(KCOUNTER)
template <class T>
inline BinaryOutArchive& operator&(BinaryOutArchive& ar, const KC_UC<T>& uc) {
  uc.save(ar, 0);
  return ar;
}

template <class T>
inline BinaryInArchive& operator&(BinaryInArchive& ar, KC_UC<T>& uc) {
  uc.load(ar, 0);
  return ar;
}

template <class T>
inline BinaryOutArchive& operator&(BinaryOutArchive& ar, const KC_VERTEX<T>& vertex) {
  vertex.save(ar, 0);
  return ar;
}

template <class T>
inline BinaryInArchive& operator&(BinaryInArchive& ar, KC_VERTEX<T>& vertex) {
  vertex.load(ar, 0);
  return ar;
}
#else
inline BinaryOutArchive& operator&(BinaryOutArchive& ar, const KC_UC& uc) {
  uc.save(ar, 0);
  return ar;
}

inline BinaryInArchive& operator&(BinaryInArchive& ar, KC_UC& uc) {
  uc.load(ar, 0);
  return ar;
}

inline BinaryOutArchive& operator&(BinaryOutArchive& ar, const KC_VERTEX& vertex) {
  vertex.save(ar, 0);
  return ar;
}

inline BinaryInArchive& operator&(BinaryInArchive& ar, KC_VERTEX& vertex) {
  vertex.load(ar, 0);
  return ar;
}
#endif

}  // namespace kc_io

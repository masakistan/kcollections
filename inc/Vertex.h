#pragma once

#if defined(PYTHON)
#include <pybind11/pybind11.h>
#endif

#include <functional>
#include "UContainer.h"
#include "globals.h"
//#include <jemalloc/jemalloc.h>
#include "uint256_t.h"

#if defined(PYTHON)
namespace py = pybind11;
#endif

#if defined(KDICT) || defined(KCOUNTER)
template <class T>
#endif
class Vertex {
private:
#if defined(KDICT) || defined(KCOUNTER)
  //Vertex<T>** vs;
  std::vector<Vertex<T>> vs;
#else
  //Vertex** vs;
  std::vector<Vertex> vs;
#endif
  uint256_t pref_pres;
#if defined(KDICT) || defined(KCOUNTER)
  UC<T> uc;
#else
  UC uc;
#endif
  //uint16_t vs.size();
  int k;
  
public:
  Vertex(int k) : uc(k), k(k) {
    pref_pres = 0;
  }

  Vertex(const Vertex& in_v) :
    k(in_v.k),
    uc(in_v.uc),
    vs(in_v.vs)
  {
  }

  ~Vertex() {
    pref_pres = 0;
    vs.clear();
  }

  int calc_vidx(uint256_t vertices, uint8_t bts) {
    vertices <<= (256 - (unsigned) bts);
    uint64_t* t = (uint64_t*) &vertices;
    int vidx = __builtin_popcountll(t[0]);
    vidx += __builtin_popcountll(t[1]);
    vidx += __builtin_popcountll(t[2]);
    vidx += __builtin_popcountll(t[3]);
    return vidx;
  }

#if defined(KDICT) || defined(KCOUNTER)
  std::vector<Vertex<T>>& get_vs() { return vs; }
#else
  std::vector<Vertex>& get_vs() { return vs; }
#endif

#if defined(KDICT) || defined(KCOUNTER)
  UC<T>* get_uc() { return &uc; }
#else
  UC* get_uc() { return &uc; }
#endif

  uint256_t get_pref_pres() { return pref_pres; }
  void set_pref_pres(uint256_t pref_pres) { this->pref_pres = pref_pres; }
  size_t get_vs_size() { return vs.size(); }

  void vertex_remove(uint8_t* bseq, int k) {
    uint8_t prefix = bseq[0];

    if((pref_pres >> (unsigned) prefix) & 0x1) {
      int vidx = calc_vidx(pref_pres, prefix);
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
    for(size_t i = 0; i < vs.size(); i++) {
      c += vs[i].get_vertex_size();
    }

    return c;
  }

  bool vertex_contains(uint8_t* bseq, int k)
  {
    uint8_t prefix = bseq[0];
    if((pref_pres >> (unsigned) prefix) & 0x1) {
      int vidx = calc_vidx(pref_pres, prefix);
      return vs[vidx].vertex_contains(&bseq[1], k - 4);
    }

    std::pair<bool, int> sres = uc.uc_find(k, bseq);
    if( sres.first ) {
      return true;
    }

    return false;
  }

#if defined(KDICT) || defined(KCOUNTER)
  T& vertex_get(uint8_t* bseq, int k) {
    uint8_t prefix = bseq[0];

    if((pref_pres >> (unsigned) prefix) & 0x1) {
      int vidx = calc_vidx(pref_pres, prefix);
      return vs[vidx].vertex_get(&bseq[1], k - 4);
    }

    std::pair< bool, int> sres = uc.uc_find(k, bseq);
    int uc_idx = sres.second;
    if(sres.first) {
      return uc.get_obj(uc_idx);
    }

    // TODO: fix for kcounter, catch exception
#if defined(PYTHON)
    throw pybind11::key_error( "Key not in dictionary!" );
#endif
  }
#endif

#if defined(KDICT) || defined(KCOUNTER)
  void vertex_insert(uint8_t* bseq, int k, T& obj, std::function<T(T&, T&)>& merge_func)
#elif KSET
  void vertex_insert(uint8_t* bseq, int k)
#endif
  {
    uint8_t prefix = bseq[ 0 ];
    if((pref_pres >> (unsigned) prefix) & 0x1) {
      int vidx = calc_vidx(pref_pres, prefix);
#if defined(KDICT) || defined(KCOUNTER)
      vs[vidx].vertex_insert(&bseq[1], k - 4, obj, merge_func);
#elif KSET
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
	//std::cout << "merging\t" << std::endl;
	//std::cout << uc.get_obj(uc_idx) << std::endl;
	//std::cout << obj << std::endl;
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
#elif KSET
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

#if defined(KDICT) || defined(KCOUNTER)
  void burst_uc(int k, std::function<T(T&, T&)>& merge_func)
#else
  void burst_uc(int k)
#endif
  {
    //std::cout << "bursting" << std::endl;
    int suffix_size = calc_bk( k );

    uint8_t* suffixes = uc.get_suffixes();
#if defined(KDICT) || defined(KCOUNTER)
    std::vector<T>& objs = uc.get_objs();
#endif
    int idx;
    for(size_t i = 0; i < uc.get_size(); i++) {
      idx = i * suffix_size;

      uint8_t* bseq = &suffixes[ idx ];
      uint8_t prefix = bseq[ 0 ];
      uint8_t* suffix = &bseq[ 1 ];
      uint8_t bits_to_shift = (unsigned) prefix;
      int vidx = calc_vidx(pref_pres, prefix);

      // check if there is already a vertex that represents this prefix
      if(!((pref_pres >> (unsigned) bits_to_shift) & 0x1)) {
	pref_pres |= ((uint256_t) 0x1 << (unsigned) bits_to_shift);
	vs.reserve(vs.size() + 1);
	vs.insert(vs.begin() + vidx, Vertex(k - 4));
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

#pragma once

#include <pybind11/pybind11.h>
#include "UContainer.h"
#include "globals.h"
#include <jemalloc/jemalloc.h>
#include "uint256_t.h"

namespace py = pybind11;

template <class T>
class Vertex {
private:
  Vertex<T>** vs;
  uint256_t pref_pres;
  UC<T>* uc;
  uint16_t vs_size;
public:
  Vertex() : vs(NULL), vs_size(0) {
    pref_pres = 0;
    uc = new UC<T>();
  }

  ~Vertex() {
    delete uc;
    pref_pres = 0;

    if(vs != NULL) {
      for(uint16_t i = 0; i < vs_size; i++) {
	delete vs[i];
      }
      free(vs);
      vs_size = 0;
    }
  }

  Vertex<T>** get_vs() { return vs; }
  UC<T>* get_uc() { return uc; }
  uint256_t get_pref_pres() { return pref_pres; }
  void set_pref_pres(uint256_t pref_pres) { this->pref_pres = pref_pres; }
  uint16_t get_vs_size() { return vs_size; }
  void set_vs_size(uint16_t vs_size) { this->vs_size = vs_size; }
  void set_vs(Vertex<T>** vs) { this->vs = vs; }
  
  void vertex_remove(uint8_t* bseq, int k) {
    uint8_t prefix = bseq[0];
    
    if((pref_pres >> (unsigned) prefix) & 0x1) {
      int vidx = calc_vidx(pref_pres, prefix);
      vs[vidx]->vertex_remove(&bseq[1], k - 4);
    }

    std::pair< bool, int > sres = uc->uc_find(k, bseq);
    int uc_idx = sres.second;
    if(sres.first) {
      uc->uc_remove(calc_bk(k), uc_idx);
      return;
    }

    throw pybind11::key_error( "Key not found!" );
  }
  
  uint64_t get_vertex_size()
  {
    uint64_t c = uc->get_size();
    for(int i = 0; i < vs_size; i++) {
      c += vs[i]->get_vertex_size();
    }
    
    return c;
  }
  
  bool vertex_contains(uint8_t* bseq, int k)
  {
    uint8_t prefix = bseq[0];
    if((pref_pres >> (unsigned) prefix) & 0x1) {
      int vidx = calc_vidx(pref_pres, prefix);
      return vs[vidx]->vertex_contains(&bseq[1], k - 4);
    }
    
    std::pair<bool, int> sres = uc->uc_find(k, bseq);
    if( sres.first ) {
      return true;
    }
    
    return false;
  }

#if defined(KDICT) || defined(KCOUNTER)
  int vertex_get(uint8_t* bseq, int k) {
    uint8_t prefix = bseq[0];
    
    if((pref_pres >> (unsigned) prefix) & 0x1) {
      int vidx = calc_vidx(pref_pres, prefix);
      return vs[vidx]->vertex_get(&bseq[1], k - 4);
    }
    
    std::pair< bool, int> sres = uc->uc_find(k, bseq);
    int uc_idx = sres.second;
    if(sres.first) {
      return uc->get_obj(uc_idx);
    }
#if KCOUNTER
    // add a vertex to be 0 if key is not found
    //int default_count = 0;
    //vertex_insert(bseq, k, default_count);
    return 0;
#endif

    throw pybind11::key_error( "Key not in dictionary!" );
  }
#endif

#if defined(KDICT) || defined(KCOUNTER)
  void vertex_insert(uint8_t* bseq, int k, int obj, std::function<int(int, int)>& merge_func)
#elif KSET
  void vertex_insert(uint8_t* bseq, int k)
#endif
  {
    uint8_t prefix = bseq[ 0 ];
    if((pref_pres >> (unsigned) prefix) & 0x1) {
      int vidx = calc_vidx(pref_pres, prefix);
#if defined(KDICT) || defined(KCOUNTER)
      vs[vidx]->vertex_insert(&bseq[1], k - 4, obj, merge_func);
#elif KSET
      vs[vidx]->vertex_insert(&bseq[1], k - 4);
#endif
      return;
    }
    
    // NOTE: if the key already exisrts, we need to update it somehow
    std::pair< bool, int > sres = uc->uc_find(k, bseq);
    int uc_idx = sres.second;
    if( sres.first ) {
      // replace object here
#if defined(KDICT) || defined(KCOUNTER)
      if(merge_func != NULL) {
	int merged_obj = merge_func(uc->get_obj(uc_idx), obj);
	uc->set_obj(uc_idx, merged_obj);
      } else {
	uc->set_obj(uc_idx, obj);
      }
#endif
      return;
    }
    
#if defined(KDICT) || defined(KCOUNTER)
    uc->uc_insert(bseq, k, uc_idx, obj);
#elif KSET
    uc->uc_insert(bseq, k, uc_idx);
#endif
    
    if(uc->get_size() == CAPACITY) {
#if defined(KDICT) || defined(KCOUNTER)
      burst_uc(k, merge_func);
#else
      burst_uc(k);
#endif
    }
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
  void burst_uc(int k, std::function<int(int, int)>& merge_func)
#else
  void burst_uc(int k)
#endif
  {
    int suffix_size = calc_bk( k );
    
    uint8_t* suffixes = uc->get_suffixes();
#if defined(KDICT) || defined(KCOUNTER)
    std::vector<int> objs = uc->get_objs();
#endif
    int idx;
    for(size_t i = 0; i < uc->get_size(); i++) {
      idx = i * suffix_size;
      
      uint8_t* bseq = &suffixes[ idx ];
      uint8_t prefix = bseq[ 0 ];
      uint8_t* suffix = &bseq[ 1 ];
      uint8_t bits_to_shift = (unsigned) prefix;
      int vidx = calc_vidx(pref_pres, prefix);
      
      // check if there is already a vertex that represents this prefix
      if(!((pref_pres >> (unsigned) bits_to_shift) & 0x1)) {
	
	if(vs == NULL) {
	  vs = (Vertex<T>**) calloc(1, sizeof(Vertex<T>*));
	} else {
	  vs = (Vertex<T>**) realloc(vs, (vs_size + 1) * sizeof(Vertex<T>*));
	}
	// move any previous vertices if necessary
	if(vidx < vs_size) {
	  //std::cout << "moving " << vidx << " to " << vidx + 1 << " (size of array " << vs_size + 1 << ") moving " << (vs_size - vidx) << " total items" << std::endl;
	  std::memmove(&vs[vidx + 1],
		       &vs[vidx],
		       (vs_size - vidx) * sizeof(Vertex<T>*)
		       );
	}
	
	pref_pres |= ((uint256_t) 0x1 << (unsigned) bits_to_shift);
	// insert a vertex at vidx
	vs[vidx] = new Vertex<T>();
	
	// increment size
	vs_size++;
      } else {
	//std::cout << "\t" << vidx << "\tusing existing vertex" << std::endl;
      }
      
#if defined(KDICT) || defined(KCOUNTER)
      vs[vidx]->vertex_insert(suffix, k - 4, objs[ i ], merge_func);
#elif KSET
      vs[vidx]->vertex_insert(suffix, k - 4);
#endif
    }
    
    delete uc;
    uc = new UC<T>();
  }
};


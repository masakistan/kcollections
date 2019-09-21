#pragma once

#include <utility>
#include <stdio.h>
#include <iostream>
#include <vector>

#include "globals.h"
#include "helper.h"
//#include <pybind11/pybind11.h>
//#include <jemalloc/jemalloc.h>

//namespace py = pybind11;

#if defined(KDICT) || defined(KCOUNTER)
template <class T>
#endif

class UC {
private:
#if defined(KDICT) || defined(KCOUNTER)
  std::vector<T> objs;
#else
  size_t size;
#endif
  uint8_t* suffixes;

public:
  UC() : suffixes(NULL) {
#if defined(KSET)
    size = 0;
#endif
  }

  ~UC() {
    if(suffixes != NULL) {
      free(suffixes);
#if defined(KDICT) || defined(KCOUNTER)
      objs.clear();
#else
      size = 0;
#endif
    }
  }

  /*
  void print(int k) {
    int len = calc_bk(k);
    int idx;
    for(size_t i = 0; i < get_size(); i++) {
      idx = i * len;
      char* dseq = deserialize_kmer(k, len, suffixes[idx]);
      std::cout << "kmer: " << dseq << std::endl;
      free(dseq);
    }
  }
  */

#if defined(KDICT) || defined(KCOUNTER)
void uc_insert(uint8_t* bseq, int k, int idx, T obj)
#else
void uc_insert(uint8_t* bseq, int k, int idx)
#endif
{
  int len = calc_bk(k);
  if(suffixes == NULL) {
    suffixes = (uint8_t*) calloc(len, sizeof(uint8_t));
  }
  else {
    suffixes = (uint8_t*) realloc(
				  suffixes,
				  len * (get_size() + 1) * sizeof(uint8_t)
				  );
  }

  if(get_size() < CAPACITY) {
    int bytes_to_move = (get_size() - idx) * len;
    int suffix_idx = idx * len;
    if(bytes_to_move > 0) {
      std::memmove(
		   &suffixes[suffix_idx + len],
		   &suffixes[suffix_idx],
		   bytes_to_move
		   );
    }

    std::memcpy(&suffixes[suffix_idx], bseq, len);

#if defined(KDICT) || defined(KCOUNTER)
    objs.insert(objs.begin() + idx, obj);
#else
    size++;
#endif
  } else {
    std::cout << "this is a mistake!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  }
}

  std::pair< bool, int > uc_find(int k, uint8_t* bseq) {
    if(suffixes == NULL) {
      return std::make_pair(false, get_size());
    }

    return binary_search(suffixes, get_size(), calc_bk(k), bseq);
  }

  void uc_remove(int bk, int idx) {
    int suffix_idx = idx * bk;
    int bytes_to_move = (get_size() - (idx + 1)) * bk;
    std::memmove(
		 &suffixes[ suffix_idx ],
		 &suffixes[ suffix_idx + bk ],
		 bytes_to_move
		 );
#if defined(KDICT) || defined(KCOUNTER)
    objs.erase(objs.begin() + idx);
#else
    size--;
#endif
}

  size_t get_size() {
#if defined(KDICT) || defined(KCOUNTER)
    return objs.size();
#else
    return size;
#endif
  }


  uint8_t* get_suffixes() {
    return suffixes;
  }

  uint8_t* get_suffix(int k, int idx) {
    int bk = calc_bk( k );
    int suffix_idx = bk * idx;
    return &suffixes[suffix_idx];
  }

#if defined(KDICT) || defined(KCOUNTER)
  T get_obj(int obj_idx) {
    return objs[obj_idx];
  }

  std::vector<T>& get_objs() {
    return objs;
  }

  void set_obj(int obj_idx, T obj) {
    objs[obj_idx] = obj;
  }
#endif
};

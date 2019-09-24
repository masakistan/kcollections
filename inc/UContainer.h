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
  int k;

public:
  UC(int k) : suffixes(NULL), k(k) {
#if defined(KSET)
    size = 0;
#endif
  }

  UC(const UC& in_uc) {
#if defined(KDICT) || defined(KCOUNTER)
    objs = std::vector<T>(in_uc.objs);
    //objs = std::vector<int>(in_uc.objs.size(), 1); //in_uc.objs;
#else
    size = in_uc.size;
#endif
    
    k = in_uc.k;
    if(get_size() > 0) {
      int len = calc_bk(k);
      //std::cout << len << "\t" << get_size() << std::endl;
      suffixes = (uint8_t*) malloc(
				   sizeof(uint8_t) * len * get_size()
				   );
      memcpy(
	     suffixes,
	     in_uc.suffixes,
	     sizeof(uint8_t) * len * get_size()
	     );
    } else {
      //std::cout << "setting suffixes to null" << std::endl;
      suffixes = NULL;
    }
  }

  ~UC() {
    clear();
  }

  void clear() {
    if(get_size() > 0) {
      free(suffixes);
      suffixes = NULL;
#if defined(KDICT) || defined(KCOUNTER)
      //objs = std::vector<T>();
      objs.clear();
      objs.shrink_to_fit();
#else
      size = 0;
#endif
    }
  }

#if defined(KDICT) || defined(KCOUNTER)
void uc_insert(uint8_t* bseq, int k, int idx, T& obj)
#else
void uc_insert(uint8_t* bseq, int k, int idx)
#endif
{
  assert(this->k == k);
  int len = calc_bk(k);
  if(get_size() == 0) {
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
    //objs = std::vector<int>(objs.size(), 1);
    objs.reserve(get_size() + 1);


    //std::cout << objs.capacity() << std::endl;
    //std::cout << "inserting at " << idx << "\t" << obj << "\t" << CAPACITY<< std::endl;
    //assert(idx < objs.capacity());
    objs.insert(objs.begin() + idx, obj);
    //objs.push_back(T(obj));
#else
    size++;
#endif
  } else {
    std::cout << "this is a mistake!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  }
}

  std::pair< bool, int > uc_find(int k, uint8_t* bseq) {
    if(get_size() == 0) {
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
    objs.shrink_to_fit();
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
  T& get_obj(int obj_idx) {
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

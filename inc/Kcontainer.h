#pragma once

#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <pthread.h>
#include <semaphore.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <vector>
#include <stdlib.h>
#include <string>
#include <stdexcept>
#include "Vertex.h"
#include "helper.h"
#include <jemalloc/jemalloc.h>
#include <math.h>

namespace py = pybind11;


typedef struct {
  int k;
  Vertex<int>* v;
} Kcontainer;

/*
class Kcontainer {
public:
  Kcontainer(const int k);
  ~Kcontainer();
  int k;
  Vertex v;
};
*/

inline void init_kcontainer(Kcontainer* kd, int k)
{
  kd->k = k;
  kd->v = new Vertex<int>();
}

inline Kcontainer* create_kcontainer(int k)
{
  srand(1337);
  Kcontainer* kd = ( Kcontainer* ) malloc( sizeof( Kcontainer ) );
  init_kcontainer(kd, k);
  return kd;
}

inline void free_kcontainer( Kcontainer* kd )
{
  delete kd->v;
  free(kd);
}

inline std::string kcontainer_get_child_suffix(Vertex<int>* v, int idx) {
  uint256_t verts = v->get_pref_pres();
  uint8_t j = 0, i = 0;
  while(true) {
    if(verts & 0x1) {
      j++;
    }

    if(j > idx)
      break;

    if((unsigned) i == 255) {
      break;
    }

    verts >>= 1;
    i++;
  }
  char* kmer = deserialize_kmer(4, &i);
  std::string skmer(kmer);
  free(kmer);
  return skmer;
}

inline bool kcontainer_contains( Kcontainer* kd, const char* kmer )
{
  uint8_t* bseq = ( uint8_t* ) calloc( kd->k, sizeof( uint8_t ) );
  serialize_kmer( kmer, kd->k, bseq );
  bool res = kd->v->vertex_contains(bseq, kd->k);
  free( bseq );
  return res;
  return false;
}

#if defined(KDICT) || defined(KCOUNTER)
inline void kcontainer_add(Kcontainer* kd, const char* kmer, int obj, std::function<int(int, int)>& merge_func)
#elif KSET
  inline void kcontainer_add(Kcontainer* kd, const char* kmer)
#endif
{
  uint8_t* bseq = ( uint8_t* ) calloc( kd->k, sizeof( uint8_t ) );
  serialize_kmer( kmer, kd->k, bseq );
  //std::cout << deserialize_kmer(kd->k, bseq) << std::endl;
#if defined(KDICT) || defined(KCOUNTER)
  kd->v->vertex_insert(bseq, kd->k, obj, merge_func);
#elif KSET
  kd->v->vertex_insert(bseq, kd->k);
#endif
  free(bseq);
}

#if defined KDICT || defined KCOUNTER
inline int kcontainer_get( Kcontainer* kd, char* kmer )
{
  uint8_t* bseq = ( uint8_t* ) calloc( kd->k, sizeof( uint8_t ) );
  serialize_kmer( kmer, kd->k, bseq );
#if defined(KDICT) || defined(KCOUNTER)
  int res = kd->v->vertex_get(bseq, kd->k);
#endif
  free( bseq );
  return res;
}
#endif

inline uint64_t kcontainer_size( Kcontainer* kd )
{
  return kd->v->get_vertex_size();
}

inline void kcontainer_remove( Kcontainer* kd, const char* kmer )
{
  uint8_t* bseq = ( uint8_t* ) calloc( kd->k, sizeof( uint8_t ) );
  serialize_kmer( kmer, kd->k, bseq );
  kd->v->vertex_remove(bseq, kd->k);
  free( bseq );
}

#if defined(KSET)
void parallel_kcontainer_add_init(Kcontainer* kd, int threads);
#elif defined(KDICT) || defined(KCOUNTER)
void parallel_kcontainer_add_init(Kcontainer* kd, int threads, const std::function<int(int, int)> &f );
#endif

#if defined(KSET) || defined(KCOUNTER)
void parallel_kcontainer_add(Kcontainer* kd, const char* kmer);
#elif defined(KDICT)
void parallel_kcontainer_add(Kcontainer* kd, const char* kmer, int value);
#endif

void* parallel_kcontainer_add_consumer(void* bin_ptr);
void parallel_kcontainer_add_join(Kcontainer* kc);
#if defined(KSET) || defined(KCOUNTER)
void parallel_kcontainer_add_seq(Kcontainer* kd, const char* seq, uint32_t length);
#elif defined(KDICT)
void parallel_kcontainer_add_seq(Kcontainer* kd, const char* seq, uint32_t length, py::iterable values);
#endif

#if defined(KSET)
void parallel_kcontainer_add_bseq(uint8_t* bseq);
#elif defined(KDICT) || defined(KCOUNTER)
void parallel_kcontainer_add_bseq(uint8_t* bseq, int value);
#endif

#if defined(KDICT)
inline void kcontainer_add_seq(Kcontainer* kd, const char* seq, uint32_t length, py::iterable values, std::function<int(int, int)> &f)
#elif defined(KCOUNTER)
inline void kcontainer_add_seq(Kcontainer* kd, const char* seq, uint32_t length, std::function<int(int, int)> &f)
#else
inline void kcontainer_add_seq(Kcontainer* kd, const char* seq, uint32_t length)
#endif
{
    int size64 = kd->k / 32;
    if(kd->k % 32 > 0) {
      size64++;
    }

    uint64_t* bseq64 = (uint64_t*) calloc(size64, sizeof(uint64_t));
    uint8_t* bseq8 = (uint8_t*) bseq64;

    uint bk = calc_bk(kd->k);
    uint8_t last_index = (kd->k - 1) % 4;

    // serialize the first kmer
    serialize_kmer(seq, kd->k, bseq8);

#if KSET
    kd->v->vertex_insert(bseq8, kd->k);
#elif defined(KCOUNTER)
    kd->v->vertex_insert(bseq8, kd->k, 1, f);
#elif defined(KDICT)
    auto iter = py::iter(values);
    kd->v->vertex_insert(bseq8, kd->k, (*iter).cast<int>(), f);
#endif

    //std::cout << strlen(seq) << std::endl;
    for(uint32_t j = kd->k; j < length; j++) {
      //std::cout << j << "\t" << seq[j] << std::endl;
      // shift all the bits over
      //bseq8[0] <<= 2;
      bseq64[0] >>= 2;
      //std::cout << "shifting\t" << deserialize_kmer(kd->k, bseq8) << std::endl;
      //for(int i = 1; i < bk; i++) {
      for(int i = 1; i < size64; i++) {
	bseq64[i - 1] |= (bseq64[i] << 62);
	bseq64[i] >>= 2;
	//std::cout << "shifting\t" << deserialize_kmer(kd->k, bseq8) << std::endl;
      }
      
      serialize_position(j, bk - 1, last_index, bseq8, seq);
#if KSET
      kd->v->vertex_insert(bseq8, kd->k);
#elif defined(KCOUNTER)
      kd->v->vertex_insert(bseq8, kd->k, 1, f);
#elif defined(KDICT)
      std::advance(iter, 1);
      kd->v->vertex_insert(bseq8, kd->k, (*iter).cast<int>(), f);
#endif
    }

    free(bseq64);

  }

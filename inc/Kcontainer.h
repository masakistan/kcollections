#pragma once

#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
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
    Vertex v;
} Kcontainer;



inline void init_kcontainer(Kcontainer* kd, int k)
{
    kd->k = k;
    init_vertex(&(kd->v));
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
    free_vertex(&(kd->v));
    free(kd);
}

inline bool kcontainer_contains( Kcontainer* kd, char* kmer )
{
    uint8_t* bseq = ( uint8_t* ) calloc( kd->k, sizeof( uint8_t ) );
    serialize_kmer( kmer, kd->k, bseq );
    bool res = vertex_contains( &( kd->v ), bseq, kd->k, 0 );
    free( bseq );
    return res;
    return false;
}

#if KDICT
inline void kcontainer_add( Kcontainer* kd, const char* kmer, py::handle* obj )
#elif KSET
inline void kcontainer_add( Kcontainer* kd, const char* kmer )
#elif KCOUNTER
inline void kcontainer_add( Kcontainer* kd, const char* kmer, int count )
#endif
{
    uint8_t* bseq = ( uint8_t* ) calloc( kd->k, sizeof( uint8_t ) );
    serialize_kmer( kmer, kd->k, bseq );
#if KDICT
    vertex_insert( &( kd->v ), bseq, kd->k, 0, obj );
#elif KSET
    vertex_insert( &( kd->v ), bseq, kd->k, 0 );
#elif KCOUNTER
    vertex_insert( &( kd->v ), bseq, kd->k, 0, count );
#endif
    free( bseq );
}

#if defined KDICT || defined KCOUNTER
#if KDICT
inline py::handle* kcontainer_get( Kcontainer* kd, char* kmer )
#elif KCOUNTER
inline int kcontainer_get( Kcontainer* kd, char* kmer )
#endif
{
    uint8_t* bseq = ( uint8_t* ) calloc( kd->k, sizeof( uint8_t ) );
    serialize_kmer( kmer, kd->k, bseq );
#if KDICT
    py::handle* res = vertex_get( &kd->v, bseq, kd->k, 0 );
#elif KCOUNTER
    int res = vertex_get_counter( &kd->v, bseq, kd->k, 0 );
#endif
    //std::cout << "kcontainer get start" << std::endl;
    //py::print( py::str( *res ) );
    //std::cout << "kcontainer get end" << std::endl;
    free( bseq );
    return res;
}
#endif

inline uint64_t kcontainer_size( Kcontainer* kd )
{
    return vertex_size( &kd->v );
}

inline void kcontainer_remove( Kcontainer* kd, char* kmer )
{
    uint8_t* bseq = ( uint8_t* ) calloc( kd->k, sizeof( uint8_t ) );
    serialize_kmer( kmer, kd->k, bseq );
    vertex_remove( &kd->v, bseq, kd->k, 0 );
    free( bseq );
}

#if defined KSET || defined KCOUNTER
void parallel_kcontainer_add_init(Kcontainer* kd, int threads);
void parallel_kcontainer_add(Kcontainer* kd, const char* kmer);
void* parallel_kcontainer_add_consumer(void* bin_ptr);
void parallel_kcontainer_add_join(Kcontainer* kc);
void parallel_kcontainer_add_seq(Kcontainer* kd, const char* seq, uint32_t length);
void parallel_kcontainer_add_bseq(Kcontainer* kd, uint8_t* bseq);

inline void kcontainer_add_seq(Kcontainer* kd, const char* seq, uint32_t length) {
    int size64 = kd->k / 32;
    if(kd->k % 32 > 0) {
        size64++;
    }

    uint64_t* bseq64 = (uint64_t*) calloc(size64, sizeof(uint64_t));
    uint8_t* bseq8 = (uint8_t*) bseq64;

    uint bk = calc_bk(kd->k);
    uint8_t holder;
    uint8_t last_index = (kd->k - 1) % 4;

    // serialize the first kmer
    serialize_kmer(seq, kd->k, bseq8);

#if KSET
    vertex_insert(&(kd->v), bseq8, kd->k, 0);
#elif KCOUNTER
    // get current count of bseq
    int count = vertex_get_counter(&(kd->v), bseq8, kd->k, 0);
    vertex_insert(&(kd->v), bseq8, kd->k, 0, ++count);
#endif

    for(int j = kd->k; j < length; j++) {
        //std::cout << j << "\t" << seq[j] << std::endl;
        // shift all the bits over
        bseq64[0] >>= 2;
        for(int i = 1; i < size64; i++) {
            bseq64[i - 1] |= bseq64[i] << 62;
            bseq64[i] >>= 2;
        }
      
        serialize_position(j, bk - 1, last_index, bseq8, seq);
        vertex_insert(&(kd->v), bseq8, kd->k, 0);
#elif KCOUNTER
        // get current count of bseq
        count = vertex_get_counter(&(kd->v), bseq8, kd->k, 0);
        vertex_insert(&(kd->v), bseq8, kd->k, 0, ++count);
#endif
    }

    free(bseq64);

}
#endif

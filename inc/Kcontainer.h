#pragma once

#include <utility>
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
#include <roaring/roaring.h>

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

inline char* kcontainer_get_child_suffix(Vertex* v, int idx) {
    uint256_t verts = v->pref_pres;
    uint8_t j = 0, i = 0;
    for(i = 0; i < 256; i++) {
        if(verts & 0x1) {
            //i++;
            j++;
        }
        if(j > idx)
            break;
        verts >>= 1;
    }
    return deserialize_kmer(4, 1, &i);
}

inline bool kcontainer_contains( Kcontainer* kd, const char* kmer )
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
#elif KCOLOR
inline void kcontainer_add( Kcontainer* kd, const char* kmer, uint32_t color )
#endif
{
    uint8_t* bseq = ( uint8_t* ) calloc( kd->k, sizeof( uint8_t ) );
    serialize_kmer( kmer, kd->k, bseq );
    //std::cout << deserialize_kmer(kd->k, calc_bk(kd->k), bseq) << std::endl;
#if KDICT
    vertex_insert( &( kd->v ), bseq, kd->k, 0, obj );
#elif KSET
    vertex_insert( &( kd->v ), bseq, kd->k, 0 );
#elif KCOUNTER
    vertex_insert( &( kd->v ), bseq, kd->k, 0, count );
#elif KCOLOR
    roaring_bitmap_t* r = roaring_bitmap_of(1, color);
    /*uint64_t c = roaring_bitmap_get_cardinality(r);
    uint32_t* a = (uint32_t*) malloc(c * sizeof(uint32_t));
    roaring_bitmap_to_uint32_array(r, a);
    std::cout << "insert color:\t" << color << std::endl;
    std::cout << "before insert" << std::endl;
    for(int i = 0; i < c; i++) {
        std::cout << a[i] << std::endl;
    }
    std::cout << "done" << std::endl;*/
    vertex_insert( &( kd->v ), bseq, kd->k, 0, r);
#endif
    free( bseq );
}

#if defined KDICT || defined KCOUNTER || defined KCOLOR
#if KDICT
inline py::handle* kcontainer_get( Kcontainer* kd, const char* kmer )
#elif KCOUNTER
inline int kcontainer_get( Kcontainer* kd, const char* kmer )
#elif KCOLOR
inline uint32_t* kcontainer_get(Kcontainer* kd, const char* kmer)
#endif
{
    uint8_t* bseq = ( uint8_t* ) calloc( kd->k, sizeof( uint8_t ) );
    serialize_kmer( kmer, kd->k, bseq );
#if KDICT
    py::handle* res = vertex_get( &kd->v, bseq, kd->k, 0 );
#elif KCOUNTER
    int res = vertex_get_counter( &kd->v, bseq, kd->k, 0 );
#elif KCOLOR
    uint32_t* res = vertex_get_colors(&kd->v, bseq, kd->k, 0);
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

inline void kcontainer_remove( Kcontainer* kd, const char* kmer )
{
    uint8_t* bseq = ( uint8_t* ) calloc( kd->k, sizeof( uint8_t ) );
    serialize_kmer( kmer, kd->k, bseq );
    vertex_remove( &kd->v, bseq, kd->k, 0 );
    free( bseq );
}

#if defined KSET || defined KCOUNTER || defined KCOLOR
void parallel_kcontainer_add_init(Kcontainer* kd, int threads);
void* parallel_kcontainer_add_consumer(void* bin_ptr);
void parallel_kcontainer_add_join(Kcontainer* kc);

#if defined KSET || defined KCOUNTER
void parallel_kcontainer_add(Kcontainer* kd, const char* kmer);
void parallel_kcontainer_add_seq(Kcontainer* kd, const char* seq, uint32_t length);
void parallel_kcontainer_add_bseq(Kcontainer* kd, uint8_t* bseq);
void kcontainer_add_seq(Kcontainer* kd, const char* seq, uint32_t length);
#elif defined KCOLOR
void parallel_kcontainer_add(Kcontainer* kd, const char* kmer, uint32_t color);
void parallel_kcontainer_add_seq(Kcontainer* kd, const char* seq, uint32_t length, uint32_t color);
void parallel_kcontainer_add_bseq(Kcontainer* kd, uint8_t* bseq, uint32_t color);
void kcontainer_add_seq(Kcontainer* kd, const char* seq, uint32_t length, uint32_t color);
#endif

#endif

#pragma once

#include <stdlib.h>
#include "Vertex.h"
#include "helper.h"
#include <jemalloc/jemalloc.h>

namespace py = pybind11;

typedef struct {
    int k;
    Vertex v;
} Kcontainer;

inline void init_kcontainer( Kcontainer* kd, int k )
{
    kd->k = k;
    init_vertex( &( kd->v ) );
}

inline Kcontainer* create_kcontainer( int k )
{
    Kcontainer* kd = ( Kcontainer* ) malloc( sizeof( Kcontainer ) );
    init_kcontainer( kd, k );
    return kd;
}

inline void free_kcontainer( Kcontainer* kd )
{
    free_vertex( &( kd->v ) );
    free( kd );
}

#if KDICT
inline void kcontainer_insert( Kcontainer* kd, char* kmer, py::object* obj )
#elif KSET
inline void kcontainer_insert( Kcontainer* kd, char* kmer )
#endif
{
    uint8_t* bseq = ( uint8_t* ) calloc( kd->k, sizeof( uint8_t ) );
    serialize_kmer( kmer, kd->k, bseq );
#if KDICT
    vertex_insert( &( kd->v ), bseq, kd->k, 0, obj );
#elif KSET
    vertex_insert( &( kd->v ), bseq, kd->k, 0 );
#endif
    free( bseq );
}

#if KDICT
inline py::object* kcontainer_get( Kcontainer* kd, char* kmer )
{
    uint8_t* bseq = ( uint8_t* ) calloc( kd->k, sizeof( uint8_t ) );
    serialize_kmer( kmer, kd->k, bseq );
    return vertex_get( &kd->v, bseq, kd->k, 0 );
}
#endif

inline bool kcontainer_contains( Kcontainer* kd, char* kmer )
{
    uint8_t* bseq = ( uint8_t* ) calloc( kd->k, sizeof( uint8_t ) );
    serialize_kmer( kmer, kd->k, bseq );
    bool res = vertex_contains( &( kd->v ), bseq, kd->k, 0 );
    free( bseq );
    return res;
}

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



#pragma once

#include <stdlib.h>
#include "Vertex.h"
#include "helper.h"

typedef struct {
    int k;
    Vertex v;
} Kdict;

inline void init_kdict( Kdict* kd, int k )
{
    kd->k = k;
    init_vertex( &( kd->v ) );
}

inline Kdict* create_kdict( int k )
{
    Kdict* kd = ( Kdict* ) malloc( sizeof( Kdict ) );
    init_kdict( kd, k );
    return kd;
}

inline void free_kdict( Kdict* kd )
{
    free_vertex( &( kd->v ) );
    free( kd );
}

inline void insert( Kdict* kd, char* kmer )
{
    uint8_t* bseq = ( uint8_t* ) calloc( kd->k, sizeof( uint8_t ) );
    serialize_kmer( kmer, kd->k, bseq );
    vertex_insert( &( kd->v ), bseq, kd->k, 0 );
    free( bseq );
}

inline bool contains( Kdict* kd, char* kmer )
{
    uint8_t* bseq = ( uint8_t* ) calloc( kd->k, sizeof( uint8_t ) );
    serialize_kmer( kmer, kd->k, bseq );
    bool res = vertex_contains( &( kd->v ), bseq, kd->k, 0 );
    free( bseq );
    return res;
}



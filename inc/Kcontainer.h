#pragma once

#include <stdlib.h>
#include <string>
#include <stdexcept>
#include "Vertex.h"
#include "helper.h"
#include <jemalloc/jemalloc.h>

namespace py = pybind11;

#define CHECK_KMER_LENGTH(kmer, k, type) ({\
    if( strlen(kmer) != k )\
    {\
        char buffer[1000];\
        sprintf( buffer, "kmer %s of length %d does not match the %s length of %d",\
            kmer, strlen(kmer), type, k );\
        throw std::length_error(std::string(buffer));\
    }\
})

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

inline bool kcontainer_contains( Kcontainer* kd, char* kmer )
{
    uint8_t* bseq = ( uint8_t* ) calloc( kd->k, sizeof( uint8_t ) );
    serialize_kmer( kmer, kd->k, bseq );
    bool res = vertex_contains( &( kd->v ), bseq, kd->k, 0 );
    free( bseq );
    return res;
}

#if KDICT
inline void kcontainer_insert( Kcontainer* kd, char* kmer, py::handle* obj )
#elif KSET
inline void kcontainer_insert( Kcontainer* kd, char* kmer )
#elif KCOUNTER
inline void kcontainer_insert( Kcontainer* kd, char* kmer, int count )
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
inline void kcontainer_add_seq(Kcontainer* kd, char* seq) {
    uint8_t* bseq = (uint8_t*) calloc(kd->k, sizeof(uint8_t));
    uint char_pos = kd->k - 1;
    uint bk = calc_bk(kd->k);
    int count = 0;
    uint8_t holder;
    uint8_t last_index = (kd->k - 1) % 4;

    // serialize the first kmer
    serialize_kmer(seq, kd->k, bseq);
    //std::cout << deserialize_kmer(kd->k, bk, bseq) << std::endl;

    while(true) {
#if KSET
        vertex_insert(&(kd->v), bseq, kd->k, 0);
#elif KCOUNTER
        // get current count of bseq
        int count = vertex_get_counter(&(kd->v), bseq, kd->k, 0);
        vertex_insert(&(kd->v), bseq, kd->k, 0, ++count);
#endif
        //std::cout << seq[char_pos] << std::endl;

        // increment char_pos
        char_pos++;
        /*count++;

        if(count > 100) {
            break;
        }*/
        if(seq[char_pos] == '\0') {
            break;
        }


        // shift all the bits over
        bseq[0] >>= 2;
        for(int i = 1; i < bk; i++) {
            holder = bseq[i] & 0x3;
            holder <<= 6;
            bseq[i - 1] |= holder;
            bseq[i] >>= 2;
        }

        switch( seq[char_pos] )
        {
            case 'a': break;
            case 'A': break;
            case 'c': bseq[bk - 1] |= MASK_INSERT[0][last_index]; break;
            case 'C': bseq[bk - 1] |= MASK_INSERT[0][last_index]; break;
            case 'g': bseq[bk - 1] |= MASK_INSERT[1][last_index]; break;
            case 'G': bseq[bk - 1] |= MASK_INSERT[1][last_index]; break;
            case 't': bseq[bk - 1] |= MASK_INSERT[2][last_index]; break;
            case 'T': bseq[bk - 1] |= MASK_INSERT[2][last_index]; break;
            default:
                std::cout << seq[char_pos] << std::endl;
                throw std::runtime_error( "Could not serialize kmer." );
        }
        //std::cout << (char_pos - kd->k) << "\t" << deserialize_kmer(kd->k, bk, bseq) << std::endl;
    }

}
#endif

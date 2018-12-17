#pragma once

#include "BloomFilter.h"
#include "helper.h"
#include "globals.h"
#include <pybind11/pybind11.h>
#include <jemalloc/jemalloc.h>

namespace py = pybind11;



struct __attribute__ ((__packed__)) CS;

#if KDICT
struct __attribute__ ((__packed__)) dVertex;
using Vertex = dVertex;
#elif KSET
struct __attribute__ ((__packed__)) sVertex;
using Vertex = sVertex;
#elif KCOUNTER
struct __attribute__ ((__packed__)) cVertex;
using Vertex = cVertex;
#endif


struct __attribute__ ((__packed__)) CC {
    BloomFilter bf;
    uint16_t suffix_size;
    uint16_t size;
    CS* child_suffixes;
    uint32_t pref[ 8 ];
};

void init_cc( CC* cc, int suffix_size );
Vertex* cc_insert( CC* cc, int k, int depth, uint8_t* sfpx );
bool cc_may_contain( CC* cc, uint8_t* bseq );
void free_cc( CC* cc );
int cc_contains_prefix( CC*cc, uint8_t* sfpx );
int rank( CC* cc, int clust_pos );
int hamming_weight( CC* cc, int index );
Vertex* insert_cluster( CC* cc, int idx, uint8_t prefix, Vertex* v, uint8_t suffix );
uint8_t* get_suffix( CC* cc, int idx );
int index_of( CC* cc, uint8_t* sfpx );
Vertex* get_child_of( CC* cc, uint8_t* sfpx, int idx );
int pref_index_from_hamming_weight( CC* cc, int clust_num );
int clust_num_from_rank( CC* cc, int clust_pos );



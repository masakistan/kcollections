#pragma once

#include "BloomFilter.h"
#include "helper.h"
#include "globals.h"
#include <pybind11/pybind11.h>
#include <jemalloc/jemalloc.h>

namespace py = pybind11;
typedef struct __attribute__ ((__packed__)) Vertex;
typedef struct __attribute__ ((__packed__)) CS;

typedef struct {
    BloomFilter bf;
    uint16_t suffix_size;
    uint16_t size;
    CS* child_suffixes;
    uint32_t pref[ 8 ];
} __attribute__ ((__packed__)) CC;

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



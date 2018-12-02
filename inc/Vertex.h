#pragma once

#include <pybind11/pybind11.h>
#include "UContainer.h"
#include "globals.h"
#include <jemalloc/jemalloc.h>
#include "uint256_t.h"

namespace py = pybind11;

#if KDICT
struct dVertex{
    dVertex* vs;
#elif KSET
struct sVertex{
    sVertex* vs;
#endif
    uint256_t pref_pres;
    UC uc;
    uint16_t vs_size;
    bool start;
};

#if KDICT
using Vertex = dVertex;
#elif KSET
using Vertex = sVertex;
#endif

struct CS {
    Vertex v;
    uint8_t suffix;
};

void vertex_remove( Vertex* v, uint8_t* bseq, int k, int depth );
uint64_t vertex_size( Vertex* v );
void init_vertex( Vertex* v );
void free_vertex( Vertex* v );
bool vertex_contains( Vertex* v, uint8_t* bseq, int k, int depth );
void burst_uc( Vertex* v, int k, int depth );

#if KDICT
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth, py::handle* obj );
#elif KSET
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth );
#endif

py::handle* vertex_get( Vertex* v, uint8_t* bseq, int k, int depth );



#pragma once

#include <pybind11/pybind11.h>
#include "CContainer.h"
#include "UContainer.h"
#include "globals.h"
#include <jemalloc/jemalloc.h>

namespace py = pybind11;

struct Vertex{
    CC* cc;
    UC uc;
    uint16_t cc_size;
    bool start;
};

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
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth, py::object* obj );
#elif KSET
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth );
#endif

py::object* vertex_get( Vertex* v, uint8_t* bseq, int k, int depth );



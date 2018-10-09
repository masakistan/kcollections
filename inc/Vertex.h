#pragma once

#include "CContainer.h"
#include "UContainer.h"
#include "globals.h"
#include <jemalloc/jemalloc.h>

/*typedef struct {
    CC* cc;
    UC uc;
} __attribute__ ((__packed__)) Vertex;*/

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

void init_vertex( Vertex* v );
void free_vertex( Vertex* v );
bool vertex_contains( Vertex* v, uint8_t* bseq, int k, int depth );
void burst_uc( Vertex* v, int k, int depth );
void vertex_insert( Vertex* v, uint8_t* bseq, int k, int depth );



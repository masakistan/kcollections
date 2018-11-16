#pragma once

#include <utility>
#include <stdio.h>
#include <iostream>
#include "globals.h"
#include "helper.h"
#include <pybind11/pybind11.h>
#include <jemalloc/jemalloc.h>

namespace py = pybind11;

typedef struct {
    uint8_t* suffixes;
#if KDICT
    py::handle* objs;
#elif KCOUNTER
    int* counts;
#endif
    uint16_t size;
} UC;

void print( UC* uc, int k, int depth );

#if KDICT
void uc_insert( UC* uc, uint8_t* bseq, int k, int depth, int idx, py::handle* obj );
#elif KSET
void uc_insert( UC* uc, uint8_t* bseq, int k, int depth, int idx );
#elif KCOUNTER
void uc_insert( UC* uc, uint8_t* bseq, int k, int depth, int idx, int count );
#endif

void free_uc( UC* uc );
int uc_contains( UC* uc, int k, int depth, uint8_t* bseq );
std::pair< bool, int > uc_find( UC* uc, int k, int depth, uint8_t* bseq );
void init_uc( UC* uc );
void uc_remove( UC* uc, int bk, int idx );


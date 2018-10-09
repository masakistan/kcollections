#pragma once

#include <stdio.h>
#include <iostream>
#include "globals.h"
#include "helper.h"
#include <jemalloc/jemalloc.h>

typedef struct {
    uint8_t* suffixes;
    uint16_t size;
} __attribute__ ((__packed__)) UC;

void print( UC* uc, int k, int depth );
void uc_insert( UC* uc, uint8_t* bseq, int k, int depthi, int idx );
void free_uc( UC* uc );
int uc_contains( UC* uc, int k, int depth, uint8_t* bseq );
void init_uc( UC* uc );



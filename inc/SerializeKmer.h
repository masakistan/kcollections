#pragma once

#include <iostream>
#include <string>
#include <sstream>

#include "Bkmer.h"

static const uint8_t MASK_INSERT[ 3 ][ 4 ] = {
        { 1, 4, 16, 64 },
        { 2, 8, 32, 128 },
        { 3, 12, 48, 192 }
    };
static const char COMP_TO_ASCII[4] = {'A', 'C', 'G', 'T'};

Bkmer* serialize_kmer( char* kmer, int k, int bk );
char* deserialize_bkmer( Bkmer* bkmer, int k, int bk );


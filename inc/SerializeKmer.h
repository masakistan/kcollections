#pragma once

#include <string>
#include <sstream>

static const uint8_t MASK_INSERT[ 3 ][ 4 ] = {
        { 1, 4, 16, 64 },
        { 2, 8, 32, 128 },
        { 3, 12, 48, 192 }
    };

uint8_t* serializeKmer( char* kmer, int k, int bk );



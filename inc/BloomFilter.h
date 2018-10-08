#pragma once

#include <vector>
#include <array>
#include "MurmurHash3.h"
#include "globals.h"
#include "helper.h"


inline uint64_t nthHash( uint8_t n,
                         uint64_t hashA,
                         uint64_t hashB,
                         uint64_t filterSize )
{
    return (hashA + n * hashB) % filterSize;
}

typedef struct {
    uint32_t m_bits[ 8 ];
    int m_nHashes;
} __attribute__ ((__packed__)) BloomFilter;

inline void init_bf( BloomFilter* bf )
{
    bf->m_nHashes = NHASHES;
    memset( &bf->m_bits, 0, 8 );
}

static std::array< uint64_t, 2 > hash( const uint8_t* data, const std::size_t len )
{
    std::array< uint64_t, 2 > hashValue;
    MurmurHash3_x64_128( data, len, 0, hashValue.data() );
    return hashValue;
}

inline void add_to_bloom_filter( BloomFilter* bf, uint8_t* data, int size )
{
    auto hashValues = hash( data, size );

    for ( int n = 0; n < bf->m_nHashes; n++ )
    {
        // 32 * 8 = 256, number of bits in bloom filter
        setbit( bf->m_bits, nthHash( n, hashValues[ 0 ], hashValues[ 1 ], 32 * 8 ) );
    }
}

inline bool bf_may_contain( BloomFilter* bf, uint8_t* data, int size )
{
    auto hashValues = hash( data, size );
    for ( int n = 0; n < bf->m_nHashes; n++ )
    {
        if ( !testbit( bf->m_bits, nthHash( n, hashValues[ 0 ], hashValues[ 1 ], 32 * 8 ) ) )
        {
            return false;
        }
    }

    return true;
}

inline void free_bf( BloomFilter* bf )
{
}



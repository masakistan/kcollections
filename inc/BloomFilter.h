#pragma once

#include <vector>
#include <array>
#include "MurmurHash3.h"
#include "Bkmer.h"


inline uint64_t nthHash( uint8_t n,
                         uint64_t hashA,
                         uint64_t hashB,
                         uint64_t filterSize )
{
    return (hashA + n * hashB) % filterSize;
}


class __attribute__ ((__packed__)) BloomFilter
{
    private:
        std::vector<bool>* m_bits;
        uint8_t m_nHashes;

        std::array< uint64_t, 2 > hash( const uint8_t* data, std::size_t len ) const;

    public:
        BloomFilter( uint64_t size, uint8_t nHashes);
        BloomFilter( const BloomFilter& bf );
        ~BloomFilter();
        void add( const Bkmer* data );
        bool may_contain( const uint8_t* data, std::size_t len ) const;
        std::vector<bool>* get_bits() const { return m_bits; }
};



#pragma once

#include <vector>
#include <array>
#include "MurmurHash3.h"
#include "Globals.h"

class BloomFilter
{
    private:
        std::vector<bool>* m_bits;
        uint8_t m_nHashes;
        std::array< uint64_t, 2 > hash( const uint8_t* data, std::size_t len );

    public:
        BloomFilter( uint64_t size, uint8_t nHashes);
        ~BloomFilter();
        bool may_contain();
};



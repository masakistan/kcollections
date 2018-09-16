#include "BloomFilter.h"

BloomFilter::BloomFilter( uint64_t size, uint8_t nHashes ) : m_nHashes( nHashes )
{
    m_bits = new std::vector< bool >( size );
}

BloomFilter::~BloomFilter()
{
    delete m_bits;
}

bool BloomFilter::may_contain()
{
    return false;
}

std::array< uint64_t, 2 > BloomFilter::hash( const uint8_t* data, std::size_t len )
{
  std::array< uint64_t, 2 > hashValue;
  MurmurHash3_x64_128( data, len, 0, hashValue.data() );

  return hashValue;
}



#include "BloomFilter.h"

BloomFilter::BloomFilter( uint64_t size, uint8_t nHashes ) : m_nHashes( nHashes )
{
    m_bits = new std::vector< bool >( size );
}

BloomFilter::~BloomFilter()
{
    delete m_bits;
}

std::array< uint64_t, 2 > BloomFilter::hash( const uint8_t* data, std::size_t len ) const
{
  std::array< uint64_t, 2 > hashValue;
  MurmurHash3_x64_128( data, len, 0, hashValue.data() );

  return hashValue;
}

void BloomFilter::add(const uint8_t *data, std::size_t len)
{
  auto hashValues = hash( data, len );

  for ( int n = 0; n < m_nHashes; n++ )
  {
      ( *m_bits )[ nthHash( n, hashValues[ 0 ], hashValues[ 1 ], m_bits->size() ) ] = true;
  }
}

bool BloomFilter::may_contain( const uint8_t *data, std::size_t len) const
{
  auto hashValues = hash( data, len );

  for ( int n = 0; n < m_nHashes; n++ )
  {
      if ( !( *m_bits )[ nthHash( n, hashValues[ 0 ], hashValues[ 1 ], m_bits->size() ) ] )
      {
          return false;
      }
  }

  return true;
}



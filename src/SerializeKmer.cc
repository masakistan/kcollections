#include "SerializeKmer.h"

uint8_t* serializeKmer( char* kmer, int k, int bk )
{
    // allocate the binary kmer (bkmer)
    uint8_t* bkmer = ( uint8_t* ) calloc( bk, sizeof( uint8_t ) );
;
    for( int pos; pos < k; pos++ )
    {
        switch( kmer[ pos ] )
        {
            case 'a': break;
            case 'A': break;
            case 'c': bkmer[ pos / 4 ] |= MASK_INSERT[ 0 ][ pos % 4 ]; break;
            case 'C': bkmer[ pos / 4 ] |= MASK_INSERT[ 0 ][ pos % 4 ]; break;
            case 'g': bkmer[ pos / 4 ] |= MASK_INSERT[ 1 ][ pos % 4 ]; break;
            case 'G': bkmer[ pos / 4 ] |= MASK_INSERT[ 1 ][ pos % 4 ]; break;
            case 't': bkmer[ pos / 4 ] |= MASK_INSERT[ 2 ][ pos % 4 ]; break;
            case 'T': bkmer[ pos / 4 ] |= MASK_INSERT[ 2 ][ pos % 4 ]; break;
            default: std::ostringstream errMsg;
                     errMsg << "Runtime Error: bad symbol \"" << kmer[ pos ] << "\" when serializing kmer";
                     throw std::runtime_error( errMsg );
        }
    }

    return bkmer;
}



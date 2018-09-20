#include "SerializeKmer.h"

Bkmer* serialize_kmer( char* kmer, int k, int bk )
{
    std::cout << "serializing " << kmer << std::endl;
    // allocate the binary kmer (bkmer)
    //uint8_t* bkmer = ( uint8_t* ) calloc( bk, sizeof( uint8_t ) );
    Bkmer* bkmer = new Bkmer( k, bk );
    uint8_t* bseq = bkmer->get_bseq();

    for( int pos = 0; pos < k; pos++ )
    {
        switch( kmer[ pos ] )
        {
            case 'a': break;
            case 'A': break;
            case 'c': bseq[ pos / 4 ] |= MASK_INSERT[ 0 ][ pos % 4 ]; break;
            case 'C': bseq[ pos / 4 ] |= MASK_INSERT[ 0 ][ pos % 4 ]; break;
            case 'g': bseq[ pos / 4 ] |= MASK_INSERT[ 1 ][ pos % 4 ]; break;
            case 'G': bseq[ pos / 4 ] |= MASK_INSERT[ 1 ][ pos % 4 ]; break;
            case 't': bseq[ pos / 4 ] |= MASK_INSERT[ 2 ][ pos % 4 ]; break;
            case 'T': bseq[ pos / 4 ] |= MASK_INSERT[ 2 ][ pos % 4 ]; break;
            default: std::ostringstream errMsg;
                     errMsg << "Runtime Error: bad symbol \"" << kmer[ pos ] << "\" when serializing kmer";
                     throw std::runtime_error( errMsg.str() );
        }
    }

    return bkmer;
}

char* deserialize_bkmer( Bkmer* bkmer, int k, int bk )
{
    int bases_processed = 0, j, pos, bases_to_process;
    char* kmer = ( char* ) malloc( sizeof( char ) * k );
    uint8_t tbkmer;
    
    for( int i = 0; i < bk; i ++ )
    {
        if( i == bk - 1 )
        {
            bases_to_process = k % 4;
        }
        else
        {
            bases_to_process = 4;
        }

        tbkmer = bkmer->get_bseq()[ i ];
        for( j = 0; j < bases_to_process; j ++ )
        {
            pos = ( i * 4 ) + j;
            kmer[ pos ] = COMP_TO_ASCII[ tbkmer & 0x3 ];
            tbkmer >>= 2;
        }
    }

    return kmer;
}



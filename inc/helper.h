#pragma once

#include <stdexcept>
#include <utility>
#include <stdlib.h>
#include <sstream>
#include <stdint.h>
#include <iostream>
#include <cstring>

#define CHECK_KMER_LENGTH(kmer, k, type) ({     \
    if( strlen(kmer) != k )\
    {\
        char buffer[1000];\
        sprintf( buffer, "kmer %s of length %d does not match the %s length of %d",\
            kmer, strlen(kmer), type, k );\
        throw std::length_error(std::string(buffer));\
    }\
})


static int testbit( uint32_t A[],  unsigned int k )
   {
      return ( (A[k/32] & (1 << (k%32) )) != 0 ) ;
   }

static void  clearbit( uint32_t A[],  unsigned int k )                
   {
      A[k/32] &= ~(1 << (k%32));
   }

static void  setbit( uint32_t A[],  unsigned int k )
   {
      A[k/32] |= 1 << (k%32);  // Set the bit at the k-th position in A[i]
   }

static int next_set_bit( uint32_t* array, int pos, int len )
{
    for( unsigned int i = pos; i < len; i++ )
    {
        if( testbit( array, i ) )
        {
            return i;
        }
    }
    return -1;
}

inline int calc_bk( int k )
{
    int bk = k / 4;
    if( k % 4 > 0 )
    {
        bk++;
    }

    return bk;
}

static const uint8_t MASK_INSERT[ 4 ][ 4 ] = {
        {0, 0, 0, 0, },
        //{ 1, 4, 16, 64 },
        {64, 16, 4, 1},
        //{ 2, 8, 32, 128 },
        {128, 32, 8, 2},
        //{ 3, 12, 48, 192 }
        {192, 48, 12, 3}
    };

static const char COMP_TO_ASCII[4] = {'A', 'C', 'G', 'T'};


static void serialize_position(uint32_t kmerPos, int arrPos, int bitPos, uint8_t* bseq, const char* kmer) {
    switch( kmer[kmerPos] )
    {
        case 'a': break;
        case 'A': break;
        case 'c': bseq[ arrPos ] |= MASK_INSERT[ 1 ][ bitPos ]; break;
        case 'C': bseq[ arrPos ] |= MASK_INSERT[ 1 ][ bitPos ]; break;
        case 'g': bseq[ arrPos ] |= MASK_INSERT[ 2 ][ bitPos ]; break;
        case 'G': bseq[ arrPos ] |= MASK_INSERT[ 2 ][ bitPos ]; break;
        case 't': bseq[ arrPos ] |= MASK_INSERT[ 3 ][ bitPos ]; break;
        case 'T': bseq[ arrPos ] |= MASK_INSERT[ 3 ][ bitPos ]; break;
        default:
            bseq[ arrPos ] |= MASK_INSERT[ rand() % 4 ][ bitPos ]; break;
            //throw std::runtime_error( "Could not serialize kmer." );
    }
}

static void serialize_kmer( const char* kmer, int k, uint8_t* bseq )
{
    for( int pos = 0; pos < k; pos++ )
    {
        serialize_position(pos, pos / 4, pos % 4, bseq, kmer);
    }
}

static char* deserialize_kmer( int k, int bk, uint8_t* bseq )
{
    int bases_processed = 0, j, pos, bases_in_byte, bases_to_process = k;
    char* kmer = ( char* ) malloc( sizeof( char ) * ( k + 1 ) );
    uint8_t tbkmer;
    
    for( int i = 0; i < bk; i ++ )
    {
        tbkmer = bseq[ i ];
        if( bases_to_process >= 4 )
        {
            bases_in_byte = 4;
        }
        else
        {
            bases_in_byte = bases_to_process;
            tbkmer >>= 8 - (bases_in_byte * 2);
        }

        for( j = 0; j < bases_in_byte; j++ )
        {
            pos = ( i * 4 ) + (bases_in_byte - j - 1);
            kmer[ pos ] = COMP_TO_ASCII[ tbkmer & 0x3 ];
            tbkmer >>= 2;
        }
        bases_to_process -= bases_in_byte;
    }
    kmer[ k ] = '\0';

    return kmer;
}

static int compare_seqs(uint8_t* seq1, uint8_t* seq2, int len) {
    //std::cout << "compare\n" << deserialize_kmer(len * 4, len, seq1) << "\n";
    //std::cout << deserialize_kmer(len*4, len, seq2) << std::endl;
    for(int i = 0; i < len; i++) {
        if((unsigned) seq1[i] < (unsigned) seq2[i]) {
            return -1;
        } else if((unsigned) seq1[i] > (unsigned) seq2[i]) {
            return 1;
        }
    }
    return 0;
}

static std::pair< bool, int > binary_search( uint8_t* suffixes, int max, int len, uint8_t* bseq )
{
    if( max == 0 )
    {
        return std::make_pair( false, 0 );
    }

    int min = 0, mid, idx, cmp;
    while( min < max )
    {
        mid = min + ( max - min ) / 2;
        idx = mid * len;

        //cmp = std::memcmp( bseq, &suffixes[ idx ], len );
        cmp = compare_seqs(bseq, &suffixes[idx], len);
        if( cmp < 0 )
        {
            max = mid;
        }
        else if( cmp == 0 )
        {
            return std::make_pair( true, mid );
        }
        else
        {
            min = mid + 1;;
        }
    }

    return std::make_pair( false, min );
}



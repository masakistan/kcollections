#pragma once

#include <stdlib.h>
#include <sstream>
#include <stdint.h>
#include <iostream>
#include <cstring>

//#define setbit(A,k)     ( A[(k/32)] |= (1 << (k%32)) )
//#define clearbit(A,k)   ( A[(k/32)] &= ~(1 << (k%32)) )
//#define testbit(A,k)    ( A[(k/32)] & (1 << (k%32)) )

static int testbit( uint32_t A[],  int k )
   {
      return ( (A[k/32] & (1 << (k%32) )) != 0 ) ;
   }

static void  clearbit( uint32_t A[],  int k )                
   {
      A[k/32] &= ~(1 << (k%32));
   }

static void  setbit( uint32_t A[],  int k )
   {
      A[k/32] |= 1 << (k%32);  // Set the bit at the k-th position in A[i]
   }

static int next_set_bit( uint32_t* array, int pos, int len )
{
    //std::cout << "nsb start at: " << pos << std::endl;
    //std::cout << "\n\n" << std::endl;
    for( int i = pos; i < len; i++ )
    {
        //std::cout << "\t\t\tnsb check: " << i << "\t" << testbit( array, i ) << std::endl;
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

static const uint8_t MASK_INSERT[ 3 ][ 4 ] = {
        /*{ 64, 16, 4, 1 },
        { 128, 32, 8, 2 },
        { 192, 48, 12, 3 }*/
        { 1, 4, 16, 64 },
        { 2, 8, 32, 128 },
        { 3, 12, 48, 192 }
    };

static const char COMP_TO_ASCII[4] = {'A', 'C', 'G', 'T'};

static void serialize_kmer( char* kmer, int k, uint8_t* bseq )
{
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
}

static char* deserialize_kmer( int k, int bk, uint8_t* bseq )
{
    int bases_processed = 0, j, pos, bases_in_byte, bases_to_process = k;
    char* kmer = ( char* ) malloc( sizeof( char ) * ( k + 1 ) );
    uint8_t tbkmer;
    
    for( int i = 0; i < bk; i ++ )
    {
        if( bases_to_process > 4 )
        {
            bases_in_byte = 4;
        }
        else
        {
            bases_in_byte = bases_to_process;
        }

        tbkmer = bseq[ i ];
        for( j = 0; j < bases_in_byte; j++ )
        {
            pos = ( i * 4 ) + j;
            kmer[ pos ] = COMP_TO_ASCII[ tbkmer & 0x3 ];
            tbkmer >>= 2;
        }
        bases_to_process -= bases_in_byte;
    }
    kmer[ k ] = '\0';

    return kmer;
}

struct compare_bseq
{
    inline bool operator() ( const uint8_t* left, const uint8_t* right )
    {
        /*std::cout << "comparing " << left->get_seq() << "\t" << right->get_seq() << std::endl;
        if( left < right )
            std::cout << "\tleft < right" << std::endl;
        else
            std::cout << "\tleft >= right" << std::endl;*/
        return *left < *right;
    }
};

static int binary_search( uint8_t* suffixes, int max, int len, uint8_t* bseq )
{
    if( max == 0 )
    {
        return 0;
    }

    int min = 0, mid, idx, cmp;
    while( min < max )
    {
        mid = min + ( max - min ) / 2;
        idx = mid * len;

        cmp = std::memcmp( bseq, &suffixes[ idx ], len );
        if( cmp < 0 )
        {
            max = mid;
        }
        else
        {
            min = mid + 1;;
        }
    }

    return min;
}

static bool binary_search_contains( uint8_t* suffixes, int max, int len, uint8_t* bseq )
{
    if( max == 0 )
    {
        return 0;
    }

    int min = 0, mid, idx, cmp;
    while( min < max )
    {
        mid = min + ( max - min ) / 2;
        idx = mid * len;

        cmp = std::memcmp( bseq, &suffixes[ idx ], len );
        if( cmp < 0 )
        {
            max = mid;
        }
        else if( cmp == 0 )
        {
            return true;
        }
        else
        {
            min = mid + 1;;
        }
    }

    return false;
}



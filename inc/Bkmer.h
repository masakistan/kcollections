#pragma once

#include <stdlib.h>
#include <stdint.h>
#include "Globals.h"

struct Bkmer
{
    uint8_t* bseq;
    Bkmer()
    {
        bseq = ( uint8_t* ) calloc( BK, sizeof( uint8_t ) );
    }

    Bkmer( const Bkmer& other )
    {
        bseq = ( uint8_t* ) calloc( BK, sizeof( uint8_t ) );
        for( int i = 0; i < BK; i++ )
        {
            bseq[ i ] = other.bseq[ i ];
        }
    }

    ~Bkmer()
    {
        free( bseq );
    }

    bool operator<( Bkmer other ) const
    {
        for( int i = 0; i < BK; i++ )
        {
            if(  unsigned( bseq[ i ] ) < unsigned( other.bseq[ i ] ) )
            {
                return true;
            }
            else if( unsigned( bseq[ i ] ) > unsigned( other.bseq[ i ] ) )
            {
                return false;
            }
            else
            {
                // continue search
            }
        }
        return false;
    }
};



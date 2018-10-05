#pragma once

#include <bitset>

inline int calc_bk( int k )
{
    int bk = k / 4;
    if( k % 4 > 0 )
    {
        bk++;
    }

    return bk;
}

inline void shift( std::bitset< 256 >& b, int idx )
{
    for( int i = 255; i > idx; i-- )
    {
        b[ i ] = b[ i - 1 ];
    }
}



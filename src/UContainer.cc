#include "UContainer.h"

void init_uc( UC* uc )
{
    uc->suffixes = NULL;
    uc->size = 0;
}

void print( UC* uc, int k, int depth )
{
    int len = calc_bk( k );
    int idx;
    for( int i = 0; i < uc->size; i++ )
    {
        idx = i * len;
        char* dseq = deserialize_kmer( k, len, &uc->suffixes[ idx ] );
        std::cout << "kmer: " << dseq << std::endl;
        free( dseq );
    }
}

void uc_insert( UC* uc, uint8_t* bseq, int k, int depth, int idx )
{
    int len = calc_bk( k );
    if( uc->suffixes == NULL )
    {
        uc->suffixes = ( uint8_t* ) calloc( len , sizeof( uint8_t ) );
    }
    else
    {
        uc->suffixes = ( uint8_t* ) realloc(
                uc->suffixes,
                len * ( uc->size + 1 ) * sizeof( uint8_t )
                );
    }

    if( uc->size < CAPACITY )
    {
        //int idx = binary_search( uc->suffixes, uc->size, len, bseq );

        int bytes_to_move = ( uc->size - idx ) * len;
        idx = idx * len;
        if( bytes_to_move > 0 )
        {
            std::memmove(
                    &uc->suffixes[ idx + len ],
                    &uc->suffixes[ idx ],
                    bytes_to_move
                    );
        }

        std::memcpy( &uc->suffixes[ idx ], bseq, len );
        uc->size++;
    }
}

void free_uc( UC* uc )
{
    if( uc->suffixes != NULL )
    {
        free( uc->suffixes );
    }
}

int uc_contains( UC* uc, int k, int depth, uint8_t* bseq )
{
    if( uc->suffixes == NULL )
    {
        return 0;
    }

    int idx = binary_search_contains( uc->suffixes, uc->size, calc_bk( k ), bseq );
    if( idx > -1 )
    {
        return idx;
    }
    else
    {
        return -1;
    }
}



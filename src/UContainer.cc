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

void uc_insert( UC* uc, uint8_t* bseq, int k, int depth )
{
    int len = calc_bk( k );
    //std::cout << "\tuc insert: " << deserialize_kmer( k, len, bseq ) << std::endl;
    if( uc->suffixes == NULL )
    {
        uc->suffixes = ( uint8_t* ) calloc( len * CAPACITY, sizeof( uint8_t ) );
    }

    if( uc->size < CAPACITY )
    {
        int idx = binary_search( uc->suffixes, uc->size, len, bseq );

        //std::cout << "inserting at idx: " << idx << std::endl;
        int bytes_to_move = ( uc->size - idx ) * len;
        idx = idx * len;
        if( bytes_to_move > 0 )
        {
            //std::cout << "\tbytes to move: " << bytes_to_move << std::endl;
            std::memmove( &uc->suffixes[ idx + len ], &uc->suffixes[ idx ], bytes_to_move );
        }

        std::memmove( &uc->suffixes[ idx ], bseq, len );
        uc->size++;
    }

    /*std::cout << "size: " << uc->size << std::endl;
    std::cout << "***************************************************************" << std::endl;
    print( uc, k, depth );
    std::cout << "***************************************************************" << std::endl;*/
}

void free_uc( UC* uc )
{
    if( uc->suffixes != NULL )
    {
        free( uc->suffixes );
    }
}

bool uc_contains( UC* uc, int k, int depth, uint8_t* bseq )
{
    //std::cout << "search index: " << binary_search_contains( uc->suffixes, uc->size, len, bseq ) << std::endl;
    if( uc->suffixes == NULL )
    {
        return false;
    }

    if( binary_search_contains( uc->suffixes, uc->size, calc_bk( k ), bseq ) )
    {
        return true;
    }
    else
    {
        return false;
    }
}



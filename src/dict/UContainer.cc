#include "dict/UContainer.h"

void init_uc( UC* uc )
{
    uc->suffixes = NULL;
    uc->objs = NULL;
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

void uc_insert( UC* uc, uint8_t* bseq, int k, int depth, int idx, py::object* obj )
{
    int len = calc_bk( k );
    if( uc->suffixes == NULL )
    {
        uc->suffixes = ( uint8_t* ) calloc( len , sizeof( uint8_t ) );
        uc->objs = ( py::object* ) calloc( len , sizeof( py::object ) );
    }
    else
    {
        uc->suffixes = ( uint8_t* ) realloc(
                uc->suffixes,
                len * ( uc->size + 1 ) * sizeof( uint8_t )
                );
        uc->objs = ( py::object* ) realloc(
                uc->suffixes,
                len * ( uc->size + 1 ) * sizeof( py::object )
                );
    }

    if( uc->size < CAPACITY )
    {
        //int idx = binary_search( uc->suffixes, uc->size, len, bseq );

        int objs_bytes_to_move = ( uc->size - idx ) * sizeof( py::object );
        int suffix_bytes_to_move = ( uc->size - idx ) * len;
        int suffix_idx = idx * len;
        if( suffix_bytes_to_move > 0 )
        {
            std::memmove(
                    &uc->suffixes[ suffix_idx + len ],
                    &uc->suffixes[ suffix_idx ],
                    suffix_bytes_to_move
                    );
            std::memmove(
                    &uc->objs[ idx + 1 ],
                    &uc->objs[ idx ],
                    objs_bytes_to_move
                    );
        }

        std::memcpy( &uc->objs[ idx ], obj, sizeof( py::object ) );
        std::memcpy( &uc->suffixes[ suffix_idx ], bseq, len );
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

void uc_remove( UC* uc, int bk, int idx )
{
    int suffix_idx = idx * bk;
    int bytes_to_move = ( uc->size - ( idx + 1 ) ) * bk;
    std::memmove(
            &uc->suffixes[ idx ],
            &uc->suffixes[ idx + bk ],
            bytes_to_move
            );
    uc->size--;
}

int uc_find( UC* uc, int k, int depth, uint8_t* bseq )
{
    if( uc->suffixes == NULL )
    {
        return uc->size;
    }

    int idx = binary_search_contains( uc->suffixes, uc->size, calc_bk( k ), bseq );
    if( idx > -1 )
    {
        return uc->size;
    }
    else
    {
        idx = ( idx + 1 ) * -1;
        return idx;
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


